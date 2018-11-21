
"""
call_peaks uses the uncontinuous score from a bigWig file to estimate peaks
@author: Anastasiia Petrova
@contact: anastasiia.petrova(at)mpi-bn.mpg.de
"""

import argparse
import sys
import os
import re
import time
#import multiprocessing
import logging
#import subprocess
#from Bio import SeqIO
#import Bio.SeqIO.FastaIO as bio
import numpy as np
from collections import defaultdict
#from scipy import stats
import pyBigWig
#from statsmodels.sandbox.stats.multicomp import multipletests #for bonfferoni
#import matplotlib.pyplot as plt
#import random
import textwrap

logger = logging.getLogger('call_peaks')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

def parse_args():
	
	parser = argparse.ArgumentParser(prog = '', description = textwrap.dedent('''                           

		This script produces a file with peaks in .bed format from the file with scores in .bigWig format.
		'''), epilog='That is what you need to make this script work for you. Enjoy it')
	
	required_arguments = parser.add_argument_group('required arguments')
	required_arguments.add_argument('--bigwig', help='a bigWig-file with scores', required=True)
	required_arguments.add_argument('--bed', help='provide a file with peaks in .bed format', required=True)

	#all other arguments are optional
	#parser.add_argument('--output_directory',  default='output', const='output', nargs='?', help='output directory, by default ./output/')
	parser.add_argument('--output_file', default='call_peaks_output.bed', type=str, help='enter a name for the output file, by default ./call_peaks_output.bed')
	parser.add_argument('--window_length', default=200, type=int, help='enter the length for a window, by defauld 200 bp')
	parser.add_argument('--step', default=100, type=int, help='enter a step to move the window, by default 100 bp')
	parser.add_argument('--percentage',  default=10, type=int, help='enter a percentage to be added to background while searching for footprints, by default 10%')
	parser.add_argument('--silent', action='store_true', help='while working with data write the information only into ./call_peaks_log.txt')
	args = parser.parse_args()

	return args

def check_directory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
		#logger.info('a new directory ' + directory + ' was created')
		print('a new directory ' + directory + ' was created')

#if there are chars that are not allowed, they will be replaced with '_', to the possibly invalid names there will be added '_' at the beginning of the name
def check_name(name_to_test):
	badchars= re.compile(r'[^A-Za-z0-9_. ]+|^\.|\.$|^ | $|^$')
	badnames= re.compile(r'(aux|com[1-9]|con|lpt[1-9]|prn)(\.|$)')

	#replace all the chars that are not allowed with '_'
	name = badchars.sub('_', name_to_test)

	#check for the reserved by the os names
	if badnames.match(name):
		name = '_' + name
	return name

def remove_file(file):
	if os.path.isfile(file):
		os.remove(file)

def make_name_from_path(full_path, output_directory, ending):
	return os.path.join(output_directory, get_name_from_path(full_path) + ending)
	
def get_name_from_path(full_path):
	return os.path.splitext(os.path.basename(full_path))[0]

def check_existing_input_files(args):

	#check if the bigWig file exists
	if not os.path.isfile(args.bigwig):
		print('there is no such bigWig file ' + args.bigwig + ', the exit is forced')
		sys.exit()
	#check if the bed file exists
	if not os.path.isfile(args.bed):
		#logger.info('there is no such bed file ' + args.bed_file + ', the exit is forced')
		print('there is no such bed file ' + args.bed + ', the exit is forced')
		sys.exit()

def make_bed_dictionary(bed_file):

	logger.info('reading of the bed file')

	bed_dictionary = {}

	with open(bed_file) as read_bed_file:
		for bed_line in read_bed_file:
			bed_line_array = re.split(r'\t', bed_line.rstrip('\n'))
			if bed_line_array[1].isdigit() and bed_line_array[2].isdigit() and int(bed_line_array[1]) <= int(bed_line_array[2]): #in the real bedfile the second column is a start position, and the third column is an end position, so we are checking if these are integers and if the start position is smaller than the end one
				key = bed_line_array[0] + ":" + bed_line_array[1] + "-" + bed_line_array[2]
				value = []
				for i in range(3, len(bed_line_array)):
					value.append(bed_line_array[i]) 

				bed_dictionary[key] = value
			else: #this is not a bed file, force exit
				logger.info('please make sure the input bed file has a right format, the problem occured on the line ' + bed_line)
				sys.exit()

	read_bed_file.close()

	return bed_dictionary

def save_footprint(footprint_count, footprint_scores, peak_footprints, chromosom, footprint_start, footprint_end, bonus_info_from_bed):

	save_current_footprint = False

	if len(footprint_scores) > 2: #exclude small regions to not work with them

		#-------------- max score position

		first_max_pos = footprint_scores.index(max(footprint_scores))
		last_max_pos = first_max_pos #assume that there is only one pos with max score

		#find the region with the highest score
		for j in range(first_max_pos, len(footprint_scores)):
			if footprint_scores[j] < first_max_pos:
				last_max_pos = j
			else:
				last_max_pos = len(footprint_scores)

		if first_max_pos != last_max_pos:
			#find a pos in the middle of these both
			max_pos = int((last_max_pos - first_max_pos) / 2)
		else:
			max_pos = first_max_pos
		#-------------------------------------

		footprint_score = np.mean(footprint_scores)

		search_key = str(chromosom) + ":" + str(footprint_start) + "-" + str(footprint_end)

		#--------------- checking for existing and overlapping footprints
		
		if len(peak_footprints.keys()) == 0: # hey, this is the first footprint!
			save_current_footprint = True

		else: # there are already footprints in this peak, so look for overlaps
			for existing_footprint_name in peak_footprints.keys():

				old_start = peak_footprints[existing_footprint_name]['start']
				old_end = peak_footprints[existing_footprint_name]['end']
				old_score = peak_footprints[existing_footprint_name]['score']

				if footprint_start >= old_start and footprint_start <= old_end: #the start of the new footprint is between the start and end of an old footprint
					if footprint_end > old_end: #the new footprint is not completely inside the old one
						#update the information about the existing footprint
						footprint_score = (peak_footprints[existing_footprint_name]['score'] + footprint_score) / 2 #find the average of both scores

						peak_footprints[existing_footprint_name]['end'] = footprint_end
						peak_footprints[existing_footprint_name]['score'] = footprint_score
						#we can not update the max_pos as we do not have the information about scores array of the existing footprint
					#else: the new footprint is completely inside the old one, do nothing
					save_current_footprint = False
					break

				elif footprint_end >= old_start and footprint_end <= old_end: #the end of the new footprint is between the start and end of an old footprint
					if footprint_start < old_start: #the new footprint is not completely inside the old one
						#update the information about the existing footprint
						footprint_score = (peak_footprints[existing_footprint_name]['score'] + footprint_score) / 2 #find the average of both scores

						peak_footprints[existing_footprint_name]['start'] = footprint_start
						peak_footprints[existing_footprint_name]['score'] = footprint_score
					#else do nothing
					save_current_footprint = False
					break

				elif footprint_start <= old_start and footprint_end >= old_end: #the old footprint is exactly between the start and end positions of the new footprint
					#update the information about the existing footprint
					peak_footprints[existing_footprint_name] = {'chromosom': chromosom, 'start': footprint_start, 'end': footprint_end, 'score': footprint_score, 'len': len(footprint_scores), 'bonus': bonus_info_from_bed, 'max_pos': max_pos}
					save_current_footprint = False
					break

				else: #so this is a new footprint that has no overlaps with the others

					save_current_footprint = True

		if save_current_footprint == True:
			#make sure to save this footprint
			
			footprint_name = "footprint_" + str(footprint_count)

			peak_footprints[footprint_name] = peak_footprints.get(footprint_name, {})
			peak_footprints[footprint_name] = {'chromosom': chromosom, 'start': footprint_start, 'end': footprint_end, 'score': footprint_score, 'len': len(footprint_scores), 'bonus': bonus_info_from_bed, 'max_pos': max_pos}

			footprint_count = footprint_count + 1

	#else do nothing

	return footprint_count, peak_footprints

def search_in_window(peak_footprints, footprint_count, chromosom, peak_start, peak_end, scores_in_peak, window_length, bed_dictionary_entry, step, percentage):

	peak_len = len(scores_in_peak)
	parts = []
	parts_positions = []

	#if necessary divide the peak with help of a sliding window
	if peak_len <= window_length:
		window_length = peak_len
		parts.append(scores_in_peak)
		parts_positions.append(0)
	else:
		pos = 0
		while pos < (peak_len - step):
			if (pos + window_length) >= len(scores_in_peak):
				part = scores_in_peak[pos:]
			else:
				part = scores_in_peak[pos:pos + window_length]

			if len(part) != 1: #otherwise it makes no sense to look on the mean within this part and look for footprints
				parts.append(part)
				parts_positions.append(pos)

			pos = pos + step

	#look in each window and save the footprints
	for j in range(len(parts)):

		window = parts[j]

		#-------------- add some percent to the background to avoid the noice we don't want to have
		bw_peak_background = np.mean(window) #find the mean of all scores within one peak
		part = (percentage*bw_peak_background)/100 #x procent of the background
		bw_peak_background = bw_peak_background + part
		#----------------------------------------------------------
		
		check_position = parts_positions[j] #the start position not within the window, but within the peak!!!
		footprint_start = check_position #for each footprint
		footprint_scores = [] #for each footprint

		for i in range(len(window)):
			position = i + 1 #calculate the relative position for a score
			score = window[i] #extract one score from the list
			if score >= bw_peak_background:
				if position != (check_position + 1): #if this position is not the next with the last position we have checked
					#save the last footprint
					if check_position != 0: #if this is not the start of the first footprint within this peak 

						footprint_count, peak_footprints = save_footprint(footprint_count, footprint_scores, peak_footprints, chromosom, footprint_start + peak_start + parts_positions[j], check_position + peak_start + parts_positions[j], bed_dictionary_entry)

					#start a new footprint
					footprint_start = position
					footprint_scores = []
						
					check_position = position

				footprint_scores.append(score) #save the current score
				check_position = position

		footprint_count, peak_footprints = save_footprint(footprint_count, footprint_scores, peak_footprints, chromosom, footprint_start + peak_start + parts_positions[j], check_position + peak_start + parts_positions[j], bed_dictionary_entry) #save the last footprint

	return peak_footprints, footprint_count

def find_peaks_from_bw(bed_dictionary, bw_file, window_length, step, percentage):

	logger.info('looking for footprints within peaks')

	footprint_count = 1
	all_footprints = {}

	bw_open = pyBigWig.open(bw_file)

	for header in bed_dictionary:

		peak_footprints = {}

		header_splitted = re.split(r':', header)
		chromosom = header_splitted[0]
		positions = re.split(r'-', header_splitted[-1])
		peak_start = int(positions[0])
		peak_end = int(positions[1])

		scores_in_peak = np.nan_to_num(np.array(list(bw_open.values(chromosom, peak_start, peak_end)))) #save the scores to an array

		peak_len = len(scores_in_peak)

		peak_footprints, footprint_count = search_in_window(peak_footprints, footprint_count, chromosom, peak_start, peak_end, scores_in_peak, window_length, bed_dictionary[header], step, percentage)

		for footprint_name in peak_footprints.keys():
			all_footprints[footprint_name] = all_footprints.get(footprint_name, {})
			all_footprints[footprint_name] = peak_footprints[footprint_name]

	all_footprints = sorted(all_footprints.items(), key = lambda x : (x[1]['chromosom'], x[1]['start']), reverse = False) 

	return all_footprints

def write_to_bed_file(all_footprints, sorted_output_file_name):
	output_file_name = "not_sorted_" + sorted_output_file_name #save in the working directory

	header = ["#chr", "start", "end", "name", "score", "len", "max_pos", "bonus_info"] #a header to know what is in the columns

	output_file = open(output_file_name, 'w') #open a file to write

	logger.info("print to the output file")

	output_file.write('\t'.join(header) + '\n') #write the header

	for footprint in all_footprints:
		output_file.write('\t'.join([footprint[1]['chromosom'], str(footprint[1]['start']), str(footprint[1]['end']), footprint[0], str(round(footprint[1]['score'], 6)), str(footprint[1]['len']), str(footprint[1]['max_pos']), '\t'.join(footprint[1]['bonus'])]) + '\n')

	output_file.close()

	#sort the bed file
	logger.info('sorting the output file')

	os.system("(head -n 2 " + output_file_name + " && tail -n +3 " + output_file_name + " | sort -k1,1V -k2,2n -k3,3n) > " + sorted_output_file_name)

	logger.info('remove the non-sorted file')

	remove_file(output_file_name)

def main():

	start = time.time()

	args = parse_args()

	check_existing_input_files(args)
	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	#check_directory(args.output_directory)

	#fh = logging.FileHandler(os.path.join(args.output_directory, "call_peaks_log.txt"))
	fh = logging.FileHandler("call_peaks_log.txt")
	fh.setLevel(logging.INFO)
	fh.setFormatter(formatter)
	logger.addHandler(fh)

	ch = logging.StreamHandler()
	ch.setLevel(logging.INFO)
	ch.setFormatter(formatter)
	logger.addHandler(ch)

	#if user do not want to see the information about the status of jobs, remove the handler, that writes to the terminal
	if args.silent:
		logger.removeHandler(ch)

	logger.info("call_peaks.py was called using these parameters:")
	logger.info(vars(args))

	bed_dictionary = make_bed_dictionary(args.bed)
	all_footprints = find_peaks_from_bw(bed_dictionary, args.bigwig, args.window_length, args.step, args.percentage)
	write_to_bed_file(all_footprints, args.output_file)

	logger.info("the number of peaks: " + str(len(bed_dictionary)))
	logger.info("the number of footprints: " + str(len(all_footprints)))

	logger.info("call_peaks needed %s minutes to generate the output" % (round((time.time() - start)/60, 2)))
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
	main()