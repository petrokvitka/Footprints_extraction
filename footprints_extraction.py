#!/usr/bin/env python

"""
footprints_extraction script uses the uncontinuous score from a bigWig file to estimate footpints within peaks of interest.
@author: Anastasiia Petrova
@contact: anastasiia.petrova(at)mpi-bn.mpg.de
"""

import argparse #for parsing the parameters
import sys #for example to exit if a problem occured
import os #for example to check the existing path of a file
import re #for example to split a string
import time #to calculate time needed to proceed the data
import logging #to write informaiton about the programm run
import numpy as np #to calculate mean for example
from collections import defaultdict #to create nested dictionaries
import pyBigWig #to work with bigWig files
import textwrap #to print the help message nice

#the logger should be callable from all functions, so create it here
logger = logging.getLogger('footprints_extraction')
logger.setLevel(logging.INFO)

formatter = logging.Formatter("%(asctime)s : %(message)s", "%Y-%m-%d %H:%M")

#the function parse_args() reads the input from user, including required and optional parameters
#if user defines no values for required parameters, the error will be printed and the script exits
#if user defines no values for optional parameters, the dafault values will be used
#this function returns a structure called args containing values for all parameters
def parse_args():
	
	parser = argparse.ArgumentParser(prog = '', description = textwrap.dedent('''                           

		This script produces a file with footprints in .bed format from the file with scores in .bigWig format and a corresponding .bed file with peaks of interest.
		'''), epilog='That is what you need to make this script work for you. Enjoy it')
	
	required_arguments = parser.add_argument_group('required arguments')
	required_arguments.add_argument('--bigwig', help='Please provide a path to the bigWig-file with scores.', required=True)
	required_arguments.add_argument('--bed', help='Please provide a path to the file with peaks in .bed format.', required=True)

	#all other arguments are optional
	#parser.add_argument('--output_directory',  default='output', const='output', nargs='?', help='The output directory, by default ./output/')
	parser.add_argument('--output_file', default='footprints_extraction.bed', type=str, help='Please enter a name for the output file, by default ./footprints_extraction.bed.')
	parser.add_argument('--window_length', default=200, type=int, help='Please enter the length for a window, by defauld 200 bp.')
	parser.add_argument('--step', default=100, type=int, help='Please enter a step to move the window, by default 100 bp.')
	parser.add_argument('--percentage',  default=0, type=int, help='Please enter a percentage to be added to background while searching for footprints, by default 0%%.')
	parser.add_argument('--min_gap', default=6, type=int, help='Please enter the number of bp allowed to be in between two footprints, by default 6 bp.')
	parser.add_argument('--silent', action='store_true', help='While working with data write the information only into ./footprints_extraction.log.')
	args = parser.parse_args()

	return args

#this function checks if the directory for the output files exists and if not, the new directory with default name will be created
#the input for this function is a directory to check
def check_directory(directory):
	if not os.path.exists(directory):
		os.makedirs(directory)
		logger.info('The new directory ' + directory + ' was created!')

#check if the file to remove exists and if so, delete it
#if the file does not exist, nothing happens
def remove_file(file):
	if os.path.isfile(file):
		os.remove(file)

#check if the input files exist, we will check for the right format of input files later on
#the input parameter args is a structure from parse_args, containing the input parameters
#if input files dont exist, the short message will be printed and the exit is forced
def check_existing_input_files(args):

	#check if the bigWig file exists
	if not os.path.isfile(args.bigwig):
		print('Error: there is no such bigWig file ' + args.bigwig + ', the exit is forced!')
		sys.exit()
	#check if the bed file exists
	if not os.path.isfile(args.bed):
		print('Error: there is no such bed file ' + args.bed + ', the exit is forced!')
		sys.exit()

#make a dictionary out of the .bed file for easy access to the information from the .bed file
#the input for this function is a path to the bed file
#the output is a dictionary containing key and values from the given bed file
#key looks like chr1:123-234, the value contains all the additional information from the bed file
def make_bed_dictionary(bed_file):

	logger.info('Reading of the bed file...')

	try: #if i can't proceede this file like so, the input was not a .bed file!
		bed_dictionary = {}

		with open(bed_file) as read_bed_file:
			for bed_line in read_bed_file:
				bed_line_array = re.split(r'\t', bed_line.rstrip('\n'))
				if bed_line_array[1].isdigit() and bed_line_array[2].isdigit() and int(bed_line_array[1]) <= int(bed_line_array[2]): #in the real bedfile the second column is a start position, and the third column is an end position, so we are checking if these are integers and if the start position is smaller than the end one
					key = bed_line_array[0] + ":" + bed_line_array[1] + "-" + bed_line_array[2]
					value = []
					for i in range(3, len(bed_line_array)):
						value.append(bed_line_array[i]) 

					if not value: #if there is no bonus information in the original bed file, add a "." to the coulmn in the output bed file
						value = ["."]

					bed_dictionary[key] = value
				else: #this is not a bed file, force exit
					logger.info('Error: please make sure the input bed file has a right format, the problem occured on the line ' + bed_line)
					sys.exit()

		read_bed_file.close()

		return bed_dictionary

	except UnicodeDecodeError: #force exit, if input is not a .bed file
		logger.info('Error: please make sure that the .bed file has a right format! The exit is forced!')
		sys.exit()

#to save a footprint, find the score and max_pos for this footprint and check for overlapping with other footprints within the current peak
#the input for this function is the count for footprint to ensure unique name for each saved footprint;
#an array containing scores (signals) from the bigwig file; the footprints already saved for the particular peak; information for the bed file:
#chromosom, start and end position, as well as additional information from the original bed file
#the function returns the count for footprints, as well as footprints for the current peak
def save_footprint(footprint_count, footprint_scores, peak_footprints, chromosom, footprint_start, footprint_end, bonus_info_from_bed): 

	save_current_footprint = False

	#make sure we are working only with regions bigger than 2, in other case we have not enough information to find the footprint
	if len(footprint_scores) > 2: 
		
		#calculate the position with the max score. if there are many positions with the same score, save one from the middle
		first_max_pos = footprint_scores.index(max(footprint_scores))
		#assume that there is only one pos with max score
		last_max_pos = first_max_pos

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

		max_pos = max_pos + 1 #as the index of an array starts with 0

		#calculate the score for the current footprint as mean of all scores from the bigwig file
		footprint_score = np.mean(footprint_scores)

		#checking for existing and overlapping footprints		
		if len(peak_footprints.keys()) == 0: # hey, this is the first footprint!
			save_current_footprint = True

		else: # there are already footprints in this peak, so look for overlaps
			for existing_footprint_name in peak_footprints.keys():

				old_start = peak_footprints[existing_footprint_name]['start']
				old_end = peak_footprints[existing_footprint_name]['end']
				old_score = peak_footprints[existing_footprint_name]['score']
				old_length = peak_footprints[existing_footprint_name]['len']

				if footprint_start >= old_start and footprint_start <= old_end: #the start of the new footprint is between the start and end of an old footprint
					if footprint_end > old_end: #the new footprint is not completely inside the old one
						#update the information about the existing footprint
						#find the average of both scores
						footprint_score = (peak_footprints[existing_footprint_name]['score'] + footprint_score) / 2

						peak_footprints[existing_footprint_name]['end'] = footprint_end
						peak_footprints[existing_footprint_name]['score'] = footprint_score
						peak_footprints[existing_footprint_name]['len'] = footprint_end - old_start
						#we can not update the max_pos as we do not have the information about scores array of the existing footprint
					#else: the new footprint is completely inside the old one, do nothing
					save_current_footprint = False
					break

				elif footprint_end >= old_start and footprint_end <= old_end: #the end of the new footprint is between the start and end of an old footprint
					if footprint_start < old_start: #the new footprint is not completely inside the old one						
						#update the information about the existing footprint
						#find the average of both scores
						footprint_score = (peak_footprints[existing_footprint_name]['score'] + footprint_score) / 2

						peak_footprints[existing_footprint_name]['start'] = footprint_start
						peak_footprints[existing_footprint_name]['score'] = footprint_score
						peak_footprints[existing_footprint_name]['len'] = old_end - footprint_start
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
			footprint_name = "f_" + str(footprint_count)

			peak_footprints[footprint_name] = peak_footprints.get(footprint_name, {})
			peak_footprints[footprint_name] = {'chromosom': chromosom, 'start': footprint_start, 'end': footprint_end, 'score': footprint_score, 'len': len(footprint_scores), 'bonus': bonus_info_from_bed, 'max_pos': max_pos}

			footprint_count = footprint_count + 1

	#else do nothing

	return footprint_count, peak_footprints

#this function uses a slide window algorithm to estimate regions where the signal is higher than the background and can possibly be footprints
#as an input this function takes: the footprints already saved for the current peak; the count for footprints to ensure unique name;
#the information for the bed file (chromosom, start- and end-positions of peak); an array with scores within the peak; the length of the window;
#as well as the dictionary containing the information from the original bed file and the optional parameters such as step to slide the window and percentage to change the background
#this function returns the array containing the found footprints within a peak as well as the footprint count to look for footprints in the next peak and ensure the unique names for footprints
def search_in_window(peak_footprints, footprint_count, chromosom, peak_start, peak_end, scores_in_peak, window_length, bed_dictionary_entry, step, percentage):
	
	#find the length of the current peak
	peak_len = len(scores_in_peak)
	parts = [] #an array to save scores of parts
	parts_positions = [] #an array to save where each part begins in the original peak

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

	#look in each window and try to save the footprints
	for j in range(len(parts)):
		window = parts[j]

		#add some percent to the background to avoid the noice we don't want to have
		bw_peak_background = np.mean(window) #find the mean of all scores within one peak
		part = (percentage*bw_peak_background)/100 #x procent of the background
		bw_peak_background = bw_peak_background + part
		
		check_position = parts_positions[j] #the start position not within the window, but within the peak!!!
		footprint_start = check_position #for each footprint
		footprint_scores = [] #for each footprint

		#look on each position inside the window
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

#this function checks if there are footprints with a small amount of bp in between, and if so, merges them together
#if at this point there are still overlaps of footprints, this function will delete them
#the input parameter are: dictionary with footprints within one peak, and the max number of bp allowed to be in between the footprints
#the output is the renewed dictionary containing only the best footprints for the output file
def check_and_merge(peak_footprints, min_gap):
	#to ensure the merging works well, sort the footprints first arter start and end positions
	#the sort can not be applied to a dictionary, we are making a list out of peak_footprints_dict
	peak_footprints_list = sorted(peak_footprints.items(), key = lambda x : (x[1]['start'], x[1]['end']), reverse = False)

	peak_footprints_new = {}
	merged_footprints = {}

	#we need to check each footprint within this peak with the other footprints for possible merging
	#for footprint_to_check in peak_footprints.keys():
	for footprint in peak_footprints_list: #work with sorted footprints
		footprint_to_check = footprint[0] #save the name of the footprint which we are working with now
		start_to_check = peak_footprints[footprint_to_check]['start']
		end_to_check = peak_footprints[footprint_to_check]['end']

		merge_footprints_left = None

		for compared_footprint in peak_footprints.keys():

			if start_to_check > peak_footprints[compared_footprint]['start'] and start_to_check - peak_footprints[compared_footprint]['end'] < min_gap:
				#make compared_footprint longer: compared_footprint + footprint_to_check
				merge_footprints_left = False
				break
			elif end_to_check < peak_footprints[compared_footprint]['end'] and peak_footprints[compared_footprint]['start'] - end_to_check < min_gap:
				#make footprint_to_check longer: footprint_to_check + compared footprint
				merge_footprints_left = True
				break

		if merge_footprints_left: #the left merging is enabled, start and end of compared_footprint should be smaller than the start of the footprint_to_check
			if start_to_check < peak_footprints[compared_footprint]['start']:
				if footprint_to_check not in peak_footprints_new.keys():
					if any(footprint_to_check in merged_footprints[x] for x in merged_footprints.keys()): #true if footprint_to_check was already merged with someone
						for k, v in merged_footprints.items():
							if footprint_to_check in v:
								main_footprint = k
						#make merging using the information from the merged_footprints and peak_footprints_new
						#UPDATE
						peak_footprints_new[main_footprint] = footprint_update(peak_footprints_new[main_footprint], peak_footprints[compared_footprint]['start'], peak_footprints[main_footprint]['end'], peak_footprints[compared_footprint]['score'])
						merged_array = merged_footprints[main_footprint]
						merged_array.append(compared_footprint)
						merged_footprints[main_footprint] = merged_array
					#there are no merged footprints with the footprint_to_check yet, so make a new one
					else:
						#add the compared footprint and footprint_to_check to the merged_footprints
						merged_footprints[footprint_to_check] = merged_footprints.get(footprint_to_check, [])
						merged_footprints[footprint_to_check] = [compared_footprint]

						peak_footprints_new[footprint_to_check] = peak_footprints.get(footprint_to_check, {})
						peak_footprints_new[footprint_to_check] = peak_footprints[footprint_to_check] #<-- update
						#UPDATE
						peak_footprints_new[footprint_to_check] = footprint_update(peak_footprints_new[footprint_to_check], peak_footprints_new[footprint_to_check]['start'], peak_footprints[compared_footprint]['end'], peak_footprints[compared_footprint]['score'])
				else: #the footprint_to_check is in peak_footprints_new already
					#the footprint_to_check can only be the main part of merging before, check it
					if footprint_to_check in merged_footprints.keys():
						#footprint_to_check was as main for merging already
						#UPDATE
						peak_footprints_new[footprint_to_check] = footprint_update(peak_footprints_new[footprint_to_check], peak_footprints_new[footprint_to_check]['start'], peak_footprints[compared_footprint]['end'], peak_footprints[compared_footprint]['score'])
						#add it to the merged_footprints as well
						merged_array = merged_footprints[footprint_to_check]
						merged_array.append(compared_footprint)
						merged_footprints[footprint_to_check] = merged_array
					else:
						#the footprint_to check was not merged with anything yet
						merged_footprints[footprint_to_check] = merged_footprints.get(footprint_to_check, [])
						merged_footprints[footprint_to_check] = [compared_footprint]
						#UPDATE
						peak_footprints_new[footprint_to_check] = footprint_update(peak_footprints_new[footprint_to_check], peak_footprints_new[footprint_to_check]['start'], peak_footprints[compared_footprint]['end'], peak_footprints[compared_footprint]['score'])
		#the right merging is enabled, start and end of compared footprint should be bigger than the start of the footprint_to_check
		elif merge_footprints_left == False: 
			if end_to_check > peak_footprints[compared_footprint]['end']:
				if compared_footprint not in peak_footprints_new.keys():				
					if any(compared_footprint in merged_footprints[x] for x in merged_footprints.keys()):
						for k, v in merged_footprints.items():
							if compared_footprint in v:
								main_footprint = k
						#make merging using the information from the merged_footprints and peak_footprints_new
						#UPDATE
						peak_footprints_new[main_footprint] = footprint_update(peak_footprints_new[main_footprint], peak_footprints[main_footprint]['start'], peak_footprints[footprint_to_check]['end'], peak_footprints[footprint_to_check]['score'])
						#add to the merged_footprints
						merged_array = merged_footprints[main_footprint]
						merged_array.append(footprint_to_check)
						merged_footprints[main_footprint] = merged_array
					else:
						#make normal update, using data from peak footprints
						merged_footprints[compared_footprint] = merged_footprints.get(compared_footprint, [])
						merged_footprints[compared_footprint] = [footprint_to_check]

						peak_footprints_new[compared_footprint] = peak_footprints.get(compared_footprint, {})
						peak_footprints_new[compared_footprint] = peak_footprints[compared_footprint]
						#UPDATE
						peak_footprints_new[compared_footprint] = footprint_update(peak_footprints_new[compared_footprint], peak_footprints_new[compared_footprint]['start'], peak_footprints[footprint_to_check]['end'], peak_footprints[footprint_to_check]['score'])
				else:
					if compared_footprint in merged_footprints.keys():
						#compared_footprint was as main for merging already
						#UPDATE
						peak_footprints_new[compared_footprint] = footprint_update(peak_footprints_new[compared_footprint], peak_footprints_new[compared_footprint]['start'], peak_footprints[footprint_to_check]['end'], peak_footprints[footprint_to_check]['score'])

						merged_array = merged_footprints[compared_footprint]
						merged_array.append(footprint_to_check)
						merged_footprints[compared_footprint] = merged_array

					else:
						#merge now and add compared_footprint to the merged_footprints
						merged_footprints[compared_footprint] = merged_footprints.get(compared_footprint, [])
						merged_footprints[compared_footprint] = [compared_footprint]
						#UPDATE
						peak_footprints_new[compared_footprint] = footprint_update(peak_footprints_new[compared_footprint], peak_footprints_new[compared_footprint]['start'], peak_footprints[footprint_to_check]['end'], peak_footprints[footprint_to_check]['score'])

		else: #save the current footprint, as it should not be merged
			peak_footprints_new[footprint_to_check] = peak_footprints_new.get(footprint_to_check, [])
			peak_footprints_new[footprint_to_check] = peak_footprints[footprint_to_check]

	return peak_footprints_new

#this function is used to update the footprint that should to be merged with another one
#as input the footprint, needed the update, as well as new start, new end and the score of the merged footprint are passed
#the output of this function is a dictionary containing the new information about the footprint
def footprint_update(footprint, start, end, score):
	new_len = end - start
	new_score = (footprint['score'] + score) / 2

	footprint['start'] = start
	footprint['end'] = end
	footprint['score'] = new_score
	footprint['len'] = new_len

	return footprint

#this function uses the information provided from the .bed file to look for footprints within the peaks of interest
#as input the information from the original bed file, as well as bigwig file is needed
#the optional parameters window_length, step and percentage are needed as well to use the sliding window algorithm and work with the "background" score
#the output of this function is a dictionary contains all the found footprints ready to write to the output file
def find_peaks_from_bw(bed_dictionary, bw_file, window_length, step, percentage, min_gap):

	logger.info('Looking for footprints within peaks...')

	try: #this try considers opening the bigwig file
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
			peak_footprints, footprint_count = search_in_window(peak_footprints, footprint_count, chromosom, peak_start, peak_end, scores_in_peak, window_length, bed_dictionary[header], step, percentage)

			#double check for overlaps and possibly merging of footprints having up to 5 bp in between
			peak_footprints = check_and_merge(peak_footprints, min_gap)

			for footprint_name in peak_footprints.keys():
				all_footprints[footprint_name] = all_footprints.get(footprint_name, {})
				all_footprints[footprint_name] = peak_footprints[footprint_name]

		all_footprints = sorted(all_footprints.items(), key = lambda x : (x[1]['chromosom'], x[1]['start']), reverse = False) 

		return all_footprints

	except RuntimeError: #if i can't work with the bigwig file like so, it was not a bigwig file!
		logger.info('Error: please make sure that the .bigWig file has a right format! The exit is forced!')
		sys.exit()

#this function writes the found footprints to the .bed file
#as input the dictionary with all footprints and the name for the output file are needed
#the function outputs first an unsorted file and finally sorts it, removing the unsorted one
def write_to_bed_file(all_footprints, sorted_output_file_name):
	
	#extract the name of directory where we want to save the file
	output_directory = os.path.dirname(sorted_output_file_name)
	#check if there is a directory provided
	if output_directory:
		check_directory(output_directory)
	
	#make sure the file is in .bed format
	output_name = "not_sorted_" + os.path.splitext(os.path.basename(sorted_output_file_name))[0] + ".bed"

	#save the output file in the working directory or in the directory provided by the user
	output_file_name = (os.path.join(output_directory, output_name))
	
	#a header to know what is in the columns
	header = ["#chr", "start", "end", "name", "score", "strand", "len", "max_pos", "bonus_info"]

	#open a file to write
	output_file = open(output_file_name, 'w')

	logger.info("Printing to the output file...")

	#write the header
	output_file.write('\t'.join(header) + '\n')

	#write each footprint line for line to the output file
	for footprint in all_footprints:
		#validation of the footprints, if there is a problem with some of them, write which one it is
		#first check the start and end positions
		if footprint[1]['start'] >= footprint[1]['end']:
			logger.info("The problem occured with start and end positions. This footprint will not be printed to the output file:")
			logger.info(footprint)
		#then check the max_pos
		elif footprint[1]['max_pos'] == 0:
			logger.info("The problem occured with max_pos of the footprint. This footprint will not be printed to the output file:")
			logger.info(footprint)
		#otherwise everything is fine, write to the output
		else:
			output_file.write('\t'.join([footprint[1]['chromosom'], str(footprint[1]['start']), str(footprint[1]['end']), footprint[0], str(round(footprint[1]['score'], 6)), '.', str(footprint[1]['len']), str(footprint[1]['max_pos']), ';'.join(footprint[1]['bonus'])]) + '\n')

	output_file.close()

	#sort the bed file
	logger.info('Sorting the output file...')

	os.system("(head -n 2 " + output_file_name + " && tail -n +3 " + output_file_name + " | sort -k1,1V -k2,2n -k3,3n) > " + sorted_output_file_name)

	logger.info('Removing the non-sorted file...')

	remove_file(output_file_name)

def main():

	start = time.time()

	args = parse_args()

	check_existing_input_files(args)
	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	#check_directory(args.output_directory)

	#fh = logging.FileHandler(os.path.join(args.output_directory, "footprints_extraction.log"))
	fh = logging.FileHandler("footprints_extraction.log")
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

	logger.info("The script footprints_extraction.py was called using these parameters: " + str(vars(args)))
	
	bed_dictionary = make_bed_dictionary(args.bed)
	all_footprints = find_peaks_from_bw(bed_dictionary, args.bigwig, args.window_length, args.step, args.percentage, args.min_gap)
	write_to_bed_file(all_footprints, args.output_file)

	logger.info("the number of peaks: " + str(len(bed_dictionary)))
	logger.info("the number of footprints: " + str(len(all_footprints)))

	logger.info("It took footprints_extraction.py %s minutes to generate the output." % (round((time.time() - start)/60, 2)))
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
	main()
