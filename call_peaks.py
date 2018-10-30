
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
import multiprocessing
import logging
import subprocess
from Bio import SeqIO
import Bio.SeqIO.FastaIO as bio
import numpy as np
from collections import defaultdict
from scipy import stats
import pyBigWig
from statsmodels.sandbox.stats.multicomp import multipletests #for bonfferoni
import matplotlib.pyplot as plt
import random

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
	parser.add_argument('--output_directory',  default='output', const='output', nargs='?', help='output directory, by default ./output/')
	parser.add_argument('--window_length', default='100', type=int, help='enter the length for a window, by defauld 100 bp')
	parser.add_argument('--threshold',  default=0.3, type=float, help='enter the threshold for peaks searching, by defauld 0.3')
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

#check if the bw score is already saved, if so check if it is bigger than the new one
def save_bw_score(key_for_bed_dict, matches_dict, bw_score, moods_score, which_score):

	if np.isnan(bw_score): bw_score = 0.0

	#bw_score = moods_score * bw_score #apply moods score as well 

	if key_for_bed_dict in matches_dict:

		if which_score == "mean":
			#save the mean of both scores
			matches_dict[key_for_bed_dict] = np.mean([matches_dict[key_for_bed_dict], bw_score])
		
		elif which_score == "greater":
			#save the biggest of scores
			if matches_dict[key_for_bed_dict] < bw_score:
				matches_dict[key_for_bed_dict] = bw_score

	else:
		matches_dict[key_for_bed_dict] = bw_score

	return matches_dict

def remove_file(file):
	if os.path.isfile(file):
		os.remove(file)

def clean_directory(cleans, output_directory, motif):
	
	for clean in cleans:
		if clean == 'all' or clean == 'cut_motifs':
			remove_file(motif)
			#the output_merge.fa will be deleted after processing of the multiprocessing

def tool_make_output(motif, genome, output_directory, cleans, p_value, bed_dictionary, moods_bg, condition1, condition2, global_mean, global_std, which_score):

	standard_moods_bg = 0.25, 0.25, 0.25, 0.25

	control_dict = {}
	overexpression_dict = {}

	differences = []
	control_dict, overexpression_dict, differences = call_moods(motif, genome, output_directory, p_value, standard_moods_bg, condition2, condition1, control_dict, overexpression_dict, differences, which_score) 
	
	#make arrays of the dictionaries
	control_array = []
	overexpression_array = []

	for key in control_dict:
		control_array.append(control_dict[key])
		overexpression_array.append(overexpression_dict[key])

	motif_name = motif.replace(output_directory, '')
	motif_name = motif_name.replace('/', '')

	#make the wilcoxon signed-rank test
	my_wilcoxon_pvalue, direction, differences, differences_normalized, motif_std = my_wilcoxon(condition2, condition1, control_array, overexpression_array, global_mean, global_std, motif_name, correction = False)

	clean_directory(cleans, output_directory, motif)

	return my_wilcoxon_pvalue, direction, differences, differences_normalized, motif_std

def make_name_from_path(full_path, output_directory, ending):
	return os.path.join(output_directory, get_name_from_path(full_path) + ending)
	
def get_name_from_path(full_path):
	return os.path.splitext(os.path.basename(full_path))[0]

def compute_differences(bed_dictionary, condition1, condition2):
	logger.info("the mean and standard deviation for the differences in peaks will be count now")

	bw_cond1 = pyBigWig.open(condition1)
	bw_cond2 = pyBigWig.open(condition2)

	global_differences = {} #dict
	differences_array = [] #to compute the mean at the end
	cond1_array = []
	cond2_array = []

	for header in bed_dictionary:
		header_splitted = re.split(r':', header)
		chromosom = header_splitted[0]
		positions = re.split(r'-', header_splitted[-1])

		#compute the background difference for this peak
		bw_global_score_cond1 = np.mean(np.nan_to_num(np.array(list(bw_cond1.values(chromosom, int(positions[0]), int(positions[1]))))))
		bw_global_score_cond2 = np.mean(np.nan_to_num(np.array(list(bw_cond2.values(chromosom, int(positions[0]), int(positions[1]))))))
		bw_global_difference = bw_global_score_cond2 - bw_global_score_cond1

		global_differences[header] = bw_global_difference

		cond1_array.append(bw_global_score_cond1)
		cond2_array.append(bw_global_score_cond2)

		differences_array.append(bw_global_difference)

	bw_cond1.close()
	bw_cond2.close()

	mu = np.mean(differences_array)
	std = np.std(differences_array, ddof = 1)

	mu_cond1 = np.mean(cond1_array)

	mu_cond2 = np.mean(cond2_array)

	cond1_name = get_name_from_path(condition1)
	cond2_name = get_name_from_path(condition2)

	return mu, std

def is_fasta(check_fasta):
	if not os.path.isfile(check_fasta):
		#logger.info('there is no file with genome, the exit is forced')
		print('there is no file with genome, the exit is forced')
		sys.exit()
	else:
		# modified code from https://stackoverflow.com/questions/44293407/how-can-i-check-whether-a-given-file-is-fasta
	    with open(check_fasta, "r") as handle:
	        fasta = SeqIO.parse(handle, "fasta")
	        return any(fasta)  # False when `fasta` is empty, i.e. wasn't a FASTA file

def check_existing_input_files(args):

	if not is_fasta(args.genome):
		#logger.info('please make sure the input genome file has a fasta format')
		print('please make sure the input genome file has a fasta format')
		sys.exit()
	if not os.path.isfile(args.condition1) or not os.path.isfile(args.condition2):
		#logger.info('please make sure the both files with conditions to compare exist')
		print('please make sure the both files with conditions to compare exist')
		sys.exit()
	if not args.condition1.endswith('.bw') or not args.condition2.endswith('.bw'):
		#logger.info('please check if the both conditions files are in bigWig format')
		print('please check if the both conditions files are in bigWig format')
		sys.exit()
	#check if the file with motifs exists
	if not os.path.isfile(args.motifs):
		#logger.info('there is no file with motifs, the exit is forced')
		print('there is no file with motifs, the exit is forced')
		sys.exit()
	#check if the bed file exists
	if not os.path.isfile(args.bed_file):
		#logger.info('there is no such bed file ' + args.bed_file + ', the exit is forced')
		print('there is no such bed file ' + args.bed_file + ', the exit is forced')
		sys.exit()

def make_bed_dictionary(bed_file):

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

def find_window(bed_file):
	
	#chromosom = 1 #start with the first chromosom
	window_length = 0

	with open(bed_file) as read_bed_file:
		for bed_line in read_bed_file:
			bed_line_array = re.split(r'\t', bed_line.rstrip('\n'))
			if bed_line_array[1].isdigit() and bed_line_array[2].isdigit() and int(bed_line_array[1]) <= int(bed_line_array[2]): #in the real bedfile the second column is a start position, and the third column is an end position, so we are checking if these are integers and if the start position is smaller than the end one
				#if chromosom == int(bed_line_array[0]):
				peak_len = int(bed_line_array[2]) - int(bed_line_array[1])
				if peak_len > window_length:
					window_length = peak_len

	print(window_length)

def save_footprint(footprint_count, footprint_scores, all_footprints, chromosom, footprint_start, check_position, bonus_info_from_bed):

	footprint_name = "footprint_" + str(footprint_count)

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

	footprint_score = np.mean(footprint_scores)
	all_footprints[footprint_name] = all_footprints.get(footprint_name, {})
	all_footprints[footprint_name] = {'chromosom': chromosom, 'start': footprint_start, 'end': check_position, 'score': footprint_score, 'len': len(footprint_scores), 'bonus': bonus_info_from_bed, 'max_pos': max_pos}
	footprint_count += 1

	return footprint_count, all_footprints

def find_peaks_from_bw(bed_dictionary, bw_file):

	footprint_count = 1
	all_footprints = {}

	bw_open = pyBigWig.open(bw_file)

	for header in bed_dictionary:
		header_splitted = re.split(r':', header)
		chromosom = header_splitted[0]
		positions = re.split(r'-', header_splitted[-1])

		scores_in_peak = np.nan_to_num(np.array(list(bw_open.values(chromosom, int(positions[0]), int(positions[1]))))) #save the scores to an array

		bw_peak_background = np.mean(scores_in_peak) #find the mean of all scores within one peak
		part = bw_peak_background/10 #10 procent of the background
		bw_peak_background = bw_peak_background + part

		check_position = 0 #for the whole peak
		footprint_start = 1 #for each footprint
		footprint_scores = [] #for each footprint

		for i in range(len(scores_in_peak)):
			position = i + 1 #calculate the relative position for a score
			score = scores_in_peak[i] #extract one score from the list
			if score >= bw_peak_background:
				if position != (check_position + 1): #if this position is not the next with the last position we have checked
					#save the last footprint
					if check_position != 0: #if this is not the start of the first footprint within this peak 
						
						footprint_count, all_footprints = save_footprint(footprint_count, footprint_scores, all_footprints, chromosom, footprint_start, check_position, bed_dictionary[header])

					#start a new footprint
					footprint_start = position
					footprint_scores = []
					
					check_position = position

				footprint_scores.append(score) #save the current score
				check_position = position

		footprint_count, all_footprints = save_footprint(footprint_count, footprint_scores, all_footprints, chromosom, footprint_start, check_position, bed_dictionary[header])

		all_footprints = sorted(all_footprints.items(), key =  lambda x : (x[1]['start']), reverse = False) 

		return all_footprints

def write_to_bed_file(all_footprints):
	output_file_name = "footprints.bed" #save in the working directory

	header = ["#chr", "start", "end", "name", "len", "score", "max_pos", "bonus_info"] #a header to know what is in the columns

	output_file = open(output_file_name, 'w') #open a file to write

	output_file.write('\t'.join(header) + '\n') #write the header

	for footprint in all_footprints:
		output_file.write('\t'.join([footprint[1]['chromosom'], str(footprint[1]['start']), str(footprint[1]['end']), footprint[0], str(footprint[1]['len']), str(footprint[1]['score']), str(footprint[1]['max_pos']), '\t'.join(footprint[1]['bonus'])]) + '\n')

	output_file.close()


def main():

	start = time.time()

	peaks_bed_file = "./small_peaks.bed"
	#peaks_bed_file = "./control_peaks.bed"
	#find_window(peaks_bed_file)

	#bed_dictionary = make_bed_dictionary(peaks_bed_file)
	bed_dictionary = {}
	bed_dictionary["chr1:3062743-3063132"] = ["control_1", "test", "hihi"]
	#bed_dictionary["chr1:3343546-3344520"] = ["control_2"]
	#bed_dictionary["chr1:3062810-3063132"] = ["control1"] #the 0.position has already a score bigger than the background


	bw_file = "./control_footprints.bw"

	all_footprints = find_peaks_from_bw(bed_dictionary, bw_file)
	write_to_bed_file(all_footprints)

	#args = parse_args()

	#check_existing_input_files(args)
	#check if there is an existing directory that user gave as input, otherwise create this directory from the path provided from the user
	#check_directory(args.output_directory)

	#fh = logging.FileHandler(os.path.join(args.output_directory, "call_peaks_log.txt"))
	#fh.setLevel(logging.INFO)
	#fh.setFormatter(formatter)
	#logger.addHandler(fh)

	#ch = logging.StreamHandler()
	#ch.setLevel(logging.INFO)
	#ch.setFormatter(formatter)
	#logger.addHandler(ch)

	#if user do not want to see the information about the status of jobs, remove the handler, that writes to the terminal
	#if args.silent:
	#	logger.removeHandler(ch)

	#logger.info("call_peaks.py was called using these parameters:")
	#logger.info(vars(args))

	#blablablaaaa

	#logger.info("call_peaks needed %s seconds to generate the output" % (time.time() - start))
	
	for handler in logger.handlers:
		handler.close()
		logger.removeFilter(handler)
	
if __name__ == "__main__":
	main()