#!/usr/bin/env python

"""This module roughly follows the algorithm used by MACS to call
ChIP-Seq peaks, but applies it to calling card data.  The main peak calling
function is passed a experiment frame, a background frame, and an TTAA_frame, 
all in ccf format.  It then builds interval trees containing
all of the background and experiment hops and all of the TTAAs.  Next, 
it scans through the genome with a window of window_size and step size of step_size 
and looks for regions that have signficantly more experiment hops than background 
hops (poisson w/ pvalue_cutoff).  It merges consecutively enriched windows and 
computes the center of the peak.  Next it computes lambda, the number 
of insertions per TTAA expected from the background distribution by taking the max of 
lambda_bg, lamda_1, lamda_5, lambda_10.  It then computes a p-value based 
on the expected number of hops = lamda * number of TTAAs in peak * number of hops
in peak.  Finally, it returns a frame that has Chr,Start,End,Center,Experiment Hops,
Fraction Experiment,Background Hops,Fraction Background,Poisson pvalue

There is also a background free version of this algorithm that computes significance 
based on the number of hops in the neighboring region.

"""

import numpy as np
import pandas as pd
import scipy.stats as scistat
from bx.intervals.intersection import Intersecter, Interval
import pybedtools.bedtool as bt 
import argparse
import re
import annotate
import sys

#######################
#Master Loop ##########
#######################

def call_peaks_and_annotate(expfile,bgfile,outputfile,TTAA_file = '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt',
	annotation_file = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed',pvalue_cutoff=1e-4,
	peak_pvalue_cutoff = 1e-3,window_size = 1000, step_size = 500,pseudocounts = 0.2):

	"""This function is the wrapper that calls the various subfunctions"""

	expframe = pd.read_csv(expfile,delimiter="\t",header=None)
	expframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	bgframe = pd.read_csv(bgfile,delimiter="\t",header=None)
	bgframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	TTAAframe = pd.read_csv(TTAA_file,delimiter="\t",header=None)
	TTAAframe.columns = ["Chr","Start","End"]

	peaks_frame = find_peaks(expframe,bgframe,TTAAframe,peak_pvalue_cutoff,window_size,step_size,pseudocounts)
	peaks_frame = annotate.annotate_peaks_frame(peaks_frame,annotation_file)
	peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
	peaks_frame = peaks_frame.sort_values(["Poisson pvalue"])
	peaks_frame.to_csv(outputfile,sep="\t",index=False)
	peaksbed_frame = annotate.make_peaksbed(peaks_frame)
	
	peaksbed_frame.columns = ['Chr','Start','End','Col4']
	peaksbed_frame_with_description = pd.merge(peaksbed_frame,peaks_frame,how = "left", on = ['Chr','Start','End'])
	peaksbed_frame_with_description['TPH Experiment'] = 'TPH Experiment=' + peaksbed_frame_with_description['TPH Experiment'].astype(str)
	peaksbed_frame_with_description['TPH Background'] = 'TPH Background=' + peaksbed_frame_with_description['TPH Background'].astype(str)
	peaksbed_frame_with_description['TPH Background subtracted'] = 'TPH Background subtracted=' + peaksbed_frame_with_description['TPH Background subtracted'].astype(str)
	peaksbed_frame_with_description['Poisson pvalue'] = 'Poisson pvalue=' + peaksbed_frame_with_description['Poisson pvalue'].astype(str)
	peaksbed_frame_with_description['Lambda'] = 'Lambda=' + peaksbed_frame_with_description['Lambda'].astype(str)
	peaksbed_frame_with_description['Lambda Type'] = 'Lambda Type=' + peaksbed_frame_with_description['Lambda Type'].astype(str)

	peaksbed_frame_with_description['Description'] = list(peaksbed_frame_with_description.loc[:,['TPH Experiment','TPH Background','TPH Background subtracted','Poisson pvalue','Lambda','Lambda Type']].values)
	
	final_peaksbed_frame = peaksbed_frame_with_description.loc[:,['Chr','Start','End','Col4','Description']]

	peaksbed_frame = final_peaksbed_frame

	pattern = "^(\S+)\."
	group = re.search(pattern,outputfile)
	bedfilename = group.groups(0)[0]+".tmp.peaks.bed"
	peaksbed_frame.to_csv(bedfilename,sep="\t",index=False,header=False)

	df = pd.read_csv(bedfilename,sep = "\t",header = None)
	df.columns = ['Chr','Start','End','Col4','Description']
	df['Description'] = df['Description'].str.replace('\n','')

	bedfilename = group.groups(0)[0]+".peaks.bed"
	df.to_csv(bedfilename,sep = "\t",header = None,index = False)


######################################################
#Master Loop for Background Free Peak Calling ########
######################################################


def call_peaks_and_annotate_bf(expfile,outputfile,TTAA_file = '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt',
	annotation_file = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed',pvalue_cutoff=1e-4,
	peak_pvalue_cutoff = 1e-3,window_size = 1000, lam_win_size=100000, step_size = 500,pseudocounts = 0.2):

	"""This function is the wrapper that calls the various subfunctions"""

	expframe = pd.read_csv(expfile,delimiter="\t",header=None)
	expframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	TTAAframe = pd.read_csv(TTAA_file,delimiter="\t",header=None)
	TTAAframe.columns = ["Chr","Start","End"]

	peaks_frame = find_peaks_bf(expframe,TTAAframe,peak_pvalue_cutoff,window_size,lam_win_size,step_size,pseudocounts)
	peaks_frame = annotate.annotate_peaks_frame(peaks_frame,annotation_file)
	peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
	peaks_frame = peaks_frame.sort_values(["Poisson pvalue"])
	peaks_frame.to_csv(outputfile,sep="\t",index=False)
	peaksbed_frame = annotate.make_peaksbed(peaks_frame)

	peaksbed_frame.columns = ['Chr','Start','End','Col4']
        peaksbed_frame_with_description = pd.merge(peaksbed_frame,peaks_frame,how = "left", on = ['Chr','Start','End'])
        peaksbed_frame_with_description['TPH Experiment'] = 'TPH Experiment=' + peaksbed_frame_with_description['TPH Experiment'].astype(str)
        peaksbed_frame_with_description['Poisson pvalue'] = 'Poisson pvalue=' + peaksbed_frame_with_description['Poisson pvalue'].astype(str)
        peaksbed_frame_with_description['Lambda'] = 'Lambda=' + peaksbed_frame_with_description['Lambda'].astype(str)
        peaksbed_frame_with_description['Lambda Type'] = 'Lambda Type=' + peaksbed_frame_with_description['Lambda Type'].astype(str)

        peaksbed_frame_with_description['Description'] = list(peaksbed_frame_with_description.loc[:,['TPH Experiment','Poisson pvalue','Lambda','Lambda Type']].values)

        final_peaksbed_frame = peaksbed_frame_with_description.loc[:,['Chr','Start','End','Col4','Description']]

        peaksbed_frame = final_peaksbed_frame

        pattern = "^(\S+)\."
        group = re.search(pattern,outputfile)
        bedfilename = group.groups(0)[0]+".tmp.peaks.bed"
        peaksbed_frame.to_csv(bedfilename,sep="\t",index=False,header=False)

        df = pd.read_csv(bedfilename,sep = "\t",header = None)
        df.columns = ['Chr','Start','End','Col4','Description']
        df['Description'] = df['Description'].str.replace('\n','')

        bedfilename = group.groups(0)[0]+".peaks.bed"
        df.to_csv(bedfilename,sep = "\t",header = None,index = False)


#######################
#The actual peakcaller#
#######################

def find_peaks(experiment_frame,background_frame,TTAA_frame,
	pvalue_cutoff = 1e-3,window_size = 1000, step_size = 500,
	pseudocounts = 0.2):

	"""This function is passed an experiment frame, a background frame, 
	and an TTAA_frame,all in ccf format:
	(Chr, Start, Stop, Reads, Strand, Barcode)

	It then builds interval trees containing all of the background and 
	experiment hops and all of the TTAAs.  Next, it scans through the 
	genome with a window of window_size and step size of step_size and 
	looks for regions that have signficantly more experiment hops than
	background hops (poisson w/ pvalue_cutoff). It merges consecutively 
	enriched windows and computes the center of the peak.  Next 
	it computes lambda, the number of insertions per TTAA expected from the
	background distribution by taking the max of lambda_bg, lamda_1, lamda_5,
	lambda_10.  It then computes a p-value based on the expected number of hops 
	= lamda * number of TTAAs in peak * number of hops in peak.  Finally, it 
	returns a frame that has 
	[Chr,Start,End,Center,Experiment Hops,
	Fraction Experiment,Background Hops,Fraction Background,Poisson pvalue]

	"""
	peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Center","Experiment Hops",
		"Fraction Experiment","TPH Experiment","Background Hops","Fraction Background","TPH Background","TPH Background subtracted","Lambda Type",
		"Lambda","Poisson pvalue"])

	
	experiment_gnashy_dict = {}
	experiment_dict_of_trees = {}
	total_experiment_hops = len(experiment_frame)
	background_gnashy_dict = {} 
	background_dict_of_trees = {}
	total_background_hops = len(background_frame)
	TTAA_frame_gbChr_dict = {} 
	TTAA_dict_of_trees = {}
	list_of_l_names = ["bg","1k","5k","10k"]

	#group by chromosome and populate interval tree with TTAA positions
	for name,group in TTAA_frame.groupby('Chr'):
		 
		TTAA_frame_gbChr_dict[name] = group
		TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
		#initialize tree
		TTAA_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in TTAA_frame_gbChr_dict[name].iterrows():	
			TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))
	#group by chromosome and populate interval tree with positions of background hops
	for name,group in background_frame.groupby('Chr'):
		background_gnashy_dict[name] = group
		background_gnashy_dict[name].index = background_gnashy_dict[name]["Start"]
		#initialize tree
		background_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in background_gnashy_dict[name].iterrows():	
			background_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 

	#group by chromosome and populate interval tree with positions of experiment hops
	for name,group in experiment_frame.groupby('Chr'):
		experiment_gnashy_dict[name] = group
		experiment_gnashy_dict[name].index = experiment_gnashy_dict[name]["Start"]
		#initialize tree
		experiment_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in experiment_gnashy_dict[name].iterrows():	
			experiment_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 

	#these will eventually be the columns in the peaks frame that will be returned.
	chr_list = []
	start_list = []
	end_list = []
	center_list = []
	num_exp_hops_list = []
	num_bg_hops_list = []
	frac_exp_list = []
	tph_exp_list = []
	frac_bg_list = []
	tph_bg_list = []
	tph_bgs_list = []
	lambda_type_list =[]
	lambda_list = []
	pvalue_list = []
	l = []
	
	#group experiment gnashyfile by chomosome
	for name,group in experiment_frame.groupby('Chr'):
		max_pos = max(group["End"])
		sig_start = 0
		sig_end = 0
		sig_flag = 0
		print max_pos
		print window_size
		print step_size
		for window_start in range(1,max_pos+window_size,step_size):
			overlap = experiment_dict_of_trees[name].find(window_start,window_start+window_size - 1)
			num_exp_hops = len(overlap)
			bg_overlap = background_dict_of_trees[name].find(window_start,window_start+window_size - 1)
			num_bg_hops = len(bg_overlap)

			#is this window significant?
			if compute_cumulative_poisson(num_exp_hops,num_bg_hops,total_experiment_hops,total_background_hops,pseudocounts) < pvalue_cutoff:
				#was last window significant?
				if sig_flag:
					#if so, extend end of windows
					sig_end = window_start+window_size-1
				else:
					#otherwise, define new start and end and set flag
					sig_start = window_start
					sig_end = window_start+window_size-1
					sig_flag = 1

			else:
				#current window not significant.  Was last window significant?
				if sig_flag:
					#add full sig window to the frame of peaks 
					
					#add chr, peak start, peak end
					chr_list.append(name) #add chr to frame
					start_list.append(sig_start) #add peak start to frame
					end_list.append(sig_end) #add peak end to frame
				
					#compute peak center and add to frame
					overlap = experiment_dict_of_trees[name].find(sig_start,sig_end)
					exp_hop_pos_list = [x.start for x in overlap]
					peak_center = np.median(exp_hop_pos_list)
					center_list.append(peak_center) #add peak center to frame

					#add number of experiment hops in peak to frame
					num_exp_hops = len(overlap)
					num_exp_hops_list.append(num_exp_hops)

					#add fraction of experiment hops in peak to frame
					frac_exp_list.append(float(num_exp_hops)/total_experiment_hops)
					tph_exp_list.append(float(num_exp_hops)*100000/total_experiment_hops)


					#add number of background hops in peak to frame
					bg_overlap = background_dict_of_trees[name].find(sig_start,sig_end)
					num_bg_hops = len(bg_overlap)
					num_bg_hops_list.append(num_bg_hops)

					frac_bg_list.append(float(num_bg_hops)/total_background_hops)
					tph_bg_list.append(float(num_bg_hops)*100000/total_background_hops)		
					#find lambda and compute significance of peak
					if total_background_hops >= total_experiment_hops: #scale bg hops down
						#compute lambda bg
						num_TTAAs = len(TTAA_dict_of_trees[name].find(sig_start,sig_end))
						lambda_bg = ((num_bg_hops*(float(total_experiment_hops)/total_background_hops))/max(num_TTAAs,1)) 


						#compute lambda 1k
						num_bg_hops_1k = len(background_dict_of_trees[name].find(peak_center-499,peak_center+500))
						num_TTAAs_1k = len(TTAA_dict_of_trees[name].find(peak_center-499,peak_center+500))
						lambda_1k = (num_bg_hops_1k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_1k,1))


						#compute lambda 5k
						num_bg_hops_5k = len(background_dict_of_trees[name].find(peak_center-2499,peak_center+2500))
						num_TTAAs_5k = len(TTAA_dict_of_trees[name].find(peak_center-2499,peak_center+2500))
						lambda_5k = (num_bg_hops_5k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_5k,1))


						#compute lambda 10k
						num_bg_hops_10k = len(background_dict_of_trees[name].find(peak_center-4999,peak_center+5000))
						num_TTAAs_10k = len(TTAA_dict_of_trees[name].find(peak_center-4999,peak_center+5000))
						lambda_10k = (num_bg_hops_10k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_10k,1))
						lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])


						#record type of lambda used
						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambda_type_list.append(list_of_l_names[index])
						#record lambda
						lambda_list.append(lambda_f)
						#compute pvalue and record it

						pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
						pvalue_list.append(pvalue)
						

						tph_bgs = float(num_exp_hops)*100000/total_experiment_hops - float(num_bg_hops)*100000/total_background_hops
						tph_bgs_list.append(tph_bgs)

						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambdatype = list_of_l_names[index]
						l = [pvalue,tph_bgs,lambda_f,lambdatype]

					else: #scale experiment hops down
						#compute lambda bg
						num_TTAAs = len(TTAA_dict_of_trees[name].find(sig_start,sig_end))
						lambda_bg = (float(num_bg_hops)/max(num_TTAAs,1)) 



						#compute lambda 1k
						num_bg_hops_1k = len(background_dict_of_trees[name].find(peak_center-499,peak_center+500))
						num_TTAAs_1k = len(TTAA_dict_of_trees[name].find(peak_center-499,peak_center+500))
						lambda_1k = (float(num_bg_hops_1k)/(max(num_TTAAs_1k,1)))



						#compute lambda 5k
						num_bg_hops_5k = len(background_dict_of_trees[name].find(peak_center-2499,peak_center+2500))
						num_TTAAs_5k = len(TTAA_dict_of_trees[name].find(peak_center-2499,peak_center+2500))
						lambda_5k = (float(num_bg_hops_5k)/(max(num_TTAAs_5k,1)))



						#compute lambda 10k
						num_bg_hops_10k = len(background_dict_of_trees[name].find(peak_center-4999,peak_center+5000))
						num_TTAAs_10k = len(TTAA_dict_of_trees[name].find(peak_center-4999,peak_center+5000))
						lambda_10k = (float(num_bg_hops_10k)/(max(num_TTAAs_10k,1)))
						lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])



						#record type of lambda used
						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambda_type_list.append(list_of_l_names[index])
						#record lambda
						lambda_list.append(lambda_f)
						#compute pvalue and record it
						pvalue = 1-scistat.poisson.cdf(((float(total_background_hops)/total_experiment_hops)*num_exp_hops+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
						pvalue_list.append(pvalue)

						tph_bgs = float(num_exp_hops)*100000/total_experiment_hops - float(num_bg_hops)*100000/total_background_hops
						tph_bgs_list.append(tph_bgs)

						index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
						lambdatype = list_of_l_names[index]
						l = [pvalue,tph_bgs,lambda_f,lambdatype]

					#number of hops that are a user-defined distance from peak center
					sig_flag = 0

					

				#else do nothing.
				
				
	#make frame from all of the lists
	peaks_frame["Chr"] = chr_list
	peaks_frame["Start"] = start_list
	peaks_frame["End"] = end_list
	peaks_frame["Center"] = center_list
	peaks_frame["Experiment Hops"] = num_exp_hops_list 
	peaks_frame["Fraction Experiment"] = frac_exp_list 
	peaks_frame["TPH Experiment"] = tph_exp_list
	peaks_frame["Background Hops"] = num_bg_hops_list 
	peaks_frame["Fraction Background"] = frac_bg_list
	peaks_frame["TPH Background"] = tph_bg_list
	peaks_frame["TPH Background subtracted"] = tph_bgs_list
	peaks_frame["Lambda Type"] = lambda_type_list
	peaks_frame["Lambda"] = lambda_list
	peaks_frame["Poisson pvalue"] = pvalue_list
	return peaks_frame


################################
#The background free peakcaller#
################################

def find_peaks_bf(experiment_frame,TTAA_frame,pvalue_cutoff = 1e-3,window_size = 1000,lam_win_size=100000,step_size = 500,pseudocounts = 0.2):
	"""This function is the same as find_peaks but it calls peaks without using
	a background distribution.  To do so, it scans through the genome
	with a window of window_size and step size of step_size and looks for regions that
	have signficantly more experiment hops in the window than expected from the number of
	experimenta hops in a window of size lam_win_size.
	It merges consecutively enriched windows and computes the center of the peak.  Next 
	it computes lambda, the number of insertions per TTAA expected from the experimental
	distribution where the expected number of hops are
	estimated from the number of hops in a window of size lam_win_size around the peak center.  It then computes
	a p-value based on the expected number of hops = lamda * number of TTAAs in peak * number of hops
	in peak.  Finally, it returns a frame that has Chr,Start,End,Center,Experiment Hops,
	Fraction Experiment Lambda Type, Lambda,Poisson pvalue

	"""
	peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Center","Experiment Hops",
		"Fraction Experiment","TPH Experiment","Lambda Type",
		"Lambda","Poisson pvalue"])

	
	experiment_gnashy_dict = {}
	experiment_dict_of_trees = {}
	total_experiment_hops = len(experiment_frame)
	TTAA_frame_gbChr_dict = {} 
	TTAA_dict_of_trees = {}
	list_of_l_names = lam_win_size

	#group by chromosome and populate interval tree with TTAA positions
	for name,group in TTAA_frame.groupby('Chr'):
		 
		TTAA_frame_gbChr_dict[name] = group
		TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
		#initialize tree
		TTAA_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in TTAA_frame_gbChr_dict[name].iterrows():	
			TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))

	#group by chromosome and populate interval tree with positions of experiment hops
	for name,group in experiment_frame.groupby('Chr'):
		experiment_gnashy_dict[name] = group
		experiment_gnashy_dict[name].index = experiment_gnashy_dict[name]["Start"]
		#initialize tree
		experiment_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in experiment_gnashy_dict[name].iterrows():	
			experiment_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 

	#these will eventually be the columns in the peaks frame that will be returned.
	chr_list = []
	start_list = []
	end_list = []
	center_list = []
	num_exp_hops_list = []
	frac_exp_list = []
	tph_exp_list = []
	lambda_type_list =[]
	lambda_list = []
	pvalue_list = []
	l = []
	
	#group experiment gnashyfile by chomosome
	for name,group in experiment_frame.groupby('Chr'):
		max_pos = max(group["End"])
		sig_start = 0
		sig_end = 0
		sig_flag = 0
		for window_start in range(lam_win_size/2,max_pos+window_size+(lam_win_size/2+1),step_size):
			overlap = experiment_dict_of_trees[name].find(window_start,window_start+window_size - 1)
			num_exp_hops = len(overlap)
			overlap_lam_win_size = experiment_dict_of_trees[name].find(window_start-(lam_win_size/2-1),window_start+window_size +(lam_win_size/2) - 1)
			num_lam_win_size_hops = len(overlap_lam_win_size)
			num_TTAAs_window = len(TTAA_dict_of_trees[name].find(window_start,window_start+window_size - 1))
			num_TTAAs_lam_win_size = len(TTAA_dict_of_trees[name].find(window_start-(lam_win_size/2-1),window_start+window_size +lam_win_size/2 - 1))

			lambda_lam_win_size = (num_lam_win_size_hops/(max(num_TTAAs_lam_win_size,1))) #expected number of hops per TTAA

			#is this window significant?
			pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_lam_win_size*max(num_TTAAs_window,1)+pseudocounts)
			if pvalue < pvalue_cutoff:
				#was last window significant?
				if sig_flag:
					#if so, extend end of windows
					sig_end = window_start+window_size-1
				else:
					#otherwise, define new start and end and set flag
					sig_start = window_start
					sig_end = window_start+window_size-1
					sig_flag = 1

			else:
				#current window not significant.  Was last window significant?
				if sig_flag:
					#add full sig window to the frame of peaks 
					
					#add chr, peak start, peak end
					chr_list.append(name) #add chr to frame
					start_list.append(sig_start) #add peak start to frame
					end_list.append(sig_end) #add peak end to frame
				
					#compute peak center and add to frame
					overlap = experiment_dict_of_trees[name].find(sig_start,sig_end)
					exp_hop_pos_list = [x.start for x in overlap]
					peak_center = np.median(exp_hop_pos_list)
					center_list.append(peak_center) #add peak center to frame

					#add number of experiment hops in peak to frame
					num_exp_hops = len(overlap)
					num_exp_hops_list.append(num_exp_hops)

					#add fraction of experiment hops in peak to frame
					frac_exp_list.append(float(num_exp_hops)/total_experiment_hops)
					tph_exp_list.append(float(num_exp_hops)*100000/total_experiment_hops)

					num_TTAAs_peak = len(TTAA_dict_of_trees[name].find(sig_start,sig_end))

					#compute lambda in lam_win_size
					num_exp_hops_lam_win_size = len(experiment_dict_of_trees[name].find(peak_center-(lam_win_size/2-1),peak_center+(lam_win_size/2)))
					num_TTAAs_lam_win_size = len(TTAA_dict_of_trees[name].find(peak_center-(lam_win_size/2-1),peak_center+(lam_win_size/2)))
					lambda_lam_win_size = float(num_exp_hops_lam_win_size)/(max(num_TTAAs_lam_win_size,1))


					lambda_f = lambda_lam_win_size


					#record type of lambda used
					lambda_type_list.append(list_of_l_names)
					#record lambda
					lambda_list.append(lambda_f)
					#compute pvalue and record it

					pvalue = 1-scistat.poisson.cdf((num_exp_hops+pseudocounts),lambda_f*max(num_TTAAs_peak,1)+pseudocounts)
					pvalue_list.append(pvalue)              
					#number of hops that are a user-defined distance from peak center
					sig_flag = 0
		
                                        
					lambdatype = lam_win_size
					l = [pvalue,float(num_exp_hops)*100000/total_experiment_hops,lambda_f,lambdatype]

				#else do nothing.
				
				
	#make frame from all of the lists
	peaks_frame["Chr"] = chr_list
	peaks_frame["Start"] = start_list
	peaks_frame["End"] = end_list
	peaks_frame["Center"] = center_list
	peaks_frame["Experiment Hops"] = num_exp_hops_list 
	peaks_frame["Fraction Experiment"] = frac_exp_list 
	peaks_frame["TPH Experiment"]= tph_exp_list
	peaks_frame["Lambda Type"] = lambda_type_list
	peaks_frame["Lambda"] = lambda_list
	peaks_frame["Poisson pvalue"] = pvalue_list
	return peaks_frame


##################
#Helper Functions#
##################

def compute_cumulative_poisson(exp_hops_region,bg_hops_region,total_exp_hops,total_bg_hops,pseudocounts):
	#usage
	#scistat.poisson.cdf(x,mu)
	#scales sample with more hops down to sample with less hops
	#tested 6/14/17
	if total_bg_hops >= total_exp_hops:
		return(1-scistat.poisson.cdf((exp_hops_region+pseudocounts),bg_hops_region * (float(total_exp_hops)/float(total_bg_hops)) + pseudocounts))
	else:
		return(1-scistat.poisson.cdf(((exp_hops_region *(float(total_bg_hops)/float(total_exp_hops)) )+pseudocounts),bg_hops_region + pseudocounts))

####################################
#Code for calling from command line#
####################################

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='cc_filter_reads.py')
	parser.add_argument('-e','--exp_file',help='experiment filename (full path)',required=True)
	parser.add_argument('-o','--output_file',help='output filename (full path)',required=True)
	parser.add_argument('-t','--TTAA_file',help='TTAA_filename',required=True)
	parser.add_argument('-a','--annotation_file',help='annotation filename (full path)',required=True)
	parser.add_argument('-b','--background_file',help='background filename (full path)',required=False,default=False)
	parser.add_argument('-pc','--pvaluecutoff',help='pvalue cutoff for significant peaks',required=False,default=1e-4)
	parser.add_argument('--peak_finder_pvalue',help='p value cutoff for identifying sig peaks',required=False,default=1e-3)
	parser.add_argument('--window',help='window size',required=False,default=1000)
	parser.add_argument('--lam_win_size',help='window size for lambda computation in background free peak-caller',required=False,default=100000)
	parser.add_argument('--step',help='step size',required=False,default=500)
	parser.add_argument('--pseudocounts',help='pseudocounts',required=False,default=0.2)

	args = parser.parse_args()
	
	if args.background_file:
		call_peaks_and_annotate(args.exp_file,args.background_file,args.output_file,args.TTAA_file,
			args.annotation_file,float(args.pvaluecutoff),float(args.peak_finder_pvalue),int(args.window),
			int(args.step),float(args.pseudocounts))
	else:
		call_peaks_and_annotate_bf(args.exp_file,args.output_file,args.TTAA_file,
			args.annotation_file,float(args.pvaluecutoff),float(args.peak_finder_pvalue),int(args.window),int(args.lam_win_size),
			int(args.step),float(args.pseudocounts))

