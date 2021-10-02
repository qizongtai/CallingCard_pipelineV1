"""
annotate.py
This module annotates peaks frames.
""" 

import pandas as pd
import pybedtools.bedtool as bt 
import argparse
import re
import random

def peaks_frame_to_bed(peaks_frame,bedfilename):
	bed_frame = peaks_frame[["Chr","Start","End"]].copy()
	bed_frame.loc[:,"Start"] = bed_frame["Start"]-1 #start coords of bed files are 0 indexed while ends are 1 indexed
	bed_frame.to_csv(bedfilename,sep = "\t",header = None,index=None)

def make_peaksbed(peaks_frame):
	"""This function converts peaks from to a bed file that can be used to
	display the peak locations on the EPCC browser.  The minimum column 
	requirement for the peaks frame is [Chr,Start,End].  This function
	will return a bed frame with the following columns:  
	[Chr in mm10 or hg19 format,Start,End,name,score]"""

	#Do I need to subtract 1 from start??
	
	bed_frame = peaks_frame[["Chr","Start","End"]].copy()
	score_list = [1000]*len(peaks_frame)
	bed_frame["Score"] = score_list
	chrlist = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	bed_frame["Chr"] = pd.Categorical(bed_frame["Chr"],chrlist)
	bed_frame = bed_frame.sort_values(['Chr','Start','End'])
	return bed_frame

#def merge_bedfile_and_peaks_frame():

def annotate_peaks_frame(peaks_frame,refGene_filename = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed'):
	"This function annotates peaks using pybedtools"

	#convert peaks frame to sorted bed
	temp_filename = "temp_peaks_"+str(random.randint(1,1000))+".bed"
	peaks_frame_to_bed(peaks_frame,temp_filename)
	peaks_bed = bt.BedTool(temp_filename)
	peaks_bed = peaks_bed.sort()	
	peaks_bed = peaks_bed.closest(refGene_filename,D="ref",t="first",k=2)
	peaks_bed.saveas(temp_filename)
	#read in gene_annotation_bedfilename
	temp_annotated_peaks = pd.read_csv(temp_filename,delimiter = "\t",header=None)
	temp_annotated_peaks.columns = ["Chr","Start","End","Feature Chr","Feature Start", "Feature End","Feature sName","Feature Name","Strand","Distance"]
	temp_annotated_peaks.loc[:,"Start"] = temp_annotated_peaks["Start"]+1 #convert start coords back to 1 indexed coordinates
	index_list = [(x,y,z) for x,y,z in zip(temp_annotated_peaks["Chr"],temp_annotated_peaks["Start"],temp_annotated_peaks["End"])]

	temp_annotated_peaks.index = pd.MultiIndex.from_tuples(index_list)
	
	index_list = [(x,y,z) for x,y,z in zip(peaks_frame["Chr"],peaks_frame["Start"],peaks_frame["End"])]
	peaks_frame.index = pd.MultiIndex.from_tuples(index_list)
	peaks_frame = peaks_frame.sortlevel(0,axis=1)
	if "Background Hops" in peaks_frame.columns:
		temp_list = ["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment",
		"Background Hops","Fraction Background","TPH Background","TPH Background subtracted","Poisson pvalue","Lambda","Lambda Type","Feature 1 sName","Feature 1 Name",
		"Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
		"Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End","Feature 2 Strand",
		"Feature 2 Distance"]
	else:
		temp_list = ["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment",
		"Poisson pvalue","Lambda","Lambda Type","Feature 1 sName","Feature 1 Name",
		"Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
		"Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End","Feature 2 Strand",
		"Feature 2 Distance"]
	peaks_frame = peaks_frame.reindex(columns=temp_list,fill_value=0)

	for idx,row in temp_annotated_peaks.iterrows():
		if not(peaks_frame.loc[idx,"Feature 1 sName"]):
			peaks_frame.loc[idx,["Feature 1 sName","Feature 1 Name","Feature 1 Start","Feature 1 End","Feature 1 Strand",
			"Feature 1 Distance"]] = list(row[["Feature sName","Feature Name","Feature Start", "Feature End","Strand","Distance"]])	
		else:
			peaks_frame.loc[idx,["Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End","Feature 2 Strand",
			"Feature 2 Distance"]] = list(row[["Feature sName","Feature Name","Feature Start", "Feature End","Strand","Distance"]])
	return peaks_frame		
	
