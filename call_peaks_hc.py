#!/usr/bin/env python
import numpy as np
import pandas as pd
import scipy.stats as scistat
from bx.intervals.intersection import Intersecter, Interval
import pybedtools.bedtool as bt 
import argparse
from scipy.cluster.hierarchy import dendrogram,linkage,fcluster
import re
import annotate
import sys

"""This module calls peaks from callng card data in ccf format.  The main peak
calling function is passed a experiment frame, a background frame, and an TTAA_frame, 
all in ccf format.  It then builds interval trees containing
all of the background and experiment hops and all of the TTAAs.  Next, 
it clusters the hops to call peaks using a few tricks described in the function.
Next it computes lambda, the number of insertions per TTAA expected from the 
background distribution by taking the max of lambda_bg, lamda_1, lamda_5, lambda_10.
It then computes a p-value based on the expected number of hops = lamda * number of TTAAs 
in peak * number of hops in peak.  Finally, it returns a frame that has 
Chr,Start,End,Center,Experiment Hops,Fraction Experiment,Background Hops,
Fraction Background,Poisson pvalue

There is also a background free version of this algorithm that computes significance 
based on the number of hops in the neighboring region.

"""

################
##Master Loop ##
################

def cluster_and_annotate(expfile,bgfile,outputfile,TTAA_file = '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt',
	annotation_file = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed',pvalue_cutoff=1e-4,
	cluster_size = 1000,pseudocounts = 0.2,macsstyle = True):

	"""This function is the wrapper that calls the various subfunctions"""

	expframe = pd.read_csv(expfile,delimiter="\t",header=None)
	expframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	bgframe = pd.read_csv(bgfile,delimiter="\t",header=None)
	bgframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	TTAAframe = pd.read_csv(TTAA_file,delimiter="\t",header=None)
	TTAAframe.columns = ["Chr","Start","End"]

	total_experiment_hops = len(expframe)
	peaks_frame = cluster_peaks_frame(expframe,cluster_size)
	peaks_frame = add_pvalues_to_peaks_frame_macs(peaks_frame,bgframe,
		TTAAframe,total_experiment_hops,pseudocounts,macsstyle)
	peaks_frame['TPH Background subtracted'] = peaks_frame['TPH Experiment']-peaks_frame['TPH Background']
	peaks_frame = annotate.annotate_peaks_frame(peaks_frame,annotation_file)
	peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
	peaks_frame = peaks_frame.sort_values(["Poisson pvalue"])
	cols = ["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment",
		"Background Hops","Fraction Background","TPH Background","TPH Background subtracted","Poisson pvalue","Lambda","Lambda Type","Feature 1 sName","Feature 1 Name",
		"Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
		"Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End","Feature 2 Strand",
		"Feature 2 Distance"]
	peaks_frame = peaks_frame[cols] #reorder columns
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

def cluster_and_annotate_bf(expfile,outputfile,TTAA_file = '/scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10.txt',
	annotation_file = '/scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed',pvalue_cutoff=1e-4,
	cluster_size = 1000,window_size=50000,pseudocounts = 0.2,macsstyle = True):

	"""This function is the wrapper that calls the various subfunctions"""

	expframe = pd.read_csv(expfile,delimiter="\t",header=None)
	expframe.columns = ["Chr","Start","End","Reads","Strand","Barcode"]

	TTAAframe = pd.read_csv(TTAA_file,delimiter="\t",header=None)
	TTAAframe.columns = ["Chr","Start","End"]

	total_experiment_hops = len(expframe)
	peaks_frame = cluster_peaks_frame(expframe,cluster_size)
	peaks_frame = add_pvalues_to_peaks_frame_macs_bf(peaks_frame,
			expframe,TTAAframe,window_size,pseudocounts,macsstyle)
	peaks_frame = annotate.annotate_peaks_frame(peaks_frame,annotation_file)
	peaks_frame = peaks_frame[peaks_frame["Poisson pvalue"] <= pvalue_cutoff]
	peaks_frame = peaks_frame.sort_values(["Poisson pvalue"])
	cols = ["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment",
		"Poisson pvalue","Lambda","Lambda Type","Feature 1 sName","Feature 1 Name",
		"Feature 1 Start","Feature 1 End","Feature 1 Strand","Feature 1 Distance", 
		"Feature 2 sName","Feature 2 Name","Feature 2 Start","Feature 2 End","Feature 2 Strand",
		"Feature 2 Distance"]
	peaks_frame = peaks_frame[cols] #reorder columns
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


#########################################
#Call Peaks via Hierarchical Clustering #
#########################################

def cluster_peaks_frame(exp_frame,distance_cutoff):
	"""This function takes a gnashy dataframe as input and clusters the insertions
	by hierarchical clustering.  It uses euclidean distance as the metric and seeks to minimize
	the average pairwise distance of cluster members.  The distance_cutoff parameter is used to define
	the different clusters.  Because hierarchical clustering is computationally inefficient, it breaks
	the chromosomes up into different pieces and clusters the individual pieces.  To ensure that the breakpoints
	don't break up potentia cluster, it searches for regions of the chromosome with a large distance between
	insertions.  If it cannot find a distance that is greater than the distance cutoff, then it prints a 
	warning.  The program outputs a frame of clusters with the following columns: 
	[Chr,Start,End,Center,Experiment Hops,Fraction Experiment]"""

	print "hello"
	MAX_ARRAY_SIZE = 5000  #maximum number of insertions to consider for cluster
	MIN_JUMP = 1000 #minimum number of insertions to consider for cluster

	peaks_frame = pd.DataFrame(columns = ["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment"])

	grouped = exp_frame.groupby(['Chr']) #group insertions by chromosome
	
	for chro in grouped.groups:
		print "Analyzing chromosome "+str(chro)
		full_hop_list = list(grouped.get_group(chro)['Start'])
		#find gap to split chromosome for clustering
		next_index = 0
		while next_index < len(full_hop_list):
			if (len(full_hop_list)-next_index) < MAX_ARRAY_SIZE:
				hop_list = full_hop_list[next_index:len(full_hop_list)]
				hop_list = [[x] for x in hop_list]
				old_index = next_index
				next_index = len(full_hop_list)
			else:
				gap_list = full_hop_list[next_index+MIN_JUMP:next_index+MAX_ARRAY_SIZE]
				difference_list = [j-i for i, j in zip(gap_list[:-1], gap_list[1:])]
				max_dist = max(difference_list)
				gap_index = difference_list.index(max(difference_list))+1
				hop_list = full_hop_list[next_index:next_index+MIN_JUMP+gap_index]
				hop_list = [[x] for x in hop_list]
				old_index = next_index
				next_index = next_index+MIN_JUMP+gap_index
				if max_dist < distance_cutoff:
					print "Warning. Chromosome "+str(chro)+" "+hop_list[0]+" jump distance of "+\
					max_dist+" is less than distance_cutoff of "+distance_cutoff+"."
				#print old_index, next_index-1, max_dist
			#cluster split chromosome
			Z = linkage(hop_list,'average','euclidean')
			#find clusters
			try:
				clusters = fcluster(Z,distance_cutoff,criterion='distance')
			except ValueError:
				pass
			#record clusters
			pos_clus_frame = pd.DataFrame(columns=['Hop Position','Cluster'])
			pos_clus_frame['Hop Position'] = hop_list
			try:
				pos_clus_frame['Cluster'] = clusters
			except ValueError:
				pass
			grouped_hops = pos_clus_frame.groupby(['Cluster'])
			for cluster in grouped_hops.groups:
				hops_in_cluster = list(grouped_hops.get_group(cluster)['Hop Position'])
				chromosome = chro
				start = min(hops_in_cluster)[0]
				end = max(hops_in_cluster)[0]+3
				center = np.median(hops_in_cluster)
				experiment_hops = len(hops_in_cluster) 
				fraction_experiment = float(experiment_hops)/len(exp_frame)
				tph_experiment = float(experiment_hops)*100000/len(exp_frame)
				peaks_frame = peaks_frame.append(pd.DataFrame({"Chr":[chromosome],"Start":[start],
					"End":[end],"Center":[center],"Experiment Hops":[experiment_hops],
					"Fraction Experiment":[fraction_experiment],"TPH Experiment":[tph_experiment]}),ignore_index = True)
	#sort cluster frame by chr then position
	chrlist = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10','chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20','chr21','chr22','chrX','chrY']
	peaks_frame["Chr"] = pd.Categorical(peaks_frame["Chr"],chrlist)
	peaks_frame = peaks_frame.sort_values(["Chr","Start"])
	peaks_frame = peaks_frame[["Chr","Start","End","Center","Experiment Hops","Fraction Experiment","TPH Experiment"]]
	peaks_frame = peaks_frame.reset_index()
	peaks_frame[["Start","End","Center","Experiment Hops"]] = peaks_frame[["Start","End","Center","Experiment Hops"]].astype(int)
	del peaks_frame["index"]
	peaks_frame.to_csv("clustered_hops.txt",sep="\t",index=False)
	return peaks_frame

############################
#Add pvalues to peaks_frame#
############################


def add_pvalues_to_peaks_frame_macs(peaks_frame,background_peaks_frame,TTAA_frame,total_experiment_hops,pseudocounts = 0.2,macs_pvalue=True):
	background_gnashy_dict = {} 
	background_dict_of_trees = {}
	total_background_hops = len(background_peaks_frame)
	TTAA_frame_gbChr_dict = {} 
	TTAA_dict_of_trees = {}
	list_of_l_names = ["bg","1k","5k","10k"]

	print "working on it..."
	#make interval tree for TTAAs
	for name,group in TTAA_frame.groupby('Chr'): 
		TTAA_frame_gbChr_dict[name] = group
		TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
		#initialize tree
		TTAA_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in TTAA_frame_gbChr_dict[name].iterrows(): 
			TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))
	print "Made TTAA interval tree"
	#group by chromosome and populate interval tree with positions of background hops
	for name,group in background_peaks_frame.groupby('Chr'):
		background_gnashy_dict[name] = group
		background_gnashy_dict[name].index = background_gnashy_dict[name]["Start"]
		#initialize tree
		background_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in background_gnashy_dict[name].iterrows():    
			background_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) #may need to be one.  Try with zero
	print "Made Background tree"
	#go through cluster frame and compute pvalues 
	num_bg_hops_list =[]
	frac_bg_list = []
	tph_bg_list = []
	lambda_type_list =[]
	lambda_list = []
	pvalue_list = []
	for idx,row in peaks_frame.iterrows():
		#add number of background hops in cluster to frame
		cluster_center = row["Center"]
		bg_overlap = background_dict_of_trees[row["Chr"]].find(row["Start"],row["End"])
		num_bg_hops = len(bg_overlap)
		num_bg_hops_list.append(num_bg_hops)
		#add fraction background 
		frac_bg_list.append(float(num_bg_hops)/total_background_hops)       
		tph_bg_list.append(float(num_bg_hops)*100000/total_background_hops)
		#find lambda and compute significance of cluster

		#replace num_exp hops with row["Experiment Hops"]

		if total_background_hops >= total_experiment_hops: #scale bg hops down
			#compute lambda bg
			num_TTAAs = len(TTAA_dict_of_trees[row["Chr"]].find(row["Start"],row["End"]))
			lambda_bg = ((num_bg_hops*(float(total_experiment_hops)/total_background_hops))/max(num_TTAAs,1)) 

			#compute lambda 1k
			num_bg_hops_1k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-499,cluster_center+500))
			num_TTAAs_1k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-499,cluster_center+500))
			lambda_1k = (num_bg_hops_1k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_1k,1))


			#compute lambda 5k
			num_bg_hops_5k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-2499,cluster_center+2500))
			num_TTAAs_5k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-2499,cluster_center+2500))
			lambda_5k = (num_bg_hops_5k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_5k,1))


			#compute lambda 10k
			num_bg_hops_10k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-4999,cluster_center+5000))
			num_TTAAs_10k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-4999,cluster_center+5000))
			lambda_10k = (num_bg_hops_10k*(float(total_experiment_hops)/total_background_hops))/(max(num_TTAAs_10k,1))
			if macs_pvalue:
				lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])
				index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
				lambda_type_list.append(list_of_l_names[index])
			else:
				lambda_f = lambda_bg
				lambda_type_list.append(list_of_l_names[0])
			lambda_list.append(lambda_f)
			#compute pvalue and record it
			pvalue = 1-scistat.poisson.cdf((row["Experiment Hops"]+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
			pvalue_list.append(pvalue)

		else: #scale experiment hops down
			#compute lambda bg
			num_TTAAs = len(TTAA_dict_of_trees[row["Chr"]].find(row["Start"],row["End"]))
			lambda_bg = (float(num_bg_hops)/max(num_TTAAs,1)) 

			#compute lambda 1k
			num_bg_hops_1k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-499,cluster_center+500))
			num_TTAAs_1k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-499,cluster_center+500))
			lambda_1k = (float(num_bg_hops_1k)/(max(num_TTAAs_1k,1)))

			#compute lambda 5k
			num_bg_hops_5k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-2499,cluster_center+2500))
			num_TTAAs_5k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-2499,cluster_center+2500))
			lambda_5k = (float(num_bg_hops_5k)/(max(num_TTAAs_5k,1)))

			#compute lambda 10k
			num_bg_hops_10k = len(background_dict_of_trees[row["Chr"]].find(cluster_center-4999,cluster_center+5000))
			num_TTAAs_10k = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-4999,cluster_center+5000))
			lambda_10k = (float(num_bg_hops_10k)/(max(num_TTAAs_10k,1)))
			
			if macs_pvalue:
				lambda_f = max([lambda_bg,lambda_1k,lambda_5k,lambda_10k])
				index = [lambda_bg,lambda_1k,lambda_5k,lambda_10k].index(max([lambda_bg,lambda_1k,lambda_5k,lambda_10k]))
				lambda_type_list.append(list_of_l_names[index])
			else:
				lambda_f = lambda_bg
				lambda_type_list.append(list_of_l_names[0])
			lambda_list.append(lambda_f)
			#compute pvalue and record it
			pvalue = 1-scistat.poisson.cdf(((float(total_background_hops)/total_experiment_hops)*row["Experiment Hops"]+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
			pvalue_list.append(pvalue)
						 
	#make frame from all of the lists 
	peaks_frame["Background Hops"] = num_bg_hops_list 
	peaks_frame["Fraction Background"] = frac_bg_list
	peaks_frame["TPH Background"] = tph_bg_list
	peaks_frame["Lambda Type"] = lambda_type_list
	peaks_frame["Lambda"] = lambda_list
	peaks_frame["Poisson pvalue"] = pvalue_list
	return peaks_frame

######################################################
#Add pvalues to peaks_frame background free algorithm#
######################################################


def add_pvalues_to_peaks_frame_macs_bf(peaks_frame,experiment_peaks_frame,TTAA_frame,lam_win_size,pseudocounts = 0.2,macs_pvalue=True):
	print "lab specific hohoho"
	experiment_gnashy_dict = {}
	experiment_dict_of_trees = {}
	TTAA_frame_gbChr_dict = {} 
	TTAA_dict_of_trees = {}
	list_of_l_names = [lam_win_size]
	print "Making interval tree for experiment hops..."
	for name,group in experiment_peaks_frame.groupby('Chr'):
		experiment_gnashy_dict[name] = group
		experiment_gnashy_dict[name].index = experiment_gnashy_dict[name]["Start"]
		#initialize tree
		experiment_dict_of_trees[name] = Intersecter()
		#populate tree with position as interval
		for idx, row in experiment_gnashy_dict[name].iterrows():
			experiment_dict_of_trees[name].add_interval(Interval(int(idx),int(idx)+3)) 
	print "Making interval tree for TTAAs..."
	#make interval tree for TTAAs
	for name,group in TTAA_frame.groupby('Chr'): 
		TTAA_frame_gbChr_dict[name] = group
		TTAA_frame_gbChr_dict[name].index = TTAA_frame_gbChr_dict[name]["Start"]
		#initialize tree
		TTAA_dict_of_trees[name] = Intersecter() 
		#populate tree with position as interval
		for idx, row in TTAA_frame_gbChr_dict[name].iterrows(): 
			TTAA_dict_of_trees[name].add_interval(Interval(int(idx),int(idx+3)))
	#go through cluster frame and compute pvalues 
	lambda_type_list =[]
	lambda_list = []
	pvalue_list = []
	for idx,row in peaks_frame.iterrows():
		#add number of background hops in cluster to frame
		cluster_center = row["Center"]
		#find lambda and compute significance of cluster
		num_TTAAs = len(TTAA_dict_of_trees[row["Chr"]].find(row["Start"],row["End"]))
		#compute lambda for window size
		num_exp_hops_lam_win_size = len(experiment_dict_of_trees[row["Chr"]].find(cluster_center-(lam_win_size/2 - 1),cluster_center+lam_win_size/2))
		num_TTAAs_lam_win_size = len(TTAA_dict_of_trees[row["Chr"]].find(cluster_center-(lam_win_size/2 - 1),cluster_center+lam_win_size/2))
		lambda_win_size = float(num_exp_hops_lam_win_size)/(max(num_TTAAs_lam_win_size,1))

		lambda_f = lambda_win_size
		lambda_type_list.append(lam_win_size)
		lambda_list.append(lambda_f)
		#compute pvalue and record it
		pvalue = 1-scistat.poisson.cdf((row["Experiment Hops"]+pseudocounts),lambda_f*max(num_TTAAs,1)+pseudocounts)
		pvalue_list.append(pvalue)

						 
	#make frame from all of the lists 
	peaks_frame["Lambda Type"] = lambda_type_list
	peaks_frame["Lambda"] = lambda_list
	peaks_frame["Poisson pvalue"] = pvalue_list
	return peaks_frame

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
	parser.add_argument('--cl_size',help='cluster size',required=False,default=1000)
	parser.add_argument('--lam_win_size',help = 'window size',required=False,default=50000)
	parser.add_argument('--pseudocounts',help='pseudocounts',required=False,default=0.2)
	parser.add_argument('--macsstyle',help='macsstyle sig calling',required=False,default=False)

	args = parser.parse_args()
        if args.macsstyle:
            if args.macsstyle in ["True","true","T",'1']:
                args.macsstyle = True
            else:
                args.macsstyle = False
        print "Hello"
	print "args.macsstyle="+str(args.macsstyle)
        if args.background_file:
		cluster_and_annotate(args.exp_file,args.background_file,args.output_file,args.TTAA_file,
			args.annotation_file,float(args.pvaluecutoff),int(args.cl_size),float(args.pseudocounts),bool(args.macsstyle))
	else:
		cluster_and_annotate_bf(args.exp_file,args.output_file,args.TTAA_file,
			args.annotation_file,float(args.pvaluecutoff),int(args.cl_size),int(args.lam_win_size),float(args.pseudocounts),bool(args.macsstyle))

