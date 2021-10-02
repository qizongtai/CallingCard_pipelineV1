#!/bin/bash
#SBATCH -n 1  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o call_peaks.out#standard out goes here
#SBATCH -e call_peaks.err # standard error goes here
#SBATCH -J PEAK

module load ccf_tools


python $CCF_TOOLS/call_peaks_macs.py -e yourfile.ccf -o outputname.sig \
	-t ~/ref/calling_card_ref/mouse/TTAA_mm10_ccf.txt -a ~/ref/calling_card_ref/mouse/refGene.mm10.Sorted.bed \
	-b yourbackground.ccf -pc 0.001 --peak_finder_pvalue 0.001 --window 500 --step 250 --pseudocounts 0.2

