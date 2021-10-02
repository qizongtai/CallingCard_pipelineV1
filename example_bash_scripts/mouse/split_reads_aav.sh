#!/bin/bash
#SBATCH -n 1  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o split_reads_aav.out#standard out goes here
#SBATCH -e split_reads_aav.err # standard error goes here
#SBATCH -J SPLIT

module load ccf_tools


$CCF_TOOLS/split_reads_aav.py -r1 ../raw/combined_R1.fastq -r2 ../raw/combined_R2.fastq \
	-b ../scripts/barcodes.txt --p3p GCGTCAATTTTACGCAGACTATCTTTCTAGGG \
	--hammp 0 --hammt 0 -o ../output_and_analysis
#The weird format for the path to barcodes.txt is required. Long story... just do it that way and I'll fix it eventually.
