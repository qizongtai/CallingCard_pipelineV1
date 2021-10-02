#!/bin/bash
#SBATCH -n 1  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o split_reads.out#standard out goes here
#SBATCH -e split_reads.err # standard error goes here
#SBATCH -J SPLIT

module load ccf_tools


$CCF_TOOLS/split_reads.py -r1 ../raw/combined_R1.fastq -r2 ../raw/combined_R2.fastq \
	-i ../raw/combined_I1.fastq -b ../scripts/barcodes.txt --p3p GCGTCAATTTTACGCAGACTATCTTTCTAGGG \
	--hammp 0 --hammt 0 -o ../output_and_analysis
