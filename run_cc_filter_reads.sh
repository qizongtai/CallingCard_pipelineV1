#!/bin/bash
#SBATCH -n 8  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o run_cc_filter_reads.out#standard out goes here
#SBATCH -e run_cc_filter_reads.err # standard error goes here
#SBATCH -J CFR

module load biopython


python cc_filter_reads.py -r1 ../raw/combined_R1.fastq -r2 ../raw/combined_R2.fastq \
	-i ../raw/combined_I1.fastq -b ../raw/barcodes.txt --hammp 0 --hammt 0 -o ../output_and_analysis

