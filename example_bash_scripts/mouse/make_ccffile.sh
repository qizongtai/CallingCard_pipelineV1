#!/bin/bash
#SBATCH -n 1
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o make_ccffile.out#standard out goes here
#SBATCH -e make_ccffile.err # standard error goes here
#SBATCH -J CCF

module load ccf_tools

python $CCF_TOOLS/make_ccffile.py -g mm10 -b "../scripts/barcodes.txt" -p "../output_and_analysis"
#The weird format for barcodes.txt path is required... long story, just do it that way.        


