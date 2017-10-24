#!/bin/bash
#SBATCH -n 8  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o run_make_gnashyfile.out#standard out goes here
#SBATCH -e run_make_gnashyfile.err # standard error goes here
#SBATCH -J CFR

module load pysam
module load pandas

python make_gnashyfile.py -g mm
        


