#!/bin/bash
#SBATCH -n 1  
#SBATCH -N 1 
#SBATCH --mem=20000
#SBATCH -o map_reads.out#standard out goes here
#SBATCH -e map_reads.err # standard error goes here
#SBATCH -J MAP

module load ccf_tools

filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p ../output_and_analysis/cc_filelist.txt)  # extracts line $SLURM_ARRAY_TASK_ID from file

srun python $CCF_TOOLS/map_reads.py -f ${filename} -g hg38 --ttrim 23 --ftrim 37 -p True -q 10 -o ../output_and_analysis

