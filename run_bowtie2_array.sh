#!/bin/bash
#SBATCH -n 8  
#SBATCH -N 1 
#SBATCH --mem=30000
#SBATCH -o run_bowtie.out#standard out goes here
#SBATCH -e run_bowtie.err # standard error goes here
#SBATCH -J BOW
#SBATCH --array=1-30   # Adding this means that you now have an enivornment variable called SLURM_ARRAY_TASK_ID

module load bowtie2/2.1.0
module load samtools
export BOWTIE2_INDEXES='/scratch/rmlab/ref/bowtie2_indexes/'

filename=$( sed -n ${SLURM_ARRAY_TASK_ID}p ../output_and_analysis/cc_filelist.txt)  # extracts line $SLURM_ARRAY_TASK_ID from file

srun python map_reads.py -f ${filename} -g hg18 --ttrim 100 --ftrim 37 -p False -o ../output_and_analysis

