# CallingCard-seq pipelineV1 (pythonv2.7 + bash)
PiggyBac transposon genomic mapping

M Yen, Z Qi et al. (2018) Genomic targets of BAP1 in uveal melanoma identified by transposase mapping. BMC Medical Genomics 11 (1), 97

Z Qi et al. An optimized, broadly applicable piggyBac transposon induction system. (2017) Nucleic Acids Research 45 (7), e55-e55



I.	Overview

This toolbox is designed to perform two major functions: 

1.  To take calling card data in fastq format, map these data to the genome and create a CCF file that can be used for visualization, and 

2.  To call significant binding peaks from a CCF file and report this data in a SIG file.  

The toolbox has been tailored for the Washington University HTCF computing environment 
(https://htcf.wustl.edu/docs/) which uses the slurm queueing system 
(https://slurm.schedmd.com/tutorials.html).    

II.  	Quick Start

A.  Making a CCF file from Illumina Reads.

1. In your experiment directory, make three subdirectories named raw/, output_and_analysis/, and scripts/.  Put your unzipped fastq data in the raw directory.  Name the files combined_R1.fastq etc. (this naming is not required, but then you will not have to change any options from the example shells)

2. In the scripts directory, make a barcodes.txt file with the format: experiment name \t primer barcode \t transposon barcode 1, transposon barcode 2, ... transposon barcode M \n. There is an example of a barcodes.txt file here: /opt/apps/labs/rmlab/code/rmitra/mammalian_cc_tools/ccf_tools_v1.0/bash_scripts/

3. From the directory 

/opt/apps/labs/rmlab/code/rmitra/mammalian_cc_tools/ccf_tools_v1.0/bash_scripts/example_human_scripts

copy map_reads.sh, split_reads.sh and make_ccffile.sh to your scripts/ folder in your experiment directory.  Open split_reads to make sure the 3' end transposon sequence matches your experiment and double check the file names.  Open map_reads.sh to make sure the mapping parameters are correct.  [For the analysis of AAV experiments, copy map_reads_aav.sh and split_reads_aav.sh]

4. ÒcdÓ into your scripts/ folder and type ÒinteractiveÓ.  Then type Òmodule load ccf_toolsÓ.  Then type Òpython $CCF_TOOLS/run_all.pyÓ [For the analysis of AAV experiments, use run_all_aav.py].  This will split the reads by barcode into different fasta files, batch them out to different nodes for mapping by bowtie2, and create a CCF file for each experiment in the barcodes.txt file.  The CCF files (with .ccf extension) are located in the output_and_analysis directory.   










usage:
python cc_filter_reads.py -r1 <read1 file> -r2 <read2 file> 
-i <index file> -b<barcode file> -o <output path>
--hammp <hamming distance for primer barcode>
--hammt <hamming distance for transposon barcode>

required fields:
    -r1 read 1 filename (full path)
    -i index filename (full path)

not required
    -r2 read 2 filename (full path)
    -b barcode file = ../raw/barcodes.txt
    -o output path = ../output_and_analysis
    -hp 0
    -tp 0
