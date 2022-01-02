# CallingCard-seq pipelineV1 (python2.7 + bash)

This pipeline was used in the following papers

- Genomic targets of BAP1 in uveal melanoma identified by transposase mapping. (2018) M Yen, Z Qi et al. BMC Medical Genomics 11 (1), 97 

  URL: https://bmcmedgenomics.biomedcentral.com/articles/10.1186/s12920-018-0424-0


- An optimized, broadly applicable piggyBac transposon induction system. (2017) Z Qi et al. Nucleic Acids Research 45 (7), e55-e55 

  URL: https://pubmed.ncbi.nlm.nih.gov/28082389/

-----
## Disclaimer
The pipeline has been tailored for the Washington University HTCF computing environment (https://htcf.wustl.edu/docs/) which uses the slurm queueing system (https://slurm.schedmd.com/tutorials.html). No guarantees are made about other systems, setups, configurations, etc. 

## Overview
This toolbox is designed to perform two major functions: 

-  To take calling card data in fastq format, map these data to the genome and create a CalingCardFormat (CCF) file that can be used for visualization, and 

-  To call significant binding peaks from a CCF file and report this data in a SIG file.  

The toolbox has been tailored for the Washington University HTCF computing environment 
(https://htcf.wustl.edu/docs/) which uses the slurm queueing system 
(https://slurm.schedmd.com/tutorials.html).    

## Quick Start

### A.  Making a CCF file from Illumina Reads.

1. In your experiment directory, make three subdirectories named raw/, output_and_analysis/, and scripts/.  Put your unzipped fastq data in the raw directory.  Name the files combined_R1.fastq etc. (this naming is not required, but then you will not have to change any options from the example shells)

2. In the scripts directory, make a barcodes.txt file with the format: experiment name \t primer barcode \t transposon barcode 1, transposon barcode 2, ... transposon barcode M \n. There is an example of a barcodes.txt file here: CallingCard_pipelineV1/example_bash_scripts/

3. From the bash shells directory copy map_reads.sh, split_reads.sh and make_ccffile.sh to the scripts/ directory. Open split_reads to make sure the 3' end transposon sequence matches your experiment and double check the file names. Open map_reads.sh to make sure the apping parameters are correct. [For the analysis of AAV experiments, copy map_reads_aav.sh and split_reads_aav.sh]

4. In the scripts directory type interactive. Then type module load ccf_tools. Then type “python $CCF_TOOLS/run_all.py” [For the analysis of AAV experiments, use run_all_aav.py]. This will split the reads by barcode into different fasta files, batch them out to different nodes for mapping by bowtie2, and create a CCF file for each experiment in the barcodes.txt file. The CCF files (with .ccf extension) are located in the output_and_analysis directory.


### B. Calling Significant Peaks from a CCF file.

1. Copy the template shell call_peaks.sh to your scripts directory. 

2. Modify the -e and -b parameters with the file names of your experiment and background file and the -o with the output parameter.

3. Modify the -t and -a paramters as follows (the TTAA locations file and the annotation file):

-  hg38:
   -t /scratch/ref/rmlab/calling_card_ref/human/TTAA_hg38_ccf.txt
   -a /scratch/ref/rmlab/calling_card_ref/human/refGene.hg38.Sorted.bed
-  hg19:
   -t /scratch/ref/rmlab/calling_card_ref/human/TTAA_hg19_ccf.txt
   -a /scratch/ref/rmlab/calling_card_ref/human/refGene.hg38.Sorted.bed
-  mm10:
   -t /scratch/ref/rmlab/calling_card_ref/mouse/TTAA_mm10_ccf.txt
   -a /scratch/ref/rmlab/calling_card_ref/mouse/refGene.mm10.Sorted.bed

4. Run the shell by typing: sbatch call_peaks.sh

5. For each experiment, two files will be created. They will end in .sig (significant peaks) and .peaks.bed (a bedgraph file that you can use to visualize peak calls on the epigenome browser at Washington Univerisity)


## The Different Types of Calling Card Data
This toolkit is designed to map raw Calling Card data for three different protocols (In addition to calling significant peaks for any .ccf file). These are the "V2" DNA calling card protocol, the AAV calling card protocol and the bulk SRT protocol. What matters for the purposes of mapping the reads is what sequences are collected in
Read1, Read 2 and Index 1. The software expects the following sequences:
### For V2 and SRT:
-  Read 1: has a 3-5 bp primer barcode followed by a sequence from the 5' TR of the piggyBac transposon. This sequence is typically TTTACGCAGACTATCTTTCTAGGG or GCGTCAATTTTACGCAGACTATCTTTCTAGGG. Then comes the TTAA and then genomic sequence.
-  Read 2: reads through the restriction enzyme bank (at most 11bp) and then into the genome. (For SRT bulk this is genome sequence right away, but for most purposes the same code can be used to analyze both experiments) 
-  Index 1: reads the transposon barcode (For SRTs this is the index read which serves as a "transposon barcode")

### For AAV:
-  Read1: 5bp primer barcode followed by GCGTCAATTTTACGCAGACTATCTTTCTAGGG TTAA then genomic sequence
-  Read2: 6bp transposon barcode followed by bank of TaqI,MspI, Cst6 restriction sites followed by genomic sequence. All samples should have the same group of transposon barcodes.

## Example Scripts

From the directory: CallingCard_pipelineV1/exmpale_bash_scripts/human

copy map_reads.sh, split_reads.sh and make_ccffile.sh to your scripts/ folder in your experiment directory.  Open split_reads to make sure the 3' end transposon sequence matches your experiment and double check the file names.  Open map_reads.sh to make sure the mapping parameters are correct.  [For the analysis of AAV experiments, copy map_reads_aav.sh and split_reads_aav.sh]

`cd` into your scripts/ folder and type `interactive`.  Then type `module load ccf_tools`.  Then type `python $CCF_TOOLS/run_all.py` [For the analysis of AAV experiments, use run_all_aav.py].  This will split the reads by barcode into different fasta files, batch them out to different nodes for mapping by bowtie2, and create a CCF file for each experiment in the barcodes.txt file.  The CCF files (with .ccf extension) are located in the output_and_analysis directory.   


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
