"""
map_reads.py
written 6/30/16 by RDM
usage
updated 7/7/16 by RDM to account for new PB primers with longer
priming sequence (default trim is now 37)
map_reads -f <base file name> -g <genome> -p <paired end flag> 
-ft <five primer bases to trim> -tt <three primer bases to trim>
-q <quality filter> -o <output path>

required
-f <base file name> -g <genome> 

not required
-p <paired end flag> default = False 
-ft <five primer bases to trim> default =37  
-tt <three primer bases to trim> default = 0
-q <quality filter>  default = 10
-o <output path> default = ../output_and_analysis



This program requires that the bowtie2 module is loaded and
that the following environment variable is set:
    export BOWTIE2_INDEXES='/scratch/rmlab/ref/bowtie2_indexes/mm10'

It also requires the samtools module is loaded.
    """

import argparse
import os

def map_reads(basefilename,genome,paired,trim_five,trim_three,qcutoff,outpath):
    outfilename = basefilename+".bam"
    outerrname = basefilename+".err"
    r1_filename = basefilename+"_R1.fastq"
    
    if paired:
        r2_filename = basefilename+"_R2.fastq"
        bowtie2_string = "bowtie2 -x "+genome+" -1 "+r1_filename+" -2 "+r2_filename+" --trim5 "+str(trim_five)+" --trim3 "+str(trim_three)+" 2> "+outerrname+" |samtools view -bS -q"+str(qcutoff)+" > "+outfilename    
    else:
        bowtie2_string = "bowtie2 -x "+genome+" -U "+r1_filename+" --trim5 "+str(trim_five)+" --trim3 "+str(trim_three)+" 2> "+outerrname+" |samtools view -bS -q"+str(qcutoff)+" > "+outfilename
    print bowtie2_string
    os.system(bowtie2_string)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='map_reads.py')
    parser.add_argument('-f','--basename',help='base filename (full path)',required=True)
    parser.add_argument('-g','--genome',help='genome name',required=True)
    parser.add_argument('-p','--paired',help='paired read flag',required=False,default=False)
    parser.add_argument('-ft','--ftrim',help='five prime bases to trim',default=37)
    parser.add_argument('-tt','--ttrim',help='three prime bases to trim',required=False,default=0)
    parser.add_argument('-q','--quality',help='quality score cutoff',required=False,default=10)
    parser.add_argument('-o','--outpath',help='output path',required=False,default='../output_and_analysis')
    args = parser.parse_args()
    if args.paired == "False":
        args.paired=False
    if not args.outpath[-1] == "/":
        args.outpath = args.outpath+"/"
    os.chdir(args.outpath)
    map_reads(args.basename,args.genome,args.paired,args.ftrim,args.ttrim,args.quality,args.outpath)

