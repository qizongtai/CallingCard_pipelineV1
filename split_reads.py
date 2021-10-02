#!/usr/bin/env python
"""
split_reads.py
written 7/25/17 by RDM

usage 
split_reads.py -r1 <read1 file> -r2 <read2 file> 
-i <index file> -b<barcode file> -o <output path>
--hammp <hamming distance for primer barcode>
--hammt <hamming distance for transposon barcode>
--p3p <3' sequence of primer>

required fields:
    -r1 read 1 filename (full path)
    -i index filename (full path)

not required
    -r2 read 2 filename (full path)
    -b barcode file = ../raw/barcodes.txt
    -o output path = ../output_and_analysis
    -hp 0
    -tp 0
    --p3p <3' sequence of primer>


This program requires that the ccf_tools module is loaded.
"""
import argparse
import ccf_tools
import os


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='cc_filter_reads.py')
    parser.add_argument('-r1','--read1',help='Read 1 filename (full path)',required=True)
    parser.add_argument('-i','--indexfn',help='index filename (full path)',required=True)
    parser.add_argument('-o','--outputpath',help='output path',required=False,default='../output_and_analysis/')
    parser.add_argument('-b','--barcodefile',help='barcode filename (full path)',required=False,default='../scripts/barcodes.txt')
    parser.add_argument('-r2','--read2',help='Read2 filename (full path)',required=False,default=False)
    parser.add_argument('--hammp',help='Primer barcode hamming distance',required=False,default=0)
    parser.add_argument('--hammt',help='Transposon barcode hamming distance',required=False,default=0)
    parser.add_argument('--p3p',help='3 prime end of primer',required=False,default="GCGTCAATTTTACGCAGACTATCTTTCTAGGG")
    args = parser.parse_args()
    if not args.outputpath[-1] == "/":
        args.outputpath = args.outputpath+"/"
    os.chdir(args.outputpath)

    ccf_tools.filter_reads_V2_transposon_index(args.read1,args.indexfn,args.read2,args.barcodefile,args.outputpath,int(args.hammp),int(args.hammt),args.p3p)



