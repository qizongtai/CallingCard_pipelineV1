#!/usr/bin/env python
"""
make_gnashyfile.py
written 7/25/17 by RDM

usage
make_ccffile -b <barcode filename> -g <genome>  mm10 or hg19 or hg38 default = mm10 -o <output path>

required
none

not required
-b <barcode filename>, default = ../output_and_analysis/barcodes.txt
-p <path>, default = ../output_and_analysis/  path for input and output

"""

import argparse
import ccf_tools
import os

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='make_gnashyfile.py')
    parser.add_argument('-b','--barcodefile',help='barcode filename (full path)',required=False,default='../scripts/barcodes.txt')
    parser.add_argument('-p','--outputpath',help='output path',required=False,default='../output_and_analysis')
    parser.add_argument('-g','--genome',help='genome to be analyzed',required=True)
    args = parser.parse_args()
    if not args.outputpath[-1] == "/":
        args.outputpath = args.outputpath+"/"
    ccf_tools.make_ccffile(args.genome,args.barcodefile,args.outputpath)


                            





