# PiggyBac_insertion_mapping
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
