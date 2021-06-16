# CallingCards piggyBac transposon genomic mapping
# M Yen, Z Qi et al. (2018) Genomic targets of BAP1 in uveal melanoma identified by transposase mapping. BMC Medical Genomics 11 (1), 97
# Z Qi et al. An optimized, broadly applicable piggyBac transposon induction system. (2017) Nucleic Acids Research 45 (7), e55-e55

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
