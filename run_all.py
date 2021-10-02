#!/bin/python
import commands, os
import argparse
import csv


def count_barcodes(barcode_filename):
#This helper function reads in the experiment name, the primer barcode, and the transposon barcodes
#from a file which is in the following format:
#expt name \t primer barcode \t transposon barcode 1,transposon barcode 2, transposon barcode 3 etc.
#It then returns a dictionary with key = a tuple (primer barcode, transposon barcode), value = expt_name
#The last line of the file should not have a return character
	reader = csv.reader(open(barcode_filename, 'r'),delimiter = '\t')
	d = {}
	count = 0
	for row in reader:
		exname,b1,b2 = row
		count = count+len(b2.split(","))
	return count

if __name__ == '__main__':
	parser = argparse.ArgumentParser(description='map_reads.py')
	parser.add_argument('-b','--barcodefilename',help='base filename (full path)',
		required=False,default='./barcodes.txt')
	args = parser.parse_args()
	total_barcodes = count_barcodes(args.barcodefilename)
	#submit the first job
	cmd = "sbatch split_reads.sh"
	print "Submitting Job1 with command: %s" % cmd
	status, jobnum = commands.getstatusoutput(cmd)
	jobnum = jobnum.split(" ")[-1]
	if (status == 0 ):
		print "Job1 is %s" % jobnum
	else:
		print "Error submitting Job1"
	# submit the second job to be dependent on the first
	cmd = "sbatch --depend=afterany:%s --array=1-%s map_reads.sh" % (jobnum,str(total_barcodes))
	print "Submitting Job2 with command: %s" % cmd
	status,jobnum = commands.getstatusoutput(cmd)
	jobnum = jobnum.split(" ")[-1]
	if (status == 0 ):
		print "Job2 is %s" % jobnum
	else:
		print "Error submitting Job2"

	# submit the third job (a swarm) to be dependent on the second
	cmd = "sbatch --depend=afterany:%s make_ccffile.sh" % jobnum
	print "Submitting swarm job  with command: %s" % cmd
	status,jobnum = commands.getstatusoutput(cmd)
	jobnum = jobnum.split(" ")[-1]
	if (status == 0 ):
		print "Job3 is %s" % jobnum
	else:
		print "Error submitting Job3"

	print "\nCurrent status:\n"
	#show the current status with 'sjobs'
	os.system("squeue -u rmitra")

