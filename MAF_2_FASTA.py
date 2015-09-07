#!/usr/bin/python

import sys
import os
import getopt
sys.path.insert(1, "alignio-maf")
from Bio.AlignIO import MafIO
from Bio import AlignIO

def usage():
	print "=============================================================================="
	print " MAF (Multiple Alignment File) to FASTA "
	print "==============================================================================\n"
	print "Ben J. Ward, 2015\n"
	print "Convert aligned genome segments in a Multiple Alignment file and generate FASTA format files."
	print "For every aligned segment.\n"
	print "Usage:\n"
	print "-h --help    | Display this message."
	print "-f --file    | Specify the Multiple Alignment File to be processed."
	print "-l --log     | Optional: Specify filename for the log."
	print "-o --outfile | Optional: specify an output file/directory."

def maf_to_fasta(inputPath, logPath, outPath):
	i = 0
	log = open(logPath, 'w')
	log.write("file, length\n")
	for multiple_alignment in AlignIO.parse(inputPath, "maf"):
		AlignIO.write(multiple_alignment, "%s/%s.fa" % (outPath, i), "fasta")
		log.write("%s.fa, %s\n" % (i, multiple_alignment[1].annotations["size"]))
		i = i+1
	log.close()

def main(argv):
	try:
		opts, args = getopt.getopt(argv, "hl:f:o:", ["help", "log", "file", "outfile"])
	except getopt.GetoptError:
		usage()
		sys.exit(2)
	autoLogFileName = True
	logFilePath = ""
	inputPath = ""
	outputDir = ""
	for opt, arg in opts:
		if opt in ("-h", "--help"):
			usage()
			sys.exit()
		elif opt in ("-l", "--log"):
			logFilePath = os.path.abspath(os.path.expanduser(arg))
			autoLogFileName = False
		elif opt in ("-f", "--file"):
			inputPath = os.path.abspath(os.path.expanduser(arg))
		elif opt in ("-o", "--outfile"):
			outputDir = os.path.abspath(os.path.expanduser(arg))
	if inputPath == "":
		print "No input file has been provided. Exiting."
		sys.exit(2)
	if outputDir == "":
		outputDir = "%s_fastaFiles" % (os.path.splitext(inputPath)[0])
	if not os.path.exists(outputDir):
		os.makedirs(outputDir)
	if autoLogFileName:
		logFilePath = "%s/log.csv" % (outputDir)
	print inputPath
	print outputDir
	print logFilePath
	maf_to_fasta(inputPath, logFilePath, outputDir)

if __name__ == "__main__":
	main(sys.argv[1:])
	sys.exit()
