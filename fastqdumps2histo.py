#! /usr/bin/env python
#
# fastqdumps2histo.py v1.2 last modified 2015-04-22
# by WRF

"""
fastqdumps2histo.py v1.2 2015-04-22

    first generate fastq.counts file from zipped reads with jellyfish count
gzip -dc reads.fastq.gz | jellyfish count -m 25 -o fastq.counts -C -U 1000 -s 1G /dev/fd/0

    add option -t for multithreads
    change -m for different kmer sizes
    change -U for maximum count, 1000 is good if genome coverage is 100+
    otherwise perhaps use a value around 500
    using a larger kmer will lower coverage, so if k > 45, set max to 500
    change -s for memory allocation (see jellyfish manual)

    fastq.counts file can be several gigabytes

    jellyfish dump output would normally be quite large, so it is piped
    then use jellyfish dump command as stdin
jellyfish dump fastq.counts | fastqdumps2histo.py - > histo.csv

    otherwise use as input file
fastqdumps2histo.py fastq.dumps > histo.csv

    if using a different kmer (other than 25, the default)
    change with -k to correspond to the -m of jellyfish
    also change -u to correspond to the -U of jellyfish
    change -u accordingly if multiple jellyfish counts are merged

    histo.csv is a comma separated matrix by GC count and coverage
    coverage is taken from the fasta header for the kmer, which appears as:
>123
ACTTGATCGTGATGCTAGTAGCTGT

    GC counts are integer counts of all Gs plus Cs
    so the max value is the kmer length

this file can be imported directly into R for generating lava lamp plots

    for alternate usage in Trinity mode to count read coverage and GC
    use the intermediate output files (.stats) from kmersorter.py
fastqdumps2histo.py -f *.fastq -s *.stats -k 100 -u 1000 -T - > histo.csv

    -k must be set to the read length, not the original kmer length
    -u is set as above
    -f are the original fastq reads (use wildcard or list multiple)
       this could also be a subset of the original reads
    -s are the .stats files for each read
    -t to use fasta read files (-t fasta)
    -T - must set both -T for Trinity mode and - for null input
"""

import sys
import argparse
import time
from itertools import izip

def get_freq(line):
	return int(line[1:].rstrip())

def get_gc(kmer):
	return kmer.count("G")+kmer.count("C")

def add_cov_to_acc(line, statDict):
	splits = line.split("\t")
	# lines of stats file look like:
	#12	17.6632	10.3014	58.3215	HISEQ:150:C5KTJANXX:4:1101:1469:1959/1	thread:5
	# or this for SRA files:
	#53	57.4143	16.8252	29.3049	SRR1032106.2000.1/H	thread:5
	medCov = int(splits[0])
	acc = splits[4].split("/")[0]
	statDict[acc] = medCov
	# nothing to return

def fastx_line_acc(line):
	return line.rstrip().split(" ")[0][1:]

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default='-', help="fasta format file")
	parser.add_argument('-d', '--delimiter', default=",", help="use alternate delimiter [,]")
	parser.add_argument('-f', '--fastx', nargs='*', help="fastx files, either fasta or fastq")
	parser.add_argument('-k', '--kmer', type=int, metavar='N', default=25, help="kmer length, or read length in Trinity mode [25]")
	parser.add_argument('-s', '--stats', nargs='*', help="Trinity stats files")
	parser.add_argument('-t', '--type', default="fastq", choices=["fasta", "fastq"], help="sequence type, fasta or [fastq]")
	parser.add_argument('-T', '--trinity', action="store_true", help="use Trinity method to extract reads")
	parser.add_argument('-u', '--upper', type=int, metavar='N', default=1000, help="upper limit for histogram [1000]")
	args = parser.parse_args(argv)

	# generate kmer x maximum matrix
	m = [[0 for x in range(args.upper+1)] for y in range(args.kmer+1)]
	kmercount = 0
	maxcount = 0
	kmertag = "kmers"

	if args.trinity:
		# checking for sequence type and getting basic parameters of sequence parsing
		if args.type=="fastq":
			LC = 4
			FH = "@"
		else: # otherwise args.type=="fasta"
			LC = 2
			FH = ">"
		# checking for matching stats and fastx files
		if len(args.stats) != len(args.fastx):
			print >> sys.stderr, "ERROR Unequal numbers of stats and %s files" % args.type, time.asctime()
			print >> sys.stderr, "Exiting", time.asctime()
			sys.exit()
		for statfile, fastxfile in izip(args.stats, args.fastx):
			statDict = {}
			print >> sys.stderr, "Reading stats file %s" % statfile, time.asctime()
			for line in open(statfile, 'r'):
				add_cov_to_acc(line, statDict)
			print >> sys.stderr, "Reading %s file %s" % (args.type, fastxfile), time.asctime()
			linecount = 0
			for line in open(fastxfile, 'r'):
				linecount += 1
				if linecount==1 and line[0]==FH:
					kmercount += 1
					acc = fastx_line_acc(line)
					freq = statDict.get(acc, 0)
				if linecount==2:
					gc = get_gc(line)
				if linecount==LC:
					# use try except in case the wrong -u is used or files were merged
					try:
						m[gc][freq] += 1
					except IndexError:
						# most likely reason for this is setting -k to kmer rather than read length
						maxcount += 1
					linecount = 0
		kmertag = "reads"
	else:
		print >> sys.stderr, "Reading fastq dumps", time.asctime()
		for line in args.input_file:
			# assumes jellyfish dump is in fasta format of kmers where header is count
			# for header, take the count
			if line[0]==">":
				kmercount += 1
				freq = get_freq(line)
			# if not a header, then it is sequence
			else:
				seq = line
				# calculate gc content for output
				gc = get_gc(seq)
				# use try except in case the wrong -u is used or files were merged
				try:
					m[gc][freq] += 1
				except IndexError:
					# any values outside of the range of the matrix are thus ignored
					maxcount += 1

	print >> sys.stderr, "Writing counts", time.asctime()
	for i in m:
		print >> wayout, args.delimiter.join([str(j) for j in i])

	print >> sys.stderr, "Counted %d %s" % (kmercount, kmertag), time.asctime()
	if maxcount:
		print >> sys.stderr, "%d %s exceeded cutoffs" % (maxcount, kmertag), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
