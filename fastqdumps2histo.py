#! /usr/bin/env python
#
# fastqdumps2histo.py v1.3
# by WRF

"""
fastqdumps2histo.py v1.3 last modified 2015-05-29

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

    because some downstream steps require the jellyfish dump, but can exclude
    entries with a value of 1, a filtered jellyfish dump can be created
    as it is read in the piped input, using the -j option

jellyfish dump fastq.counts | fastqdumps2histo.py -j fastq.dumps - > histo.csv

    otherwise use directly as an input file
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
fastqdumps2histo.py -f *.fastq -s *.stats -u 1000 -T - > histo.csv

    -r can set the read length, otherwise it is determined automatically
    -u is set as above
    -f are the original fastq reads (use wildcard or list multiple)
       this could also be a subset of the original reads
    -s are the .stats files for each read
    -t to specify fasta read files (-t fasta)
    -T - must set both -T for Trinity mode and - for null input
    -k is unused for this step

    for Trinity mode with long reads (such as sanger ESTs)
    use the -p option to calculate GC as percentage rather than raw count
fastqdumps2histo.py -f reads.fasta -s reads.stats -u 1000 -T -p - > histo.csv
"""

import sys
import argparse
import time
from itertools import izip

def get_freq(line):
	return int(line[1:].rstrip())

def get_gc(kmer):
	return kmer.count("G")+kmer.count("C")

def get_gc_perc_int(kmer):
	# N's must be substracted from total length, otherwise GC values will be too low for scaffolds
	return (kmer.count("G")+kmer.count("C")) * 100 / ( len(kmer.rstrip()) - kmer.count("N") )

def stats_to_dict(statfile):
	sd = {}
	for line in open(statfile, 'r'):
		add_cov_to_acc(line, sd)
	return sd

def add_cov_to_acc(line, statDict):
	lsplits = line.split("\t")
	# lines of stats file look like:
	#12	17.6632	10.3014	58.3215	HISEQ:150:C5KTJANXX:4:1101:1469:1959/1	thread:5
	# or this for SRA files:
	#53	57.4143	16.8252	29.3049	SRR1032106.2000.1/H	thread:5
	medCov = int(lsplits[0])
	acc = lsplits[4].split(" ")[0].rsplit("/",1)[0]
	# was previously using split("/")[0], caused errors with different split in fastx_line_acc
	statDict[acc] = medCov
	# nothing to return

def fastx_line_acc(line):
	return line.rstrip().split(" ")[0][1:]

def get_fastx_type(seqsFile, verbose):
	with open(seqsFile, 'r') as lr:
		headerLine = lr.readline()
		headerType = headerLine[0]
		seqLength = len(lr.readline().rstrip())
		if verbose:
			print >> sys.stderr, "Sequence header in format of: {}".format(headerLine.rstrip())
			print >> sys.stderr, "Header starts with: {}".format(headerType)
	return headerType, seqLength

def count_matrix(fastxfile, statDict, FH, LC, m, kmercount, maxcount, gcfunction):
	linecount = 0
	for line in open(fastxfile, 'r'):
		linecount += 1
		if linecount==1 and line[0]==FH:
			### TODO linecount should be reset when a new header line is found for fasta
			kmercount += 1
			acc = fastx_line_acc(line)
			freq = statDict.get(acc, 0)
		elif linecount==2:
			gc = gcfunction(line)
		# if it is the last line of each sequence, so cannot be elif
		if linecount==LC:
			# use try except in case the wrong -u is used or files were merged
			try:
				m[gc][freq] += 1
				#if freq == 1: # for debugging only
				#	print >> sys.stderr, acc
			except IndexError:
				# most likely reason for this is setting -k to kmer rather than read length
				maxcount += 1
			except UnboundLocalError:
				# this happened when freq was called before assignment
				print >> sys.stderr, "ERROR Check data type", time.asctime()
				# possibly wrong data type, no point continuing
				sys.exit()
			linecount = 0
	return m, kmercount, maxcount

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default='-', help="fasta format file")
	parser.add_argument('-d', '--delimiter', default=",", help="use alternate delimiter [,]")
	parser.add_argument('-f', '--fastx', nargs='*', help="fastx files, either fasta or fastq")
	parser.add_argument('-j', '--jellyfish', help="optional file for the filtered jellyfish dump")
	parser.add_argument('-k', '--kmer', type=int, metavar='N', default=25, help="kmer length [25]")
	parser.add_argument('-l', '--lower', type=int, metavar='N', default=2, help="lower limit for optional jellyfish dump [2]")
	parser.add_argument('-p', '--percentage', action="store_true", help="count GC as percentage for long reads, not by fixed read length (ignores -r)")
	parser.add_argument('-r', '--read-length', type=int, metavar='N', help="read length for Trinity mode [auto-detect]")
	parser.add_argument('-s', '--stats', nargs='*', help="Trinity stats files")
	parser.add_argument('-t', '--type', help="sequence type, fasta, fastq or [auto-detect]")
	parser.add_argument('-T', '--trinity', action="store_true", help="use Trinity method to extract reads")
	parser.add_argument('-u', '--upper', type=int, metavar='N', default=1000, help="upper limit for histogram [1000]")
	parser.add_argument('-v', '--verbose', action="store_true", help="extra output for debugging")
	args = parser.parse_args(argv)

	# gc and freq are set by default to 0, in case the reads cannot be parsed
	gc, freq = 0, 0
	kmercount = 0
	maxcount = 0
	kmertag = "kmers"

	if args.trinity:
		# checking for matching stats and fastx files
		if len(args.stats) != len(args.fastx):
			print >> sys.stderr, "ERROR Unequal numbers of stats and %s files" % args.type, time.asctime()
			print >> sys.stderr, "Exiting", time.asctime()
			sys.exit()
		# autodetect length and header type from fastx files
		FH, seqLen = get_fastx_type(args.fastx[0], args.verbose)
		# checking for sequence type and getting basic parameters of sequence parsing
		if (FH=="@" or args.type=="fastq"):
			LC = 4
		else: # otherwise args.type=="fasta"
			LC = 2
		# overwrite if read length is given
		if args.read_length:
			seqLen = args.read_length

		# set up different matrix if using percentage counting
		if args.percentage:
			# generate a matrix of 0 to 100 x cov maximum
			m = [[0 for x in range(args.upper+1)] for y in range(101)]
		else:
			# generate read length x cov maximum matrix, using read length as ymax
			m = [[0 for x in range(args.upper+1)] for y in range(seqLen+1)]
		if args.verbose:
			print >> sys.stderr, "Generating matrix of %d by %d" % (len(m), len(m[0]) )

		# iterate over pairs of fastx and stats files
		for statfile, fastxfile in izip(args.stats, args.fastx):
			print >> sys.stderr, "Reading stats file %s" % statfile, time.asctime()
			statDict = stats_to_dict(statfile)
			statCount = len(statDict)
			print >> sys.stderr, "Found stats for %d sequences" % statCount, time.asctime()
			if args.verbose:
				print >> sys.stderr, "Stats stored in format of: {}".format(statDict.iterkeys().next())
			print >> sys.stderr, "Reading file %s" % (fastxfile), time.asctime()
			if args.percentage:
				m, kmercount, maxcount = count_matrix(fastxfile, statDict, FH, LC, m, kmercount, maxcount, get_gc_perc_int)
			else:
				m, kmercount, maxcount = count_matrix(fastxfile, statDict, FH, LC, m, kmercount, maxcount, get_gc)
			if sum(x[0] for x in m)==kmercount:
				print >> sys.stderr, "ERROR No stats found for reads", time.asctime()
		kmertag = "reads"
	else:
		if args.jellyfish:
			jf = open(args.jellyfish, 'w')
			jfwroteback = 0
		# generate kmer x maximum matrix, using kmer length as ymax
		m = [[0 for x in range(args.upper+1)] for y in range(args.kmer+1)]
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
				if args.jellyfish and (args.upper >= freq >= args.lower):
					jf.write(">{}\n{}".format(freq, seq) )
					jfwroteback += 1
		if args.jellyfish:
			jf.close()
	print >> sys.stderr, "Writing counts", time.asctime()
	for i in m:
		print >> wayout, args.delimiter.join([str(j) for j in i])

	print >> sys.stderr, "Counted %d %s" % (kmercount, kmertag), time.asctime()
	if maxcount:
		print >> sys.stderr, "%d %s exceeded cutoffs" % (maxcount, kmertag), time.asctime()
	if not args.trinity and args.jellyfish:
		if jfwroteback:
			# only used in kmer counting mode
			print >> sys.stderr, "Wrote %d kmers to %s" % (jfwroteback, args.jellyfish), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
