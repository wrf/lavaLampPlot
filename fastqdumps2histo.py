#! /usr/bin/env python
#
# fastqdumps2histo.py v1.1 last modified 2015-04-02
# by WRF

"""
fastqdumps2histo.py v1.1 2015-04-02

    first generate fastq.counts file from zipped reads with jellyfish count
gzip -dc reads.fastq.gz | jellyfish count -m 25 -o fastq.counts -C -U 1000 -s 1000000000 /dev/fd/0

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

    histo.csv is a comma separated matrix by GC count and coverage
    coverage is taken from the fasta header for the kmer, which appears as:
>123
ACTTGATCGTGATGCTAGTAGCTGT

    GC counts are integer counts of all Gs plus Cs
    so the max value is the kmer length

this file can be imported directly into R for generating lava lamp plots
"""

import sys
import argparse
import time

def get_freq(line):
	return int(line[1:].rstrip())

def get_gc(kmer):
	return kmer.count("G")+kmer.count("C")

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', type = argparse.FileType('rU'), default = '-', help="fasta format file")
	parser.add_argument('-k', '--kmer', type=int, metavar='N', default=25, help="kmer length [25]")
	parser.add_argument('-s', '--separator', default=",", help="use alternate separator [,]")
	parser.add_argument('-u', '--upper', type=int, metavar='N', default=1000, help="upper limit for histogram [1000]")
	args = parser.parse_args(argv)

	# generate kmer x maximum matrix
	m = [[0 for x in range(args.upper+1)] for y in range(args.kmer+1)]
	kmercount = 0

	print >> sys.stderr, "Reading fastq dumps", time.asctime()
	for line in args.input_file:
		# assumes jellyfish dump is in fasta format of kmers where header is count
		# for header, take the count
		if line[0]==">":
			kmercount += 1
			freq = get_freq(line)
		# if not a header, then it is sequence
		else:
			seq = line.rstrip()
			# calculate gc content for output
			gc = get_gc(seq)
			m[gc][freq] += 1

	print >> sys.stderr, "Writing counts", time.asctime()
	for i in m:
		print >> wayout, args.separator.join([str(j) for j in i])

	print >> sys.stderr, "Counted %d kmers" % (kmercount), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
