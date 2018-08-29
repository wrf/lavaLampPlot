#! /usr/bin/env python
#
'''
kmersorter.py v1.2 2018-08-29
# take jellyfish kmer dumps and sort by coverage

  1) using the compressed reads, count the kmers
gzip -dc reads1.fastq.gz reads2.fastq.gz | jellyfish count -m 31 -o fastq.counts -C -s 2000000000 -U 1000 -t 8 /dev/fd/0 

  2) generate a histogram (for your own reference)
jellyfish histo fastq.counts > fastq.counts.histo

  3) dump all kmers
jellyfish dump fastq.counts > fastq.counts.dumps

  ### if using Trinity mode, continue here
  # This is recommended to get better coverage slices if you have Trinity
  # these steps make use of various Trinity plugins and scripts
  4) run kmersorter with Trinity mode and all reads and the Trinity directory
kmersorter.py -T -k 31 -p 4 -a 100 -b 200 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.counts.dumps

  use -k to specify kmer length used by jellyfish
  use -p to specify max number of threads for some parallel processes
  # the amount of memory needed may be around the size of the fastq dump file
  # this could easily be 20-50Gb, so use a high memory machine

  # the above command will generate two files of filtered reads
  # the same command can be called, changingly only -a and -b to slice
    a different set of reads, without recreating any intermediate file

  # there is a known problem with SRA headers for the merging script
  # when extracting the reads, use an alternate command to change the header
fastq-dump --split-files --defline-seq "@$sn[_$rn]/$ri" SRR1032106.sra

  ### if using normal mode, continue after step 3
  4) extract the kmers with desired coverage
kmersorter.py -a 100 -b 200 -k 31 -l fastq.counts.dumps > filtered_kmers.fasta

  # in normal mode, these options can be set
  use -l to convert sequences to lowercase for downstream searching
  use -s and -w to filter by GC content
  # such as taking GC greater than 0.5 and less than 0.8
kmersorter.py -a 100 -b 200 -s 0.5 -w 0.8 -k 31 -l fastq.counts.dumps > filtered_kmers.fasta
'''

import sys
import argparse
import time
import os
import subprocess
from itertools import izip

def get_gc(kmer):
	return kmer.count("G")+kmer.count("C")

def check_gc(gc, strong, weak):
	return (weak >= gc >= strong)

def fastq_line_acc(line):
	return line.rstrip().split(" ")[0][1:]

def stats_to_dict(statfile, above, below):
	sd = {}
	readCount, addCount = 0,0
	print >> sys.stderr, "# Reading stats file", statfile, time.asctime()
	for line in open(statfile, 'r'):
		# in order to save memory, stats are filtered before adding to dictionary
		readCount += 1
		splits = line.split("\t")
		medCov = int(splits[0])
		if below >= medCov >= above:
			add_cov_to_acc(line, sd)
			addCount += 1
	print >> sys.stderr, "Read %d stats and %d pass coverage" % (readCount, addCount), time.asctime()
	return sd

def check_cov(acc, statDict):
	# since stats should be either the coverage value or not in the dictionary
	# freq can either take 0 or the coverage
	freq = statDict.get(acc, 0)
	# so the bool is either 1, passes coverage, or 0, not in dictionary
	return bool(freq)

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

def get_header_from_seqtype(seqtype):
	# returns two part header to correspond to FH and LC
	# one for the symbol of the fastx header
	# one for the number of lines per read
	if seqtype == "fastq":
		return "@", 4
	else:
		return ">", 2

def convert_to_fasta(fastool_path, reads):
	# this is one of two functions dependent on a Trinity program, calling fastool
	print >> sys.stderr, "### Converting reads to fasta", time.asctime()
	fasta_reads = "%s.fasta" % (os.path.splitext(reads)[0])
	if os.path.isfile(fasta_reads):
		print >> sys.stderr, "# Reads already converted, skipping...", time.asctime()
	else:
		fastool_args = [fastool_path, "-I", reads]
		with open(fasta_reads, 'a') as fr:
			subprocess.call(fastool_args, stdout=fr)
	return fasta_reads

def kmer_cov_to_stats(kmertocov_path, reads, jfdump, kmer, threads):
	# this function is dependent on the Trinity Inchworm program fastaToKmerCoverageStats
	stats_filename = "%s.k%d.stats" % (reads, kmer)
	if os.path.isfile(stats_filename):
		print >> sys.stderr, "# Coverage already counted for %s, skipping..." % (reads), time.asctime()
	else:
		print >> sys.stderr, "# Counting coverage on:", reads, time.asctime()
		kmertocov_args = [kmertocov_path, "--reads", reads, "--kmers", jfdump, "--kmer_size", "%d" % kmer, "--num_threads", "%d" % threads, "--DS"]
		with open(stats_filename, 'w') as sf:
			subprocess.call(kmertocov_args, stdout=sf)
	return stats_filename

def collect_reads(statdict1, statdict2, reads1, reads2, left_filt, right_filt, above, below, strong, weak, seqtype):
	# this part was formerly a blatant rip off of Trinity's nbkc_normalize.pl script
	# script calculated read pair average coverages, which caused problems with read extraction
	FH, LC = get_header_from_seqtype(seqtype)
	# linenum should have values only of 0 through LC (which is normally 4)
	linenum = 0
	pairCount, writePairs = 0, 0
	lineCount, writeCount = 0, 0
	keepPair, passGC = False, False
	print >> sys.stderr, "Stats contain %d and %d accessions" % (len(statdict1), len(statdict2) ), time.asctime()
	with open(left_filt, 'w') as lf, open(right_filt, 'w') as rf:
		for line1, line2 in izip( open(reads1,'r'), open(reads2,'r') ):
			lineCount += 1
			linenum += 1
			# check if this is the first line in the quartet, and that the first character is '@'
			if linenum == 1 and line1[0] == FH:
				pairCount += 1
				acc1, acc2 = line1, line2
			if linenum == 2:
				gc1, gc2 = get_gc(line1), get_gc(line2)
				# if either read is in both the selected cov and gc, keep the pair
				if (check_gc(gc1, strong, weak) and check_cov(fastq_line_acc(acc1), statdict1) ) or (check_gc(gc2, strong, weak) and check_cov(fastq_line_acc(acc2), statdict2) ):
					keepPair = True
					# if writing, add to writePairs which should finally equal the length of the dictionary
					writePairs += 1
					writeCount += 1
					# write the previous line and add 1 to the count if it passes
					lf.write(acc1)
					rf.write(acc2)
					# then write the current line containing the sequence, and carry on
			if keepPair:
				writeCount += 1
				lf.write(line1)
				rf.write(line2)
			# regardless if anything was written, reset linenum and passCoverage after 4 lines
			if linenum == LC:
				linenum = 0
				keepPair, passGC = False, False
	print >> sys.stderr, "Counted %d and wrote %d lines" % (lineCount, writeCount), time.asctime()
	print >> sys.stderr, "Counted %d sequence pairs" % (pairCount), time.asctime()
	print >> sys.stderr, "Kept %d read pairs" % (writePairs), time.asctime()

def collect_unpaired(statdict1, reads1, left_filt, above, below, strong, weak, seqtype):
	# this is to collect unpaired reads, perhaps sanger sequences or singletons
	# this assumes reads are long, and strong and weak should be percentages
	FH, LC = get_header_from_seqtype(seqtype)
	linenum = 0
	readCount, writeReads = 0, 0
	lineCount, writeCount = 0, 0
	keepSeq, passGC = False, False
	print >> sys.stderr, "Stats contain %d accessions" % (len(statdict1)), time.asctime()
	with open(left_filt, 'w') as lf:
		for line1 in open(reads1,'r'):
			lineCount += 1
			linenum += 1
			# check if this is the first line in the quartet, and that the first character is '@'
			if linenum == 1 and line1[0] == FH:
				readCount += 1
				acc1= line1
			if linenum == 2:
				gc1 = float(get_gc(line1))/len(line1.rstrip()) # should be a float between 0 and 1
				# if the read is in both the selected cov and gc, keep the pair
				if (check_gc(gc1, strong, weak) and check_cov(fastq_line_acc(acc1), statdict1) ):
					keepSeq = True
					# if writing, add to writeReads
					writeReads += 1
					writeCount += 1
					# write the previous line and add 1 to the count if it passes
					lf.write(acc1)
					# then write the current line containing the sequence, and carry on
			if keepSeq:
				writeCount += 1
				lf.write(line1)
			# regardless if anything was written, reset linenum and passCoverage after 4 lines
			if linenum == LC:
				linenum = 0
				keepSeq, passGC = False, False
	print >> sys.stderr, "Counted %d and wrote %d lines" % (lineCount, writeCount), time.asctime()
	print >> sys.stderr, "Counted %d sequence seqs" % (readCount), time.asctime()
	print >> sys.stderr, "Kept %d reads" % (writeReads), time.asctime()

def get_fastx_type(seqsFile):
	with open(seqsFile, 'r') as lr:
		header = lr.readline()[0]
		seqLen = len(lr.readline().rstrip())
	return header, seqLen

def fastx_header_to_type(header):
	if header == ">":
		return "fasta"
	elif header == "@":
		return "fastq"
	else:
		print >> sys.stderr, "# Error unknown header type, %s" % (header)
		print >> sys.stderr, "# Exiting", time.asctime()
		sys.exit()

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")

	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('input_file', help="jellyfish dump of fasta format kmers")
	parser.add_argument('-a', '--above', type=int, metavar='N', default=1, help="retrieve kmers with a count greater than N [1]")
	parser.add_argument('-b', '--below', type=int, metavar='N', default=10000, help="retrieve kmers with a count less than N [10000]")
	parser.add_argument('-s', '--strong', type=float, metavar='0.F', default=0.0, help="retrieve kmers with GC percent greater than N [0.0]")
	parser.add_argument('-w', '--weak', type=float, metavar='0.F', default=1.0, help="retrieve kmers with GC percent less than N [1.0]")
	parser.add_argument('-t', '--type', help="input files and line numbering are fasta or fastq [auto-detect]", default='auto')
	parser.add_argument('-k', '--kmer', type=int, metavar='N', default=25, help="kmer length [25]")
	parser.add_argument('-l', '--lowercase', action="store_true", help="convert sequences to lowercase")
	parser.add_argument('-T', '--trinity', action="store_true", help="use Trinity method to extract reads")
	parser.add_argument('-D', '--directory', help="trinity directory")
	parser.add_argument('-1', '--left-reads', help="fastq reads 1 of a pair, where length is constant")
	parser.add_argument('-2', '--right-reads', help="fastq reads 2 of a pair, where length is constant")
	parser.add_argument('-L', '--long-reads', help="file of long reads, contigs, or single end reads")
	parser.add_argument('-c', '--stdev-cutoff', type=int, metavar='N', default=200, help="upper limit for standard deviation cutoff [200]")
	parser.add_argument('-p', '--processors', type=int, metavar='N', default=1, help="number of CPUs [1]")
	parser.add_argument('-S', '--stats', action="store_false", help="only collect left and right stats in Trinity mode")
	args = parser.parse_args(argv)

	keepcov = 0
	keepgc = 0
	kmercount = 0
	freqsum = 0

	# checks for normal data
	with open(args.input_file,'r') as kf:
		headerType = kf.readline()[0]
		if not headerType == ">":
			print >> sys.stderr, "# ERROR kmers not in fasta format"
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()
		kmerlength = len(kf.readline().rstrip())
		if not kmerlength == args.kmer:
			print >> sys.stderr, "# ERROR kmer length %d incorrect for file %s" % (args.kmer, args.input_file)
			print >> sys.stderr, "# kmer length appears to be %d, set this with -k" % (kmerlength )
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

	# check that strong and weak are either floats between 0 and 1, or integers
	if args.strong > 1 or args.weak > 1:
		if args.long_reads:
			print >> sys.stderr, "# ERROR for long reads -L, -s and -w must be a float between 0.0 and 1.0"
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()
		if int(args.strong)!=args.strong or int(args.weak)!=args.weak:
			print >> sys.stderr, "# ERROR -s and -w must be a float between 0.0 and 1.0 or integers of raw GC counts"
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

	if args.trinity:
		# PART I - checking setup before trying anything
		print >> sys.stderr, "### Running Trinity type read sorting", time.asctime()
		if args.directory and os.path.isdir(args.directory):
			print >> sys.stderr, "# Using Trinity directory %s" % (args.directory), time.asctime()
			trinity_dir = os.path.abspath(args.directory)
		else:
			print >> sys.stderr, "# Must provide a valid directory -D", time.asctime()
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

		# check if all programs are there
		print >> sys.stderr, "### Checking for sub programs", time.asctime()
		fastool_path = os.path.join(trinity_dir, "util/support_scripts/fastQ_to_fastA.pl")
		kmertocov_path = os.path.join(trinity_dir, "Inchworm/bin/fastaToKmerCoverageStats")
		allProgsFound = True
		for prog in [fastool_path, kmertocov_path]:
			if not os.path.isfile(prog):
				print >> sys.stderr, "# Cannot find program %s" % (prog), time.asctime()
				allProgsFound = False
			else:
				print >> sys.stderr, "# %s found... OK" % (prog)
		if not allProgsFound:
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

		print >> sys.stderr, "### Checking input files", time.asctime()
		if args.left_reads and args.right_reads:
			if os.path.isfile(args.left_reads) and os.path.isfile(args.right_reads):
				print >> sys.stderr, "# Reads found... OK"
				headerType, seqLength = get_fastx_type(args.left_reads)
				if not args.long_reads:
					print >> sys.stderr, "# Detected read length of %d... OK" % (seqLength)
				if args.type=="auto":
					seqType = fastx_header_to_type(headerType)
					print >> sys.stderr, "# Detected seq type %s... OK" % (seqType)
				elif args.type=="fastq":
					seqType = args.type
				elif args.type=="fasta":
					seqType = args.type
				else:
					print >> sys.stderr, "# Error unknown sequence type, %s" % (args.type)
					print >> sys.stderr, "# Exiting", time.asctime()
					sys.exit()
		if args.long_reads and os.path.isfile(args.long_reads):
			print >> sys.stderr, "# Long reads found... OK"
			longheaderType, longseqLength = get_fastx_type(args.long_reads)

		## PART II - actually running the programs
		# most of this pipeline is copied from Trinity insilico_read_normalization.pl
		print >> sys.stderr, "### Converting reads to fasta", time.asctime()
		if args.left_reads and args.right_reads:
			if seqType == "fasta":
				print >> sys.stderr, "# Input type already fasta, skipping...", time.asctime()
				left_reads = args.left_reads
				right_reads = args.right_reads
			else:
				left_reads = convert_to_fasta(fastool_path, args.left_reads)
				right_reads = convert_to_fasta(fastool_path, args.right_reads)
		if args.long_reads:
			longseqType = fastx_header_to_type(longheaderType)
			if longseqType=="fastq":
				long_reads = convert_to_fasta(fastool_path, args.long_reads)
			else:
				print >> sys.stderr, "# Input type already fasta, skipping...", time.asctime()
				long_reads = args.long_reads

		print >> sys.stderr, "### Generating coverage stats", time.asctime()
		if args.left_reads and os.path.isfile(left_reads):
			left_stats = kmer_cov_to_stats(kmertocov_path, left_reads, args.input_file, args.kmer, args.processors)
		if args.right_reads and os.path.isfile(right_reads):
			right_stats = kmer_cov_to_stats(kmertocov_path, right_reads, args.input_file, args.kmer, args.processors)
		if args.long_reads and os.path.isfile(long_reads):
			long_stats = kmer_cov_to_stats(kmertocov_path, long_reads, args.input_file, args.kmer, args.processors)

		if args.stats:
			print >> sys.stderr, "### Sorting reads by coverage and/or GC", time.asctime()
			if args.above==1 and args.below==10000 and args.strong==0.0 and args.weak==1.0:
				print >> sys.stderr, "# No cutoffs specified, -a -b -s -w", time.asctime()
			else:
				# this is a list of an empty string so that it can be string-joined later
				coverage_mods = ['']
				if not args.above==1 or not args.below==10000:
					coverage_mods.append("a%d.b%d" % (args.above,args.below) )
				print >> sys.stderr, "# Filtering coverage between %d and %d" % (args.above, args.below), time.asctime()
				# set GC limits, int is needed in case read length is not 100
				if (int(args.strong)==args.strong and int(args.weak)==args.weak):
					strong, weak = args.strong, args.weak
					print >> sys.stderr, "# Filtering GC between %d and %d" % (strong, weak), time.asctime()
				else:
					if args.long_reads:
						lstrong, lweak = args.strong, args.weak
						print >> sys.stderr, "# Filtering GC between %.3f and %.3f" % (lstrong, lweak), time.asctime()
					else:
						strong = int(seqLength * args.strong + 0.9)
						weak = int(seqLength * args.weak + 0.9)
						print >> sys.stderr, "# Filtering GC between %d and %d" % (strong, weak), time.asctime()
				if not args.strong==0.0 or not args.weak==1.0:
					if args.long_reads:
						coverage_mods.append("s%d.w%d" % (lstrong*100,lweak*100) )
					else:
						coverage_mods.append("s%d.w%d" % (strong,weak) )

				# generate output file names
				coverage_string = ".".join(coverage_mods)

				print >> sys.stderr, "### Retrieving reads", time.asctime()
				if args.left_reads and args.right_reads:
					statDict1 = stats_to_dict(left_stats, args.above, args.below)
					statDict2 = stats_to_dict(right_stats, args.above, args.below)
					left_filt = "%s.k%d%s.%s" % (os.path.splitext(args.left_reads)[0], args.kmer, coverage_string, seqType)
					right_filt = "%s.k%d%s.%s" % (os.path.splitext(args.right_reads)[0], args.kmer, coverage_string, seqType)
					collect_reads(statDict1, statDict2, args.left_reads, args.right_reads, left_filt, right_filt, args.above, args.below, strong, weak, seqType)
				if args.long_reads and long_stats:
					longstatDict = stats_to_dict(long_stats, args.above, args.below)
					long_filt = "%s.k%d%s.%s" % (os.path.splitext(args.long_reads)[0], args.kmer, coverage_string, longseqType)
					collect_unpaired(longstatDict, args.long_reads, long_filt, args.above, args.below, lstrong, lweak, longseqType)
		print >> sys.stderr, "# Process finished", time.asctime()

	# in kmer mode:
	# lowercase is needed for using grep to quickly pull out sequences,
	# as otherwise those sequences can be found randomly in the quality scores
	# for example
	# $ grep ggggggggggggggggggggggg filtered_kmers.lowercase.fasta
	#gggggggggggggggggggggggggagccga
	#ggggggggggggggggggggggggggaggga   x16 uppercase as quality score
	#gggggggggggggggggggggggggggggca   x469 all Gs with 2 others
	#gggggggggggggggggggggggggtgggga
	#ggggggggggggggggggggggggagccgaa
	#gggggggggggggggggggggggaggtggca
	#gggggggggggggggggggggggggggggaa   x610
	#aaggggggggggggggggggggggggggggg   x291
	#
	else:
		print >> sys.stderr, "### Sorting kmers from counts", time.asctime()
		print >> sys.stderr, "# Filtering coverage between %d and %d" % (args.above, args.below), time.asctime()
		print >> sys.stderr, "# Filtering GC between %.2f and %.2f" % (args.strong, args.weak), time.asctime()
		if args.lowercase:
			print >> sys.stderr, "# Converting kmers sequences to lowercase", time.asctime()
		with open(args.input_file,'r') as kf:
			for line in kf:
				# assumes jellyfish dump is in fasta format of kmers where header is count
				# for header, take the count
				if line[0]==">":
					kmercount += 1
					freq = int(line[1:].rstrip())
				# if not a header, then it is sequence
				else:
					seq = line.rstrip()
					# filter by freq, given by the options
					if args.below > freq >= args.above:
						keepcov += 1
						# calculate gc content for output
						gc = (seq.count("G")+seq.count("C"))/float(args.kmer)
						# also filter by GC content as needed
						if args.weak >= gc >= args.strong:
							keepgc += 1
							freqsum += freq
							# write out fasta header as number_frequency_gc
							print >> wayout, ">%d_%d_%.3f" % (kmercount, freq, gc)
							# convert to lowercase for certain post filtering
							if args.lowercase:
								print >> wayout, seq.lower()
							else:
								print >> wayout, seq

		print >> sys.stderr, "# Counted %d kmers" % (kmercount), time.asctime()
		print >> sys.stderr, "# Found %d within coverage, and wrote %d" % (keepcov, keepgc), time.asctime()
		print >> sys.stderr, "# Average kept frequency (kmers * coverage) is %.2f" % (float(freqsum)/kmercount ), time.asctime()

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
