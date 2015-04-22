#! /usr/bin/env python
#
# added Trinity mode for better coverage slices 2015-04-19
# kmersorter.py v1.1 2015-03-03
# sorts kmers based on coverage from jellyfish dump

'''
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
kmersorter.py -T -k 31 -m 20G -p 4 -a 100 -b 200 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.counts.dumps

  use -k to specify kmer length used by jellyfish
  use -m to specify the max memory usage for the `sort` command
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

def filter_coverage(pairstats, above, below, keptAccs, stdCutoff):
	# this part is a blatant rip off of Trinity's nbkc_normalize.pl script
	total_lines = 0
	nocovCount = 0
	aberrant = 0
	kept_count = 0
	accessionDict = {}
	with open(pairstats, 'r') as ps, open(keptAccs, 'w') as fa:
		for line in ps:
			total_lines += 1
			medCov, avgCov, stdev, pctDev, acc = line.rstrip().split('\t')
			medCov = int(medCov)
			pctDev = int(pctDev)
			if medCov < 1:
				nocovCount +=1
				continue
			if pctDev >= stdCutoff:
				aberrant += 1
				continue
			if below >= medCov >= above:
				kept_count += 1
				accessionDict[acc] = True
				print >> fa, acc
	print >> sys.stderr, "{} / {} = {:.2f} % reads kept after filtering".format(kept_count, total_lines, kept_count*100.0/total_lines )
	print >> sys.stderr, "{} / {} = {:.2f} % reads with aberrant coverage profiles".format(aberrant, total_lines, aberrant*100.0/total_lines )
	print >> sys.stderr, "{} / {} = {:.2f} % reads with no coverage".format(nocovCount, total_lines, nocovCount*100.0/total_lines )
	return accessionDict

def collect_reads(accessionDict, left, right, left_filt, right_filt, strong, weak, seqtype):
	if seqtype == "fastq":
		FH = "@"
		LC = 4
	else:
		FH = ">"
		LC = 2
	# linenum should have values only of 0 through LC (which is normally 4)
	linenum = 0
	pairCount, writePairs = 0, 0
	lineCount, writeCount = 0, 0
	covPassCount = 0
	passCoverage, passGC = False, False
	print >> sys.stderr, "Retained %d accessions" % (len(accessionDict) ), time.asctime()
	with open(left_filt, 'w') as lf, open(right_filt, 'w') as rf:
		for line1, line2 in izip( open(left,'r'), open(right,'r') ):
			lineCount += 1
			linenum += 1
			# check if this is the first line in the quartet, and that the first character is '@'
			if linenum == 1 and line1[0] == FH:
				pairCount += 1
				acc1, acc2 = line1, line2
				# if that acc is in the dictionary of reads to keep from the previous step
				# it should return True, otherwise the default is False
				passCoverage = accessionDict.get(fastq_line_acc(line1), False)
				# this will either add 0 or 1 to covPassCount
				covPassCount += int(passCoverage)
			if passCoverage and linenum == 2:
				gc1, gc2 = get_gc(line1), get_gc(line2)
				if check_gc(gc1, strong, weak) or check_gc(gc2, strong, weak):
					passGC = True
					# if writing, add to writePairs which should finally equal the length of the dictionary
					writePairs += 1
					writeCount += 1
					# write the previous line and add 1 to the count if it passes
					lf.write(acc1)
					rf.write(acc2)
					# then write the current line containing the sequence, and carry on
			if passCoverage and passGC:
				writeCount += 1
				lf.write(line1)
				rf.write(line2)
			# regardless if anything was written, reset linenum and passCoverage after 4 lines
			if linenum == LC:
				linenum = 0
				passCoverage, passGC = False, False
	print >> sys.stderr, "Counted %d and wrote %d lines" % (lineCount, writeCount), time.asctime()
	print >> sys.stderr, "Counted %d sequence pairs" % (pairCount), time.asctime()
	print >> sys.stderr, "%d sequence pairs passed coverage, %d passed GC" % (covPassCount, writeCount), time.asctime()

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
	parser.add_argument('-1', '--left-reads', help="fastq reads 1 of a pair")
	parser.add_argument('-2', '--right-reads', help="fastq reads 2 of a pair")
	parser.add_argument('-T', '--trinity', action="store_true", help="use Trinity method to extract reads")
	parser.add_argument('-D', '--directory', help="trinity directory")
	parser.add_argument('-c', '--stdev-cutoff', type=int, metavar='N', default=200, help="upper limit for standard deviation cutoff [200]")
	parser.add_argument('-m', '--memory', default="1G", help="memory usage maximum [1G]")
	parser.add_argument('-p', '--processors', type=int, metavar='N', default=1, help="number of CPUs [1]")
	parser.add_argument('-S', '--stats', action="store_false", help="only collect left and right stats in Trinity mode")
	args = parser.parse_args(argv)

	keepcov = 0
	keepgc = 0
	kmercount = 0
	freqsum = 0

	if args.strong > 1 or args.weak > 1:
		print >> sys.stderr, "# ERROR -s and -w must be a float between 0.0 and 1.0"
		print >> sys.stderr, "# Exiting", time.asctime()
		sys.exit()

	if args.trinity:
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
		fastool_path = os.path.join(trinity_dir, "trinity-plugins/fastool/fastool")
		kmertocov_path = os.path.join(trinity_dir, "Inchworm/bin/fastaToKmerCoverageStats")
		mergeLR_path = os.path.join(trinity_dir, "util/support_scripts/nbkc_merge_left_right_stats.pl")
		allProgsFound = True
		for prog in [fastool_path, kmertocov_path, mergeLR_path]:
			if not os.path.isfile(prog):
				print >> sys.stderr, "# Cannot find program %s" % (prog), time.asctime()
				allProgsFound = False
			else:
				print >> sys.stderr, "# %s found... OK" % (prog)
		if not allProgsFound:
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

		print >> sys.stderr, "### Checking input files", time.asctime()
		if os.path.isfile(args.left_reads) and os.path.isfile(args.right_reads):
			print >> sys.stderr, "# Reads found... OK"
			headerType, seqLength = get_fastx_type(args.left_reads)
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
		else:
			print >> sys.stderr, "# Cannot find input reads %s, %s" % (args.left_reads, args.right_reads)
			print >> sys.stderr, "# Exiting", time.asctime()
			sys.exit()

		# most of this pipeline is copied from Trinity insilico_read_normalization.pl
		print >> sys.stderr, "### Converting reads to fasta", time.asctime()
		if seqType == "fasta":
			print >> sys.stderr, "# Input type already fasta, skipping...", time.asctime()
			left_reads = args.left_reads
			right_reads = args.right_reads
		else:
			left_reads = "%s.left.fa" % (os.path.splitext(args.left_reads)[0])
			right_reads = "%s.right.fa" % (os.path.splitext(args.right_reads)[0])
			if os.path.isfile(left_reads) and os.path.isfile(right_reads):
				print >> sys.stderr, "# Reads already converted, skipping...", time.asctime()
			else:
				fastool_args = [fastool_path, "--illumina-trinity", "--to-fasta", args.left_reads]
				with open(left_reads, 'a') as lr:
					subprocess.call(fastool_args, stdout=lr)
				fastool_args = [fastool_path, "--illumina-trinity", "--to-fasta", args.right_reads]
				with open(right_reads, 'a') as rr:
					subprocess.call(fastool_args, stdout=rr)

		print >> sys.stderr, "### Generating coverage stats", time.asctime()
		left_stats = "%s.k%d.stats" % (left_reads, args.kmer)
		right_stats = "%s.k%d.stats" % (right_reads, args.kmer)
		if os.path.isfile(left_stats) and os.path.isfile(right_stats):
			print >> sys.stderr, "# Coverage already counted, skipping...", time.asctime()
		else:
			print >> sys.stderr, "# Counting left coverage:", left_reads, time.asctime()
			kmertocov_args = [kmertocov_path, "--reads", left_reads, "--kmers", args.input_file, "--kmer_size", "%d" % args.kmer, "--num_threads", "%d" % args.processors, "--DS"]
			with open(left_stats, 'w') as ls:
				subprocess.call(kmertocov_args, stdout=ls)
			print >> sys.stderr, "# Counting right coverage:", right_reads, time.asctime()
			kmertocov_args = [kmertocov_path, "--reads", right_reads, "--kmers", args.input_file, "--kmer_size", "%d" %  args.kmer, "--num_threads", "%d" % args.processors, "--DS"]
			with open(right_stats, 'w') as rs:
				subprocess.call(kmertocov_args, stdout=rs)

		if args.stats:
			# the value of this step is in question
			print >> sys.stderr, "### Sorting stats by read name", time.asctime()
			left_sorted = "%s.sort" % left_stats
			right_sorted = "%s.sort" % right_stats
			if os.path.isfile(left_sorted) and os.path.isfile(right_sorted):
				print >> sys.stderr, "# Stats already sorted, skipping...", time.asctime()
			else:
				print >> sys.stderr, "# Sorting left stats:", left_stats, time.asctime()
				sort_args = ["/usr/bin/sort","--parallel=%d" % args.processors, "-k5,5", "-T", ".", "-S", args.memory, left_stats]
				with open(left_sorted, 'w') as lo:
					subprocess.call(sort_args, stdout=lo)
				print >> sys.stderr, "# Sorting right stats:", right_stats, time.asctime()
				sort_args = ["/usr/bin/sort","--parallel=%d" % args.processors, "-k5,5", "-T", ".", "-S", args.memory, right_stats]
				with open(right_sorted, 'w') as ro:
					subprocess.call(sort_args, stdout=ro)

			print >> sys.stderr, "### Merging left and right stats", time.asctime()
			# this merging takes the average of the left and right, rounding up
			paired_stats = "%s_pairs.k%d.stats" % (args.left_reads.split("_")[0], args.kmer)
			if os.path.isfile(paired_stats):
				print >> sys.stderr, "# Stats already merged, skipping...", time.asctime()
			else:
				mergeLR_args = ["perl", mergeLR_path, "--left", left_sorted, "--right", right_sorted, "--sorted"]
				with open(paired_stats, 'w') as ps:
					subprocess.call(mergeLR_args, stdout=ps)
				print >> sys.stderr, "# Done merging", time.asctime()

			print >> sys.stderr, "### Sorting read pairs by coverage", time.asctime()
			print >> sys.stderr, "# Filtering coverage between %d and %d" % (args.above, args.below), time.asctime()
			core_accs = "%s.k%d.a%d.b%d.sd%d.accs" % (paired_stats, args.kmer, args.above, args.below, args.stdev_cutoff)
			accessions = filter_coverage(paired_stats, args.above, args.below, core_accs, args.stdev_cutoff)
			print >> sys.stderr, "# Done filtering", time.asctime()
			print >> sys.stderr, "### Retrieving reads", time.asctime()
			left_filt = "%s.k%d.a%d.b%d.sd%d.%s" % (os.path.splitext(args.left_reads)[0], args.kmer, args.above, args.below, args.stdev_cutoff, seqType)
			right_filt = "%s.k%d.a%d.b%d.sd%d.%s" % (os.path.splitext(args.right_reads)[0], args.kmer, args.above, args.below, args.stdev_cutoff, seqType)
			# set GC limits, int is needed in case read length is not 100
			strong = int(seqLength * args.strong + 0.9)
			weak = int(seqLength * args.weak + 0.9)
			print >> sys.stderr, "# Filtering GC between %d and %d" % (strong, weak), time.asctime()
			collect_reads(accessions, args.left_reads, args.right_reads, left_filt, right_filt, strong, weak, seqType)
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
