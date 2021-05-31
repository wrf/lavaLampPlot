#!/usr/bin/env python

'''
sort_reads_from_bam.py v0.1 2021-05-31

    can take stdin for BAM files, and allow reads as .gz:
samtools view hits.bam | sort_reads_from_bam.py -i - -c bacterial_contigs -1 reads1.fq.gz -2 reads2.fq.gz

    output files are automatically named by -o
    so would be reads1.fq.kept.gz and reads2.fq.kept.gz

    -c contigs file can be delimited names, or fasta headers by grep ">"
    for tab as a delimiter, used -d $'\\t' with SINGLE QUOTES
    -c can also be the name of a single contig, or "*" to get unaligned reads
samtools view hits.bam | sort_reads_from_bam.py -i - -c contig123 -1 reads1.fq.gz -2 reads2.fq.gz

    to count mapped reads per contig and not extract any reads:
samtools view hits.bam | sort_reads_from_bam.py -i - > hits_counts.txt

    trimmed reads were originally the same length, so counts are uniform
    but if reads are not the same length (perhaps Sanger), use -b

    to take only reads that did not align:
samtools view hits.bam | sort_reads_from_bam.py -i - -c "*" -1 reads1.fq.gz -2 reads2.fq.gz

    bowtie2 BAM must be generated without --no-al, such as:
bowtie2 -q -p 6 --no-discordant --no-sq -x contigs.fa -1 reads1.fq.gz -2 reads2.fq.gz | samtools view -bT contigs.fa - | samtools sort - hits.sorted
'''

import sys
import os
import argparse
import time
import gzip
from collections import defaultdict
from itertools import izip
from Bio import SeqIO

def get_contig_list(contigfile, delimiter):
	'''parse file of contig names and return a dictionary where keys are contig IDs and values are True'''
	if not os.path.isfile(contigfile):
		if contigfile=="*": # if star is given, then make dict of * only
			sys.stderr.write( "# Searching reads that did not align\n" )
			contigdict = {"*":True}
		else:
			sys.stderr.write( "# Searching reads mapped to {}\n".format(contigfile) )
			contigdict = {contigfile:True}
		return contigdict
	sys.stderr.write( "# Reading contigs file {}  {}\n".format(contigfile, time.asctime() ) )
	if delimiter and delimiter=="\\t": # in case the instructions are disregarded and "\t" is used
		delimiter = "\t"
	contigdict = {}
	example = ""
	for line in open(contigfile,'r'):
		line = line.rstrip()
		if delimiter:
			line = line.split(delimiter)[0]
		if line[0] == ">":
			line = line[1:]
		contigdict[line] = True
		if not example:
			example = line
			sys.stderr.write( "Contigs from file parsed as: {}\n".format(example) )
	sys.stderr.write( "# Found {} contigs  {}\n".format(len(contigdict) , time.asctime() ) )
	return contigdict

def bam_to_read_counts(samfile, nodiscord, doexclude, keepsingle, contigs, bybase, verbose):
	'''parse the SAM file, and count reads mapped to each contig as a dict, and if collecting reads, then also return a dict of read IDs to keep'''
	sys.stderr.write( "# Parsing SAM input  {}\n".format( time.asctime() ) )
	readcount = 0
	keepcount = 0
	readidstokeep = {}
	scaffoldCountsDict = defaultdict(int)
	for line in samfile:
		if line:
			readcount += 1
			if verbose:
				if not readcount % 100000:
					sys.stderr.write(".")
				if not readcount % 1000000:
					sys.stderr.write( "{}   {}\n".format( readcount, time.asctime() ) )
			# Read_12345	0	contig	86	254	100M	*	0	0	TTCTCATGG
			lsplits = line.split("\t")
			if not keepsingle and lsplits[6]=="*": # second read did not match
				continue
			# if bam was generated with --no-discordant this isnt needed
			if nodiscord and not lsplits[6]=="=": # lsplits[6] is contig of paired read
				continue
			#if int(lsplits[1])<256: # assume read is ok
			scaffold = lsplits[2]

			# if reads are not uniform length, must count total bases mapped, and later divide by contig length
			if bybase:
				readlength = len(lsplits[9])
				scaffoldCountsDict[scaffold] += readlength
			else: # otherwise just add 1 for each read
				scaffoldCountsDict[scaffold] += 1

			if contigs:
				if doexclude: # should take all IDs that do not match contigs
					if not contigs.get(scaffold,False): # thereby excluding reads matching those contigs
						keepcount += 1
						readidstokeep[lsplits[0]] = True
				else:
					if contigs.get(scaffold,False): # should take only IDs that match contigs
						keepcount += 1
						readidstokeep[lsplits[0]] = True
	else:
		sys.stderr.write( "{}   {}\n".format( readcount, time.asctime() ) )
	sys.stderr.write( "# Finished parsing {} reads  {}\n".format(readcount, time.asctime() ) )
	if keepcount:
		sys.stderr.write( "Parsed {} reads to keep for {} IDs\n".format(keepcount, len(readidstokeep) ) )
	return readidstokeep, scaffoldCountsDict

def get_read_pairs(reads1, reads2, kept1, kept2, readIDstoKeep, readformat, opentype, verbose):
	'''from the kept read IDs dict, iterate through read files and keep reads that match the dict'''
	readcounter = 0
	writecounter = 0
	with opentype(reads1,'r') as f1, opentype(reads2,'r') as f2, opentype(kept1,'w') as o1, opentype(kept2,'w') as o2:
		for sr1,sr2 in izip(SeqIO.parse(f1, readformat), SeqIO.parse(f2, readformat) ):
			readcounter += 1
			if verbose:
				if not readcounter % 100000:
					sys.stderr.write(".")
				if not readcounter % 1000000:
					sys.stderr.write("{}   {}\n".format( readcounter, time.asctime() ) )
			if readIDstoKeep.get(sr1.id, False) or readIDstoKeep.get(sr2.id, False):
				writecounter += 1
				o1.write(sr1.format(readformat))
				o2.write(sr2.format(readformat))
	sys.stderr.write( "Parsed {} read pairs  {}\n".format(readcounter, time.asctime() ) )
	sys.stderr.write( "Kept {} ({:.2f}%) read pairs\n".format(writecounter, writecounter*100.0/readcounter) )
	# no return

def get_long_reads(longreads, keptreads, readIDstoKeep, readformat, opentype, verbose):
	'''from the kept read IDs dict, iterate through read file and keep reads that match the dict, as with get_read_pairs but for only one file'''
	readcounter = 0
	writecounter = 0
	with opentype(longreads,'r') as f1, opentype(keptreads,'w') as o1:
		for sr in SeqIO.parse(f1, readformat):
			readcounter += 1
			if verbose:
				if not readcounter % 10000:
					sys.stderr.write(".")
				if not readcounter % 100000:
					sys.stderr.write("{}   {}\n".format( readcounter, time.asctime() ) )
			if readIDstoKeep.get(sr.id, False):
				writecounter += 1
				o1.write(sr.format(readformat))
	sys.stderr.write( "Parsed {} reads  {}\n".format(readcounter, time.asctime() ) )
	sys.stderr.write( "Kept {} ({:.2f}%) reads\n".format(writecounter, writecounter*100.0/readcounter ) )
	# no return

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-i','--input', type = argparse.FileType('rU'), default = '-', help="sam or bam file")
	parser.add_argument('-c','--contigs', help="list of contigs, genomic or transcriptomic")
	parser.add_argument('-d','--delimiter', help="delimiter of contig name, such as ',' or ' ' from csv or table")
	parser.add_argument('-1','--file1', help="reads one of a pair")
	parser.add_argument('-2','--file2', help="reads two of a pair")
	parser.add_argument('-l','--long', help="long reads, unpaired, sorted separately from pairs")
	parser.add_argument('-f','--format', default="fastq", help="read format [fastq]")
	parser.add_argument('-o','--outname', default="kept", help="filename tag for kept reads [kept]")
	parser.add_argument('-n','--no-discordant', action="store_true", help="only count aligned pairs, ignoring singletons and discordant alignments")
	parser.add_argument('-b','--bases', action="store_true", help="count coverage as bases mapped, instead of reads")
	parser.add_argument('-e','--exclude', action="store_true", help="remove reads that match contigs, rather than keep them")
	parser.add_argument('-s','--single', action="store_true", help="treat reads as separate, allowing for single reads to match")
	parser.add_argument('-v','--verbose', action="store_true", help="verbose output")
	args = parser.parse_args(argv)

	contigdict = get_contig_list(args.contigs, args.delimiter) if args.contigs else None

	readIDstoKeep, scaffoldCounts = bam_to_read_counts(args.input, args.no_discordant, args.exclude, args.single, contigdict, args.bases, args.verbose)

	# GET PAIRED END READS
	if args.file1 and args.file2:
		f1base, f1ext = os.path.splitext(args.file1)
		f2base, f2ext = os.path.splitext(args.file2)
		if f1ext==".gz" and f2ext==".gz":
			sys.stderr.write( "# Parsing read pairs as {} gz\n".format(args.format, time.asctime() ) )
			outfile1 = "{}.{}.gz".format(f1base, args.outname)
			outfile2 = "{}.{}.gz".format(f2base, args.outname)
			opentype = gzip.open
		else:
			sys.stderr.write( "# Parsing read pairs as {}  {}\n".format(args.format, time.asctime() ) )
			outfile1 = "{}.{}".format(args.file1, args.outname)
			outfile2 = "{}.{}".format(args.file2, args.outname)
			opentype = open
		get_read_pairs(args.file1, args.file2, outfile1, outfile2, readIDstoKeep, args.format, opentype, args.verbose)

	# GET LONG READS
	if args.long:
		lbase, lext = os.path.splitext(args.long)
		if lext==".gz": # for anything like reads1.fastq.gz
			sys.stderr.write( "# Parsing {} as {} gz  {}\n".format(args.long, args.format, time.asctime() ) )
			outfile1 = "{}.{}.gz".format(lbase, args.outname)
			opentype = gzip.open
		else:
			opentype = open
			sys.stderr.write( "# Parsing {} as {}  {}\n".format(args.long, args.format, time.asctime() ) )
			outfile1 = "{}.{}{}".format(lbase, args.outname, lext) # lext should include .
		get_long_reads(args.long, outfile1, readIDstoKeep, args.format, opentype, args.verbose)

	# if not separating reads, then output the scaffold counts
	scafcountsum = 0
	for k,v in sorted(scaffoldCounts.items(), key=lambda x: x[0]):
		if contigdict:
			if args.exclude: # should print scaffolds that do not match contigs
				if not contigdict.get(k,True):
					scafcountsum += v
					wayout.write( "{} {}\n".format(k,v) )
			else:
				if contigdict.get(k,False):
					scafcountsum += v
					wayout.write( "{} {}\n".format(k,v) )
		else: # if contigs are not preselected, then just output only contigs and counts found in the SAM
			wayout.write( "{} {}\n".format(k,v) )
	if scafcountsum:
		sys.stderr.write( "Found {} total mapped reads  {}\n".format(scafcountsum, time.asctime() ) )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
