#!/usr/bin/env python
#
# v1.n 2022-12-02
# 2023-03-14 add option for rna hits -r

'''hits_to_coverage.py  last modified 2023-03-20

    combine number of reads mapped to each contig with
    length and GC content, to make a tabular output as:
seqID  seqnumber  length  coverage  GC%  gaps

hits_to_coverage.py -f contigs.fasta -b hits_from_bam.txt > coverage.tab

    to optionally change the read length (default is 100),
    use the -l option

hits_to_coverage.py -f contigs.fasta -b hits_from_bam.txt -l 125 > coverage.tab

    bedgraph input can instead be used with -g
    when using a bedgraph as the input, to infer coverage
    also set -l to 1 as the read length is ignored for bed

hits_to_coverage.py -f contigs.fasta -g contigs.cov.bg -l 1 > contigs.gc_cov.tab
'''

import sys
import time
import argparse
import gzip
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b', '--hits-from-bam', help="hits per contig, counted from the .bam file, instead of -g")
	parser.add_argument('-g', '--bedgraph', help="bedgraph file of genome coverage, instead of -b, can be .gz")
	parser.add_argument('-r', '--rna-hits', help="RNAseq hits per contig, counted from the .bam file, to be appended if available")
	parser.add_argument('-f', '--fasta', help="fasta files of contigs or scaffolds")
	parser.add_argument('-l', '--read-length', type=float, metavar='N', default=100.0, help="read length (must be constant)")
	parser.add_argument('--rna-read-length', type=float, metavar='N', default=100.0, help="secondary read length, used with -r (must be constant)")
	parser.add_argument('-A', '--above', type=int, metavar='N', default=0, help="keep contigs with a count greater than N [0], this removes any contigs with 0 coverage by default")
	parser.add_argument('-B', '--below', type=int, metavar='N', default=10000, help="keep contigs with a count less than N [10000]")
	parser.add_argument('-S', '--strong', type=float, metavar='0.F', default=0.0, help="keep contigs with GC percent greater than N [0.0]")
	parser.add_argument('-W', '--weak', type=float, metavar='0.F', default=100.0, help="keep contigs with GC percent less than N [100.0]")
	args = parser.parse_args(argv)

	dna_hitdict = {}
	rna_hitdict = {}
	multicov = False
	use_dummy = False

	if args.hits_from_bam:
		print("# Reading counts from BAM file {}, assuming read length of {}  {}".format(args.hits_from_bam, args.read_length, time.asctime() ), file=sys.stderr)
		for line in open(args.hits_from_bam,'r'):
			line = line.strip()
			if line and line[0] != "#":
				seqid = line.split(' ')[0]
				hits = line.split(' ')[1]
				dna_hitdict[seqid] = int(hits)
	elif args.bedgraph:
		if args.read_length==100.0: # meaning default
			print("# Using bedgraph, changing assumed read length from {} to 1".format(args.read_length), file=sys.stderr)
			args.read_length = 1.0
		if args.bedgraph.rsplit('.',1)[-1]=="gz": # autodetect gzip format
			opentype = gzip.open
			print("# Inferring counts from gzipped bedgraph file {}  {}".format(args.bedgraph, time.asctime() ), file=sys.stderr)
		else: # otherwise assume normal open for fasta format
			opentype = open
			print("# Inferring counts from bedgraph file {}  {}".format(args.bedgraph, time.asctime() ), file=sys.stderr)
		for line in opentype(args.bedgraph,'rt'):
			line = line.strip()
			lsplits = line.split("\t")
			contigname, startpos, endpos, coverage = lsplits[0:4]
			span = int(endpos)-int(startpos)
			cov_span = span * int(coverage)
			dna_hitdict[contigname] = dna_hitdict.get(contigname,0) + cov_span
	elif args.rna_hits: # just to check if it is there, to not exit
		pass
	else:
		use_dummy = True
		print("# NO input given with -b, or -r, or -g, running in dummy mode", file=sys.stderr)

	if args.rna_hits:
		print("# Reading secondary counts from BAM file {}  {}".format(args.rna_hits, time.asctime() ), file=sys.stderr)
		for line in open(args.rna_hits,'r'):
			seqid, hits = line.strip().split(' ')
			rna_hitdict[seqid] = int(hits)

	seqcounter = 0
	keptcounter = 0
	baseremoval = 0

	length_total = 0
	covbp_total = 0

	if len(dna_hitdict) > 0 or len(rna_hitdict) > 0: # if both are used, write as multicoverage
		multicov = True

	if multicov:
		sys.stdout.write( "scaffold\tnumber\tlength\tcoverage\tGC\tgaps\tcoverage2\n" )
	else:
		sys.stdout.write( "scaffold\tnumber\tlength\tcoverage\tGC\tgaps\n" )

	print( "# Reading contigs from {}  {}".format(args.fasta, time.asctime() ), file=sys.stderr)
	for seqrec in SeqIO.parse(args.fasta,"fasta"):
		seqcounter += 1
		seqlength = len(seqrec.seq)
		length_total += seqlength
		gaps = (seqrec.seq.count("N")+seqrec.seq.count("n"))
		gc = (seqrec.seq.count("G")+seqrec.seq.count("C")) * 100.0 / ( len(seqrec.seq) - gaps )
		if len(dna_hitdict) > 0:
			raw_coverage = int(args.read_length * dna_hitdict.get(seqrec.id, int(use_dummy) ) )
		else: # assume to use RNA dict instead
			raw_coverage = int(args.read_length * rna_hitdict.get(seqrec.id, 0 ) )
		covbp_total += raw_coverage
		contig_coverage = 1.0 * raw_coverage / seqlength
		if contig_coverage <= args.above or contig_coverage > args.below or gc < args.strong or gc > args.weak:
			print( "# IGNORING {}: len={} GC={:.2f} COVERAGE={:.2f}".format( seqrec.id, seqlength, gc, contig_coverage ), file=sys.stderr)
			baseremoval += seqlength
			continue
		keptcounter += 1
		if multicov:
			rna_coverage = 1.0 * args.rna_read_length * rna_hitdict.get(seqrec.id,0) / seqlength
			sys.stdout.write( "{}\t{}\t{}\t{:.1f}\t{}\t{}\t{:.1f}\n".format( seqrec.id, seqcounter, seqlength, contig_coverage, gc, gaps, rna_coverage ) )
		else:
			sys.stdout.write( "{}\t{}\t{}\t{:.1f}\t{}\t{}\n".format( seqrec.id, seqcounter, seqlength, contig_coverage, gc, gaps ) )
	print( "# Counted {} sequences, kept {}  {}".format( seqcounter, keptcounter, time.asctime() ), file=sys.stderr)
	print( "# {} total bp, average coverage is {:.1f}".format( length_total, 1.0*covbp_total/length_total ), file=sys.stderr)
	if baseremoval:
		print( "# Removed {} sequences, totaling {} bases".format( seqcounter-keptcounter, baseremoval ), file=sys.stderr)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
