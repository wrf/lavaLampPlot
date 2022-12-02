#!/usr/bin/env python

'''hits_to_coverage.py  last modified 2022-12-02

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
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b', '--hits-from-bam', help="hits per contig, counted from the .bam file, instead of -g")
	parser.add_argument('-g', '--bedgraph', help="bedgraph file of genome coverage, instead of -b")
	parser.add_argument('-f', '--fasta', help="fasta files of contigs or scaffolds")
	parser.add_argument('-l', '--read-length', type=float, metavar='N', default=100.0, help="read length (must be constant)")
	parser.add_argument('-A', '--above', type=int, metavar='N', default=0, help="keep contigs with a count greater than N [0]")
	parser.add_argument('-B', '--below', type=int, metavar='N', default=10000, help="keep contigs with a count less than N [10000]")
	parser.add_argument('-S', '--strong', type=float, metavar='0.F', default=0.0, help="keep contigs with GC percent greater than N [0.0]")
	parser.add_argument('-W', '--weak', type=float, metavar='0.F', default=100.0, help="keep contigs with GC percent less than N [100.0]")
	args = parser.parse_args(argv)

	hitdict = {}
	if args.hits_from_bam:
		print("# Reading counts from BAM file {}  {}".format(args.hits_from_bam, time.asctime() ), file=sys.stderr)
		for line in open(args.hits_from_bam,'r'):
			seqid, hits = line.strip().split(' ')
			hitdict[seqid] = int(hits)
	elif args.bedgraph:
		print("# Inferring counts from bedgraph file {}  {}".format(args.bedgraph, time.asctime() ), file=sys.stderr)
		for line in open(args.bedgraph,'r'):
			line = line.strip()
			lsplits = line.split("\t")
			contigname, startpos, endpos, coverage = lsplits[0:4]
			span = int(endpos)-int(startpos)
			cov_span = span * int(coverage)
			hitdict[contigname] = hitdict.get(contigname,0) + cov_span
	else:
		sys.exit("ERROR: no input given with -b or -g, exiting")

	seqcounter = 0
	keptcounter = 0
	baseremoval = 0

	length_total = 0
	covbp_total = 0

	print( "# Reading contigs from {}  {}".format(args.fasta, time.asctime() ), file=sys.stderr)
	sys.stdout.write( "scaffold\tnumber\tlength\tcoverage\tGC\tgaps\n" )
	for seqrec in SeqIO.parse(args.fasta,"fasta"):
		seqcounter += 1
		seqlength = len(seqrec.seq)
		length_total += seqlength
		gaps = (seqrec.seq.count("N")+seqrec.seq.count("n"))
		gc = (seqrec.seq.count("G")+seqrec.seq.count("C")) * 100.0 / ( len(seqrec.seq) - gaps )
		raw_coverage = int(args.read_length * hitdict.get(seqrec.id,0) )
		covbp_total += raw_coverage
		contig_coverage = 1.0 * raw_coverage / seqlength
		if contig_coverage <= args.above or contig_coverage > args.below or gc < args.strong or gc > args.weak:
			print( "# IGNORING {}: len={} GC={:.2f} COVERAGE={:.2f}".format( seqrec.id, seqlength, gc, contig_coverage ), file=sys.stderr)
			baseremoval += seqlength
			continue
		keptcounter += 1
		sys.stdout.write( "{}\t{}\t{}\t{:.1f}\t{}\t{}\n".format( seqrec.id, seqcounter, seqlength, contig_coverage, gc, gaps ) )
	print( "# Counted {} sequences, kept {}  {}".format( seqcounter, keptcounter, time.asctime() ), file=sys.stderr)
	print( "# {} total bp, average coverage is {:.1f}".format( length_total, 1.0*covbp_total/length_total ), file=sys.stderr)
	if baseremoval:
		print( "# Removed {} sequences, totaling {} bases".format( seqcounter-keptcounter, baseremoval ), file=sys.stderr)

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
