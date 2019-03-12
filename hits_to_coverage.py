#!/usr/bin/env python

'''hits_to_coverage.py  last modified 2019-03-11

    combine number of reads mapped to each contig with
    length and GC content, to make a tabular output as:
seqID  seqnumber  length  coverage  GC%  gaps

hits_to_coverage.py -f contigs.fasta -b hits_from_bam.txt > coverage.tab

    to optionally change the read length (default is 100),
    use the -l option

hits_to_coverage.py -f contigs.fasta -b hits_from_bam.txt -l 125 > coverage.tab
'''

import sys
import argparse
from Bio import SeqIO

def main(argv, wayout):
	if not len(argv):
		argv.append("-h")
	parser = argparse.ArgumentParser(formatter_class=argparse.RawDescriptionHelpFormatter, description=__doc__)
	parser.add_argument('-b', '--hits-from-bam', help="hits per contig, counted from the .bam file")
	parser.add_argument('-f', '--fasta', help="fasta files of contigs or scaffolds")
	parser.add_argument('-l', '--read-length', type=float, metavar='N', default=100.0, help="read length (must be constant)")
	parser.add_argument('-A', '--above', type=int, metavar='N', default=0, help="keep contigs with a count greater than N [0]")
	parser.add_argument('-B', '--below', type=int, metavar='N', default=10000, help="keep contigs with a count less than N [10000]")
	parser.add_argument('-S', '--strong', type=float, metavar='0.F', default=0.0, help="keep contigs with GC percent greater than N [0.0]")
	parser.add_argument('-W', '--weak', type=float, metavar='0.F', default=100.0, help="keep contigs with GC percent less than N [100.0]")
	args = parser.parse_args(argv)

	hitdict = {}
	print >> sys.stderr, "# Reading counts from {}".format(args.hits_from_bam)
	for line in open(args.hits_from_bam,'r'):
		seqid, hits = line.strip().split(' ')
		hitdict[seqid] = int(hits)

	seqcounter = 0
	keptcounter = 0
	baseremoval = 0
	print >> sys.stderr, "# Reading contigs from {}".format(args.fasta)
	print >> sys.stdout, "scaffold\tnumber\tlength\tcoverage\tGC\tgaps"
	for seqrec in SeqIO.parse(args.fasta,"fasta"):
		seqcounter += 1
		seqlength = len(seqrec.seq)
		gaps = (seqrec.seq.count("N")+seqrec.seq.count("n"))
		gc = (seqrec.seq.count("G")+seqrec.seq.count("C")) * 100.0 / ( len(seqrec.seq) - gaps )
		coverage = int(args.read_length * hitdict.get(seqrec.id,0) / seqlength)
		if coverage <= args.above or coverage > args.below or gc < args.strong or gc > args.weak:
			print >> sys.stderr, "# IGNORING {}: len={} GC={:.2f} COVERAGE={:.2f}".format( seqrec.id, seqlength, gc, coverage )
			baseremoval += seqlength
			continue
		keptcounter += 1
		print >> sys.stdout, "{}\t{}\t{}\t{}\t{}\t{}".format( seqrec.id, seqcounter, seqlength, coverage, gc, gaps )
	print >> sys.stderr, "# Counted {} sequences, kept {}".format( seqcounter, keptcounter )
	if baseremoval:
		print >> sys.stderr, "# Removed {} sequences, totaling {} bases".format( seqcounter-keptcounter, baseremoval )

if __name__ == "__main__":
	main(sys.argv[1:],sys.stdout)
