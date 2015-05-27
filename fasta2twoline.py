#! /usr/bin/env python
#
# fasta2twoline.py 2015-05-27
'''convert padded fasta to two line fasta format:

fasta2twoline.py seqs.fasta > twoline_seqs.fasta
'''

import sys

if len(sys.argv) < 2:
	print __doc__
	sys.exit()
else:
	fastafile = sys.argv[1]
	seq = ""
	print >> sys.stderr, "Reading", fastafile
	for line in open(fastafile,'r'):
		line = line.rstrip()
		if line[0]==">":
			if seq:
				print seq
				seq = ""
			print line
		else:
			seq += line
	print >> sys.stderr, "Done"
