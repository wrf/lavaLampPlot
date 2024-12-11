#!/usr/bin/env python
# merge_uc_clusters.py created 2022-12-07
# allow gzip 2023-04-11
# check for write formats from R that cause problems 2024-12-11

"""merge_uc_clusters.py  last modified 2024-12-11
    combine values of OTU counts based on cluster output from vsearch

./merge_uc_clusters.py cluster.uc otu_table.tab > merged_table.tab

    assumes OTU table is in tab-delimited format, as:
seq_name    sample1    sample2    sample3
seq1234     0          10         100
seq5678     50         5          0

    assumes tab-delimited cluster format from vsearch, as:
S	0	503	*	*	*	*	*	ASV1210682	*
S	1	503	*	*	*	*	*	ASV1218341	*
S	2	503	*	*	*	*	*	ASV1221110	*

    after running vsearch:
~/vsearch-2.22.1-linux-x86_64-static/bin/vsearch --cluster_fast dna-sequences.fasta --id 0.97 --centroids centroids.fas --uc clusters.uc
"""

import sys
import gzip
from collections import defaultdict

if len(sys.argv) < 2:
	sys.exit(__doc__)
else:
	cluster_file = sys.argv[1]
	# 1. Record type: S, H, or C. Each fasta sequence in the input file can be either a cluster centroid (S) or a hit (H) assigned to a cluster. Cluster records (C) summarize information (size, centroid label) for each cluster.
	# 2. Cluster number (zero-based).
	# 3. Sequence length (S, H), or cluster size (C).
	# 4. Percentage of similarity with the centroid sequence (H), or set to ’*’ (S,C).
	# 5. Match orientation + or - (H), or set to ’*’ (S, C).
	# 6. Not used, always set to ’*’ (S, C) or 0 (H).
	# 7. Not used, always set to ’*’ (S, C) or 0 (H).
	# 8. Not used, always set to ’*’.
	# 9. Label of the query sequence (H), or of the centroid sequence (S, C).
	# 10. Label of the centroid sequence (H), or set to ’*’ (S, C).

	# S	3	454	*	*	*	*	*	64814cf8dbb546e9c8b75f122ad7cdc0	*
	# H	3	453	99.8	+	0	0	453MI	6564221144419fdc7114c189a9b11116	64814cf8dbb546e9c8b75f122ad7cdc0

	seqid_to_centroid = {} # key is seqid, value is centroid id
	sys.stderr.write("# Reading cluster data from {}\n".format(cluster_file) )
	for line in open(cluster_file,'r'):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			record_type = lsplits[0]
			if record_type == "S": # centroid sequence
				seqid_to_centroid[ lsplits[8] ] = lsplits[8]
			elif record_type == "H": # hit to a cluster
				seqid_to_centroid[ lsplits[8] ] = lsplits[9]
			else: # implicitly "C" meaning cluster record
				pass
	sys.stderr.write("# Counted {} seqs from cluster data\n".format( len(seqid_to_centroid) ) )

	line_counter = 0
	counts_table_file = sys.argv[2]
	merged_counts_dict = {}

	if counts_table_file.rsplit('.',1)[-1]=="gz": # autodetect gzip format
		_opentype = gzip.open
		sys.stderr.write("# Reading counts table from {} as gzipped\n".format(counts_table_file) )
	else: # otherwise assume normal open
		_opentype = open
		sys.stderr.write("# Reading counts table from {}\n".format(counts_table_file) )

	for line in _opentype(counts_table_file,'rt'):
		line = line.rstrip() # use rstrip() in case of initial empty tab from R
		if line:
			line_counter += 1
			lsplits = line.split("\t")
			if line_counter==1: # first line
				n_samples = len(lsplits[1:])
				if lsplits[0]=="": # first column is empty, give name
				    line = "sampleID" + line
				sys.stderr.write("# Detected {} samples\n".format(n_samples) )
				print(line, file=sys.stdout) # use print to add back rstripped return
				continue
			seqid = lsplits[0] # seq name should always be first column
			cluster_target = seqid_to_centroid.get(seqid, None)
			if cluster_target is None:
				sys.stderr.write("# WARNING: no cluster found for {}\n".format(seqid) )
				cluster_target = seqid
			# for each value in the line, convert to float
			seq_counts = [ float(x) for x in lsplits[1:] ]
			# if that cluster already exists, add to existing values
			if cluster_target in merged_counts_dict:
				merged_counts_dict[cluster_target] = [ i+j for i,j in zip(seq_counts, merged_counts_dict.get(cluster_target) ) ]
			else: # otherwise create cluster
				merged_counts_dict[cluster_target] = seq_counts
	sys.stderr.write("# Read {} lines from {}\n".format(line_counter, counts_table_file) )

	# write new clusters
	for k,v in merged_counts_dict.items():
		sys.stdout.write("{}\t{}\n".format(k, "\t".join(list(map(str,v))) ) )
	sys.stderr.write("# Wrote {} clusters\n".format( len(merged_counts_dict) ) )

#
