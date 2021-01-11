#!/usr/bin/env python
#
# add_mgrast_annotation_to_stats.py

'''add_mgrast_annotation_to_stats.py  last modified 2021-01-11
    add MG-RAST blast annotations to the existing GC-cov stats table

    download MG-RAST annotation table (as tab-delimited but named .csv)
    from the "Annotation Downloads" section
    selecting "organism" for "Annotation Type"
    and "RefSeq" for "Data Source"

add_mgrast_annotation_to_stats.py mg_annotation.csv scaffolds.stats.tab > scaffolds.stats_w_clade.tab

'''

data_format = '''
query sequence id	hit m5nr id (md5sum)	percentage identity	alignment length	number of mismatches	number of gap openings	query start	query end	hit start	hit end	e-value	bit score	semicolon separated list of annotations
mgm4916968.3|NODE_150_length_72604_cov_5.527009_16410_18462_-|RefSeq	00007f2b56700f4d0a6b02136ba44a5a	63.56	118	43		209	326	230	347	3.2e-35	146.0	[Clostridium botulinum A str. ATCC 19397];[Clostridium botulinum A str. Hall];[Clostridium botulinum A str. ATCC 3502]
mgm4916968.3|NODE_412_length_35560_cov_5.120406_33316_33880_-|RefSeq	00019f9bf8aac34f97a3c0dadff0607c	67.95	78	25		53	130	149	226	1.0e-23	108.0	[Planctomyces brasiliensis DSM 5305]
'''

import sys
from collections import defaultdict,Counter

if len(sys.argv)<2:
	sys.exit(__doc__)
else:
	seq_to_taxID = defaultdict( lambda: defaultdict(int) ) # key is seq ID
	linecounter = 0

	sys.stderr.write("# Reading annotations from {}\n".format( sys.argv[1] ) )
	for line in open(sys.argv[1],'r'):
		line = line.strip()
		if line and line[0]!="#":
			#                0           1                         2                   3                      4                        5                   6               7         8            9         10          11                12
			# query sequence id    hit m5nr id (md5sum)    percentage identity    alignment length    number of mismatches    number of gap openings    query start    query end    hit start    hit end    e-value    bit score    semicolon separated list of annotations
			lsplits = line.split("\t")
			raw_query_id = lsplits[0]
			if raw_query_id=="query sequence id" or line.find("Download")==0:
				continue
			linecounter += 1
			query_id = raw_query_id.split("|")[1].rsplit("_",3)[0]
			# extract annotations from semicolon delimited list
			annotation_list = lsplits[12].replace("[","").replace("]","").split(";")

			# convert to genus
			genus_list = []
			for li in annotation_list:
				if li.split(" ")[0]=="Candidatus": # catch for Candidatus species
					genus_list.append(li.split(" ")[1])
				else: # otherwise assume genus is first word
					genus_list.append(li.split(" ")[0])
			# count occurrence
			genus_counter = Counter(genus_list)
			for g,c in genus_counter.items():
				seq_to_taxID[query_id][g] += 1
	sys.stderr.write("# Found {} lines with {} annotated scaffolds\n".format( linecounter, len(seq_to_taxID) ) )

	scaffold_count = 0
	sys.stderr.write("# Reading tabular stats from {}\n".format( sys.argv[2] ) )
	for line in open(sys.argv[2]):
		line = line.strip()
		if line:
			lsplits = line.split("\t")
			scaffold_id = lsplits[0]
			if scaffold_id=="scaffold": # means header line
				sys.stdout.write("{}\t{}\n".format( line, "top_annotation" ) )
				continue
			scaffold_count += 1
			annotation_counter = seq_to_taxID.get(scaffold_id,None)
			if annotation_counter is None:
				line_w_annotation = "{}\tUnknown\n".format( line )
			else:
				top_annotation = sorted(annotation_counter.items(), key=lambda x: x[1], reverse=True)[0][0]
				line_w_annotation = "{}\t{}\n".format( line, top_annotation )
			sys.stdout.write(line_w_annotation)
	sys.stderr.write("# Found stats for {} scaffolds\n".format( scaffold_count ) )

#
