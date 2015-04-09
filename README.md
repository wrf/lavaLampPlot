# lavaLampPlot
instructions, python and R code for generating lava lamp plots of kmer coverage

## Overview
This includes a few scripts to generate a "lava lamp" plot of kmer coverage separated by GC content of the kmers. This is usually enough to reveal heterozygosity or the presence of symbionts with very different GC content than the host.

An example plot is shown for *Hydra vulgaris* [SRR1032106](http://www.ncbi.nlm.nih.gov/sra/SRX378887) using a kmer of 31. The putative symbiont is visible in the plot with a high GC content.

Steps for analysis and generation of the plots. Some of the instructions were borrowed from Joe Ryan's [estimate genome size](https://github.com/josephryan/estimate_genome_size.pl) script:

## Dependencies
Download the [jellyfish kmer counter](http://www.genome.umd.edu/jellyfish.html) (or another preferred kmer counter, I guess...). If you are using Trinity for transcriptome assembly, then you already have it since the jellyfish binary is supplied with Trinity in the `trinity-plugins/` folder.

## Operation
1. Run jellyfish on the raw genomic data.
   
   `jellyfish count -m 31 -s 1000000000 -C -o fastq.counts -U 1000 -t 8 all_reads.fastq`
  
   data can also be piped, in which case you must provide the input `/dev/fd/0` (or `/dev/stdin`):
   
   `cat *.fastq | jellyfish count -m 31 -s 1000000000 -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
   
   or unzipped in real time:
  
   `gzip -dc *.fastq.gz | jellyfish count -m 31 -s 1000000000 -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
  
   `-s` is the memory buffer, 1 billion is safe for many analyses, though the array will expand if it is maxed out. This can consume a lot of memory (30-50Gb), so it is not advisable to run on a laptop. `-U` is the max coverage to count. 1000 will catch mostly transposons. `-t` is the thread count. `-C` refers to canonical counting, meaning both + and - strands.

2. Write out kmers with coverage counts. This file can be massive (50+Gb) so it is instead piped through the python script. `-k` and `-u` are kmer size and max coverage from jellyfish. This step is a single process, so can take some time for very large files with 1 billion+ kmers (20-30 minutes). This generates a matrix of GC vs coverage across all kmers, where the value at each position is the number of unique kmers with that GC and coverage. This is run through python so that R does not have to deal with the counting in addition to the graphing.
   
   `jellyfish dump fastq.counts | fastqdumps2histo.py -k 31 -u 1000 - > fastq.gc_cov.histo.csv`

3. Configure the R script for this .csv file, and pdf output. Add text to annotate the plot as needed with the  `text(x,y,"important things")` command. The X and Y values correspond to the center position of the X and the GC count (not the percentage). Thus a GC of 25% would be something like 7 or 8 for a kmer of 31. This script can be run interactively in an R environment (such as RKWard or RCommander) or in the command line.

   `Rscript jellyfish_gc_coverage_blob_plot_v2.R`

## Usage considerations
#### Choosing k-mer length
Due to the connection between kmer length and coverage, there is necessarily a balance between longer kmers, which will resolve the y-axis better, and higher coverage, which will resolve the x-axis better. Kmers between 31 and 41 tend to perform fairly well.

#### Sequencing depth
Sequencing coverage is also important. Environmental metagenomes (not amplicon sequencing) of fairly low coverage (5 gigbases) did not display any useful complexity. For metazoan genomes of 100-200Mb, around 20-30Gb of sequence data usually produced high quality plots displaying blobs for both animal and potential bacterial symbionts (see the example plot of Hydra, using one SRA of 25Gb of raw sequence).

#### Memory usage
Jellyfish can easily max out the memory on a system. The hash buffer (`-s`) will expand if jellyfish counts more kmers than the number initially provided. Some systems have safeguards to kill any process demanding 99% of the memory, but otherwise this may crash a computer. In those cases, it is advisable to split the dataset or just run `jellyfish count` on a faction of the reads, such as the left or right reads alone. In both cases, the counts are then combined using the `jellyfish merge` command:

`jellyfish merge -o combined.counts left.counts right.counts`

followed by running `jellyfish dump` and the python script as above:

`jellyfish dump combined.counts | fastqdumps2histo.py -k 31 -u 2000 - > combined.gc_cov.histo.csv`

Note that when using merge, the `-u` option for fastqdumps2histo must be changed accordingly, probably multiplied for each file that is merged. If this was two files with `-U 1000`, then set the value to 2000 rather than 1000.

The manual from jellyfish v1 implies that if the output file is not specified with `-o`, jellyfish should create an extra output file with the default name ("output_") when the hash size is full. I have not tested if this works for sparing memory usage.

## Misc
As this is not really published work, citing is probably not necessary. Nonetheless, it may be advisable to say that any figures were created using this repo, something like "used lavaLampPlot python and R scripts by WRF".
