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
  
   data can also be piped:
   
   `cat *.fastq | jellyfish count -m 31 -s 1000000000 -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
   
   or unzipped:
  
   `gzip -dc *.fastq.gz | jellyfish count -m 31 -s 1000000000 -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
  
   `-s` is the memory buffer, 1 billion is safe for many analyses, though the array will expand if it is maxed out. This can consume a lot of memory (30-50Gb), so it is not advisable to run on a laptop. `-U` is the max coverage to count. 1000 will catch mostly transposons. `-t` is the thread count. `-C` refers to canonical counting, meaning both + and - strands.

2. Write out kmers with coverage counts. This file can be massive (50+Gb) so it is instead piped through the python script. `-k` and `-u` are kmer size and max coverage from jellyfish. This step is a single process, so can take some time for very large files with 1 billion+ kmers (20-30 minutes). This generates a matrix of GC vs coverage across all kmers, where the value at each position is the number of unique kmers with that GC and coverage. This is run through python so that R does not have to deal with the counting in addition to the graphing.
   
   `jellyfish dump fastq.counts | fastqdumps2histo.py -k 31 -u 1000 - > fastq.gc_cov.histo.csv`

3. Configure the R script for this .csv file, and pdf output. Add text to annotate the plot as needed with the  `text(x,y,"important things")` command. The X and Y values correspond to the center position of the X and the GC count (not the percentage). Thus a GC of 25% would be something like 7 or 8 for a kmer of 31.

## Usage considerations
Due to the connection between kmer length and coverage, there is necessarily a balance between longer kmers, which will resolve the y-axis better, and higher coverage, which will resolve the x-axis better. Kmers between 31 and 41 tend to perform fairly well.

## Misc
As this is not really published work, citing is probably not necessary. Nonetheless, it may be advisable to say that any figures were created using this repo, something like "used lavaLampPlot python and R scripts by WRF".
