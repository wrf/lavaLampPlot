# lavaLampPlot
instructions, python and R code for generating lava lamp plots of kmer coverage

## Overview
This includes a few scripts to generate a "lava lamp" plot of kmer coverage separated by GC content of the kmers. This is usually enough to reveal heterozygosity or the presence of symbionts with very different GC content than the host.

Two example plots are shown for *Hydra vulgaris* [SRR1032106](http://www.ncbi.nlm.nih.gov/sra/SRX378887) using a kmer of 31. The putative symbiont is visible in the plot with a high GC content. The second plot counts based on GC of the full read and median coverage kmers in that read, rather than the kmers alone.

Steps for analysis and generation of the plots. Some of the instructions were borrowed from Joe Ryan's [estimate genome size](https://github.com/josephryan/estimate_genome_size.pl) script:

## Dependencies
Download the [jellyfish kmer counter](http://www.genome.umd.edu/jellyfish.html) (or another preferred kmer counter, I guess...). If you are using Trinity for transcriptome assembly, then you already have it since the jellyfish binary is supplied with Trinity in the `trinity-plugins/` folder.

## Operation
1. Run jellyfish on the raw genomic data.
   
   `jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 all_reads.fastq`
  
   data can also be piped, in which case you must provide the input `/dev/fd/0` (or `/dev/stdin`):
   
   `cat *.fastq | jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
   
   or unzipped in real time:
  
   `gzip -dc *.fastq.gz | jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
  
   `-s` is the memory buffer, 2G (2 billion, or 2000000000) is safe for many analyses, though the array will expand if it is maxed out. This can consume a lot of memory (30-50Gb), so it is not advisable to run on a laptop. `-U` is the max coverage to count. 1000 will catch mostly transposons. `-t` is the thread count. `-C` refers to canonical counting, meaning both + and - strands. 

   Optionally, a lower limit can be set with `-L`. For generating kmer plots (step 3 below) this would affect the data and not be advised, but for the read plots (steps 4 and 5), the Trinity program fastaToKmerCoverageStats assumes that any values below 2 are 1. Since the bulk of the kmers have a count of 1 (66% in the *Hydra* sample), this makes a much smaller file and will run faster at subsequent steps.

2. Write out kmers with coverage counts. This file can be massive (easily 50Gb+, 83G for *Hydra* SRR1032106) so it is instead piped through the python script. `-k` and `-u` are kmer size and max coverage from jellyfish. This step is a single process, so can take some time for very large files with 1 billion+ kmers (20-30 minutes). This generates a matrix of GC vs coverage across all kmers, where the value at each position is the number of unique kmers with that GC and coverage. This is run through python so that R does not have to deal with the counting in addition to the graphing.
   
   `jellyfish dump fastq.counts | fastqdumps2histo.py -k 31 -u 1000 - > fastq.gc_cov.histo.csv`

3. Configure the R script for this .csv file, and pdf output. Add text to annotate the plot as needed with the  `text(x,y,"important things")` command. The X and Y values correspond to the center position of the X and the GC count (not the percentage). Thus a GC of 25% would be something like 7 or 8 for a kmer of 31. This script can be run interactively in an R environment (such as RKWard or RCommander) or in the command line.

   `Rscript jellyfish_gc_coverage_blob_plot_v2.R`

4. Slice out sections of reads based on kmer coverage. This requires finally generating the jellyfish dump file, so that will be redone first. Then `-T` (Trinity mode) is called with kmersorter. This generates a table of median coverage per read, which would be used in Trinity to randomly normalize the data based on coverage. Instead, here it uses the same method to deterministically keep reads with coverage above or below some specified value (`-a` and `-b`, respectively). Specify the kmer used with `-k`. The Trinity base directory is given with `-D`. Some of the subprocesses can use multithreads, which is specified by `-p`. The sorting step requires a memory limit, which is `-m`; note that the most memory intensive step, fastaToKmerCoverageStats, may exceed this limit as it reads the entire fastq.dumps file into memory.

   `jellyfish dump fastq.counts > fastq.dumps`
   
   `kmersorter.py -T -k 31 -m 20G -p 8 -a 100 -b 200 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.dumps`

5. Generate a more precise coverage to GC map using the entire read rather than kmers. This is run similarly as before with some alternate options in fastqdumps2histo. As above, Trinity mode is specified with `-T`. The intermediate stats files from kmersorter are then used with the raw reads to count the coverage and GC. The `-k` value here is the length of the reads, not the kmer length.

   `fastqdumps2histo.py -s reads_1.stats reads_2.stats -f reads_1.fq reads_2.fq -k 100 -u 1000 -T - > reads.gc_cov.histo.csv`

6. Run the R script on this .csv file as above.

## Usage considerations
#### Choosing k-mer length
Due to the connection between kmer length and coverage, there is necessarily a balance between longer kmers, which will resolve the y-axis better, and higher coverage, which will resolve the x-axis better. Kmers between 31 and 41 tend to perform fairly well.

#### Sequencing depth
Sequencing coverage is also important. Environmental metagenomes (not amplicon sequencing) of fairly low coverage (5 gigbases) did not display any useful complexity. For metazoan genomes of 100-200Mb, around 20-30Gb of sequence data usually produced high quality plots displaying blobs for both animal and potential bacterial symbionts (see the example plot of *Hydra*, using one SRA of 25Gb of raw sequence).

#### Memory usage
Jellyfish can easily max out the memory on a system. For a hash size of 1G (1 billion kmers), for a k-mer of 31, it takes 6Gb of memory on my system; 2G takes 12Gb and 3-4G takes 23Gb (probably scales binarily). The hash buffer (`-s`) will expand if jellyfish counts more kmers than the number initially provided. Some systems have safeguards to kill any process demanding 99% of the memory, but otherwise this may crash a computer. In those cases, it is advisable to split the dataset or just run `jellyfish count` on a faction of the reads, such as the left or right reads alone. In both cases, the counts are then combined using the `jellyfish merge` command:

`jellyfish merge -o combined.counts left.counts right.counts`

followed by running `jellyfish dump` and the python script as above:

`jellyfish dump combined.counts | fastqdumps2histo.py -k 31 -u 2000 - > combined.gc_cov.histo.csv`

Note that when using merge, the `-u` option for fastqdumps2histo must be changed accordingly, probably multiplied for each file that is merged. If this was two files with `-U 1000`, then set the value to 2000 rather than 1000. Also, the hash size (`-s`) must be the same for all files in the merge, so pick something large enough, maybe 4G (4 billion) if you have the memory.

The merging generates counts that are slightly different from counting them as a single set, on the order of 0.001% *more* when separate sets were merged than when counted as a single set. This is after excluding those above the counts cutoff for the single set (which might account for 0.01%). I have no explanation for this.

The manual from jellyfish v1 implies that if the output file is not specified with `-o`, jellyfish should create an extra output file with the default name ("output_") when the hash size is full. I have not tested if this works for sparing memory usage.

#### SRA Files
The Trinity mode in kmersorter has a problem with the SRA headers generated from the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/) fastq-dump (due to problems in the called Trinity scripts). The stats file generation will work normally, but then the stats cannot be sorted correctly and ultimately the process will fail merging the stats. Brian Haas suggested this alternate command to generate the header (which one would also use during Trinity normalization).

`fastq-dump --split-files --defline-seq '@$sn[_$rn]/$ri' SRR1032106.sra`

## Misc
As this is not really published work, citing is probably not necessary. Nonetheless, it may be advisable to say that any figures were created using this repo, something like "used lavaLampPlot python and R scripts by WRF".
