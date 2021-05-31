# lavaLampPlot #
Instructions, python and R code for generating lava lamp plots of kmer coverage as the kmers themselves, or based on kmer-coverage averages on the raw reads.

This tool was used in the analysis of the [*Hoilungia hongkongensis* genome](https://bitbucket.org/molpalmuc/hoilungia-genome/src/master/) ([Figure S2](https://doi.org/10.1371/journal.pbio.2005359.s002)). It is unlikely I will make a dedicated paper, so if these scripts are used, please cite the paper:

Eitel, M., Francis, W.R., Varoqueaux, F., Daraspe, J., Osigus, H-J., Krebs, S., Vargas, S., Blum, H., Williams, G.A., Schierwater, B., Wörheide, G. (2018) [Comparative genomics and the nature of placozoan species](https://doi.org/10.1371/journal.pbio.2005359). PLoS Biology 16(7):e2005359.

### Contents: ###

* [Dependencies](https://github.com/wrf/lavaLampPlot#dependencies)
* [Operation](https://github.com/wrf/lavaLampPlot#operation)
* [Blob plots](https://github.com/wrf/lavaLampPlot#making-a-blob-plot-of-contigs)
* [Shiny app](https://github.com/wrf/lavaLampPlot#shinyapp-test)
* [Usage tips](https://github.com/wrf/lavaLampPlot#usage-considerations)
* [Troubleshooting](https://github.com/wrf/lavaLampPlot#troubleshooting)

## Overview
This includes a few scripts to generate a "lava lamp" plot of kmer coverage separated by GC content of the kmers. This is a heat-map of frequency of reads or kmers with a particular coverage and GC content, giving a sense of the coverage of the target species and any symbionts or contaminants.

Such plots were intended for use on metagenomic data ([whole-genome shotgun](https://www.ncbi.nlm.nih.gov/genbank/metagenome/), not amplicon), and can be generated with raw data, meaning that no assembly or other processing of the reads is needed. This is usually sufficient to reveal heterozygosity (for eukaryotes) or the presence of symbionts with very different GC content than the host, and instruct whether to even use the data or not. For instance, distinct "spots" suggest an good extraction, while a smear would suggest that there were some biases in the extraction or sequencing and that assembly may be highly problematic.

Below is an example plot from *Hydra vulgaris* [SRR1032106](http://www.ncbi.nlm.nih.gov/sra/SRX378887) using a kmer of 31. The putative symbiont is visible in the plot with a high GC content. The second plot counts based on GC of the full read and median coverage kmers in that read, rather than the kmers alone.

![hydra_vulgaris_SRR1032106_k31_u300_reads.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/hydra_vulgaris_SRR1032106_k31_u300_reads.png)

This is conceptually similar to what is done in the "blob" plots by [blobtools](https://github.com/DRL/blobtools), (paper by [Kumar et al 2013](https://doi.org/10.3389/fgene.2013.00237)) though that process makes use of assembled contigs while this one considers the raw reads directly. However, instructions for generating such plots are also [below](https://github.com/wrf/lavaLampPlot#making-a-blob-plot-of-contigs). Some assemblers output contig coverage information with the contigs, and these can be used [without mapping the reads](https://github.com/wrf/lavaLampPlot#as-a-direct-output-of-some-assemblers).

![hydra_dovetail_w_meta_SRR1032106.gc_cov.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/hydra_dovetail_w_meta_SRR1032106.gc_cov.png)

When making such plots for metagenomics, be aware of the implications of contig coverage for many potential analyses. The coverage is calculated by reads per base, but for the blobs, should reflect the relative number of copies of contigs or chromosomes. For something like a bacterial sample, this could be considered to be a relative counts of cells.

![metagenome_gc_cov_schematic_v2.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/metagenome_gc_cov_schematic_v2.png)

Mapping the original reads back to the assembled contigs should make a mostly even coverage within each scaffold. The raw total counts per gene or per contig will vary by the length, and therefore must be normalized in order to be meaningful (such normalization is already expressed by the coverage). For example, if two bacterial species had equal abundance, but one has a genome 2x larger than the other, this will account for 2x the number of reads, even if the coverage is the same.

![raw_vs_normalized_metagenome_counts_v1.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/raw_vs_normalized_metagenome_counts_v1.png)

This means longer genes, contigs, or chromosomes will have more counts, but also means that total counts to homologous operons between species will vary by both coverage and length. If the objective is to make a statement about the abundance of particular genes or pathways in the sample, the contig coverage could be simplistically applied to the operon. However, it may be better to quantify the actual expression (qPCR or metatranscriptomics) or protein activity instead.

![metagenome_gc_cov_vs_read_counts_v1.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/metagenome_gc_cov_vs_read_counts_v1.png)

## Dependencies
Download the [jellyfish kmer counter](https://github.com/gmarcais/Jellyfish) (or another preferred kmer counter). This used to be supplied with [Trinity for transcriptome assembly](https://github.com/trinityrnaseq/trinityrnaseq/wiki), and was previously in the `trinity-plugins/` folder.

## Operation
1. Run jellyfish on the raw genomic data, i.e. paired end reads.
   
   `jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 all_reads.fastq`
  
   data can also be piped, in which case you must provide the input `/dev/fd/0` (or `/dev/stdin`):
   
   `cat *.fastq | jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
   
   or unzipped in real time:
  
   `gzip -dc *.fastq.gz | jellyfish count -m 31 -s 2G -C -o fastq.counts -U 1000 -t 8 /dev/fd/0`
  
   * `-s` is the memory buffer, 2G (2 billion, or 2000000000) is safe for many analyses, though the array will expand if it is maxed out. This can consume a lot of memory (30-50Gb), so it is not advisable to run on a laptop.
   * `-U` is the max coverage to count. 1000 will catch mostly transposons.
   * `-t` is the thread count / number of CPUs.
   * `-C` refers to canonical counting, meaning both + and - strands.
   * `-m` is the kmer size. If downstream steps involving Trinity will be used (Step 4 and after) then do not select a kmer greater than 32, as this is the upper limit for the `fastaToKmerCoverageStats` program.

   Optionally, a lower limit can be set with `-L`. For generating kmer plots (step 3 below) this would affect the data and not be advised, but for the read plots (steps 4 and 5), the Trinity program `fastaToKmerCoverageStats` assumes that any values below 2 are 1. Since the bulk of the kmers have a count of 1 (66% in the *Hydra* sample), this makes a much smaller file and will run faster at subsequent steps.

2. Write out kmers with coverage counts. This file can be massive (easily 50Gb+, 83G for *Hydra* SRR1032106) so it is instead piped through the python script. `-k` and `-u` are kmer size and max coverage from jellyfish. This step is a single process, so can take some time for very large files with 1 billion+ kmers (20-30 minutes). This generates a matrix of GC vs coverage across all kmers, where the value at each position is the number of unique kmers with that GC and coverage. This is run through python so that R does not have to deal with the counting in addition to the graphing.
   
   `jellyfish dump fastq.counts | fastqdumps2histo.py -k 31 -u 1000 - > fastq.gc_cov.histo.csv`

   If the Trinity mode steps (4 and 5) will eventually be carried out, then the jellyfish kmer dump is required. A filtered version where kmers with a count of 1 have been removed can be generated using the `-j` option.
   
   `jellyfish dump fastq.counts | fastqdumps2histo.py -k 31 -u 1000 -j fastq.dumps - > fastq.gc_cov.histo.csv`

3. Configure the R script for this .csv file, and pdf output. Add text to annotate the plot as needed with the  `text(x,y,"important things")` command. The X and Y values correspond to the center position of the X and the GC count (not the percentage). Thus a GC of 25% would be something like 7 or 8 for a kmer of 31. This script can be run interactively in an R environment (such as RKWard or RCommander) or in the command line. For the command line, an optional 3rd argument can be given to change the X-axis scale (default is 1000).

   `Rscript jellyfish_gc_coverage_blob_plot_v2.R fastq.gc_cov.histo.csv counts_k31_gc_cov.pdf`

4. Slice out sections of reads based on kmer coverage. This requires finally generating the jellyfish dump file, so that will be redone first if it was not done at the earlier step. 

   `jellyfish dump fastq.counts > fastq.dumps`

   Then `-T` (Trinity mode) is called with kmersorter. This generates a table of median coverage per read, which would be used in Trinity to randomly normalize the data based on coverage. Instead, here it uses the same method to deterministically keep reads with coverage above or below some specified value (`-a` and `-b`, respectively). Specify the kmer used with `-k`. The Trinity base directory is given with `-D`. Some of the subprocesses can use multithreads, which is specified by `-p`. Note that the most memory intensive step, fastaToKmerCoverageStats, has no intrinsic memory limit as it reads the entire fastq.dumps file into memory; this is why the reduced dump is generated with the `-j` option in step 2.

   `kmersorter.py -T -k 31 -p 8 -a 100 -b 200 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.dumps`

5. Generate a more precise coverage to GC map using the entire read rather than kmers. This is run similarly as before with some alternate options in fastqdumps2histo. As above, Trinity mode is specified with `-T`. The intermediate stats files from kmersorter are then used with the raw reads to count the coverage and GC. The read length should be detected automatically, though can be specified with the `-r` option.

   `fastqdumps2histo.py -s reads_1.stats reads_2.stats -f reads_1.fq reads_2.fq -u 1000 -T - > reads.gc_cov.histo.csv`

6. Run the R script on this .csv file as above.

7. With the stats files already generated from step 5, cut out sections of reads from the plot by running kmersorter again with different options. Use `-s` and `-w` to set the lower and upper bounds of GC percent cutoff. For example:

   `kmersorter.py -T -k 31 -p 8 -a 50 -b 100 -s 0.5 -w 0.7 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.dumps`
   
   `kmersorter.py -T -k 31 -p 8 -a 45 -b 80 -s 0.35 -w 0.45 -D ~/trinityrnaseq/ -1 reads_1.fq -2 reads_2.fq fastq.dumps`

8. To count kmer coverage of long reads (such as Sanger, PacBio, Moleculo, etc.) or contigs, use the `-L` option in kmersorter. This will calculate the median kmer coverage of the Illumina reads on each long read, that is, even though the coverage might be lower in the long read library, individual spots should become even more distinct on the plot because they are using the Illumina kmer coverage. In this example, the `-S` option is used to only calculate the stats at this step, however this can be excluded and used with filtering options. The original paired reads can be ignored for this step.

   `kmersorter.py -T -k 31 -p 8 -S -D ~/trinityrnaseq/ -L long_reads.fq fastq.dumps`
   
   After stats are generated, the filter options can be applied as above to select groups of reads by coverage or GC. 

   `kmersorter.py -T -k 31 -p 8 -a 45 -b 80 -s 0.35 -w 0.45 -D ~/trinityrnaseq/ -L long_reads.fq fastq.dumps`

9. Generate a long-read coverage to GC map, similar to the above step, except with the addition of the `-p` option. This option indicates that reads are not a standard length, and percentage GC should be calculated for each read.

   `fastqdumps2histo.py -s long_reads.stats -f long_reads.fq -u 1000 -T -p - > long_reads.gc_cov.histo.csv`
   
   This script assumes that fastq files are always 4 lines and fasta files are always 2 lines. If your fasta format files (say for assembled contigs) are more than 2 lines per sequence (which would be header line and sequence line), then convert them first using the fasta2twoline.py script.
   
   `fasta2twoline.py contigs.fasta > twoline_contigs.fasta`

## Making a blob-plot of contigs
The above steps 8 and 9 can be applied to assembled contigs/scaffolds, however, it may be conceptually easier to generate something like a [blob-plot](https://github.com/DRL/blobtools) to view the contigs and/or contamination directly.

![twilhelma_2014_vs_scaffolds_v1.coverage.png](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/twilhelma_2014_vs_scaffolds_v1.coverage.png)

1. Map reads with a normal read mapper, like [hisat2](http://ccb.jhu.edu/software/hisat2/index.shtml). Convert `.sam` file to a sorted `.bam`, using [samtools](https://github.com/samtools/samtools/releases)

    `~/samtools-1.8/samtools view -Su twilhelma_2014_vs_scaffolds_v1.sam | ~/samtools-1.8/samtools sort - -o twilhelma_2014_vs_scaffolds_v1.sorted.bam`

2. Extract the counts using the `sort_reads_from_bam.py` script.

    `~/samtools-1.8/samtools view twilhelma_2014_vs_scaffolds_v1.sorted.bam | ~/git/lavaLampPlot/sort_reads_from_bam.py -i - > twilhelma_2014_vs_scaffolds_v1.sorted.hits_from_bam.txt`

3. Compile the coverage, length, and GC content using the `hits_to_coverage.py` script.

    `~/git/lavaLampPlot/hits_to_coverage.py -f twilhelma_scaffolds_v1.fasta -b twilhelma_2014_vs_scaffolds_v1.sorted.hits_from_bam.txt > twilhelma_2014_vs_scaffolds_v1.coverage.tab`

4. Generate the plot with the R script `contig_gc_coverage.R`, which will automatically name the output file as `.pdf`.

    `Rscript ~/git/lavaLampPlot/contig_gc_coverage.R twilhelma_2014_vs_scaffolds_v1.coverage.tab`

### As a direct output of some assemblers
Some assemblers report the contig coverage directly within the fasta header, such as [SPAdes](https://github.com/ablab/spades). Another [script](https://bitbucket.org/wrf/sequences/src/master/spadescontigstocovgc.py) can convert this directly into a compatible table.

## ShinyApp test ##
This is the first attempt to convert the plotting into an interactive [ShinyApp](https://shiny.rstudio.com/tutorial/), with 3 slide bars to control what scaffolds are shown. Here, using the example data from [Ephydatia muelleri](https://spaces.facsci.ualberta.ca/ephybase/), the pink spot on the left shows a low-coverage bacterial partial-chromosome. Mouse-clicks displays the stats for that scaffold, and others nearby, below the plot.

![blobplot_screenshot_Emuelleri.png](https://github.com/wrf/lavaLampPlot/blob/master/Rshiny/blobplot_screenshot_Emuelleri.png)

## Usage considerations

![wrfs_tips_on_env_dna_bioinfo_v2.svg](https://github.com/wrf/lavaLampPlot/blob/master/sample_data/wrfs_tips_on_env_dna_bioinfo_v2.svg)

#### Choosing k-mer length
Due to the connection between kmer length and coverage, there is necessarily a balance between longer kmers, which will resolve the y-axis better, and higher coverage, which will resolve the x-axis better. Kmers between 31 and 41 tend to perform fairly well. Because the Trinity steps cannot use a kmer larger than 32, setting the kmer to 31 is advisable (odd numbers are preferable to prevent double counting from palindromic sequences).

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
The previous version of Trinity mode in kmersorter had a problem with the SRA headers generated from the [SRA Toolkit](http://www.ncbi.nlm.nih.gov/books/NBK158900/) fastq-dump (due to problems in the called Trinity scripts). The stats file generation would work normally, but then the stats cannot be sorted correctly and ultimately the process will fail merging the stats. Brian Haas suggested [this alternate command](https://github.com/trinityrnaseq/trinityrnaseq/wiki/Trinity-FAQ#ques_sra_fq_conversion) to generate the header (which one would also use during Trinity normalization).

`fastq-dump --split-files --defline-seq '@$sn[_$rn]/$ri' SRR1032106.sra`

This may no longer be a problem since `sort` is no longer called in the Trinity mode for kmersorter. This program was removed since the purpose was to organize the reads in a pair so they could be merged. The merging caused problems for extraction because the merge takes the average of the two coverage medians in a pair, so read coverage of 1 and 99 in a pair would make 50. This does not generate a count on the lava lamp plot, since the stats are counted separately in `fastqdumps2histo.py`. However, during extraction of regions in kmersorter, additional/junk reads outside of the defined extraction boundary were collected when neither read was in the defined region.

## Troubleshooting
Two problems have come up a few times, and sample plots are shown in the [error_plots folder](https://github.com/wrf/lavaLampPlot/tree/master/error_plots). For the case of [Acropora](https://github.com/wrf/lavaLampPlot/blob/master/error_plots/acropora_digitifera_15kb-MP_DRR001432_18Gb_k31_reads_u1000.pdf), there is a strip at the bottom. This was caused by the sequence in the read spanning multiple lines per read, instead of just one. For this situation, convert the reads to two-line fasta with the `fasta2twoline.py` script. For the other case of [Aiptasia](https://github.com/wrf/lavaLampPlot/blob/master/error_plots/aiptasia_pallida_550bp-PE_SRR1648349_11Gb_k31_reads_u1000.pdf), the plot is distorted since the reads were trimmed and are not the same length. To solve this, use the `-p` option in `fastqdumps2histo.py` to calculate the length and GC for each sequence.

