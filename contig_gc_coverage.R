# contig gc coverage  last modified 2018-05-14
# output should be tab-delimited as as:
# contigname  length  coverage  GC

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/tethya_wilhelma-genome/twilhelma_v2_dovetail/twilhelma_2014_dna_vs_dovetail_coverage.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

coveragedata = read.table(inputfile, header=FALSE, sep='\t')

names = coveragedata[,1]

covmax = 1000

contiglengths = coveragedata[,2]
magnituderange = range(log(contiglengths,base=10))
# size of each point
pchsize = log(contiglengths,base=10)-magnituderange[1]
# sizes of three reference points in the legend
legendsizes = c(0.25,0.5,0.75) * (magnituderange[2]-magnituderange[1]) + (magnituderange[1])
legendlabels = round(10^legendsizes)
legendpch = legendsizes - magnituderange[1]

gc = coveragedata[,4]
highgc = gc > 50
sum(contiglengths[highgc])
sum(contiglengths[!highgc])

# default is green
pointcolor = rep("#39bc6744", length(names))
longcontigs = contiglengths > legendlabels[2]
massivecontigs = contiglengths > legendlabels[3]
# midsize is blue
pointcolor[longcontigs] = "#386edc66"
# longest contigs are magenta
pointcolor[massivecontigs] = "#d51ea477"

# generate figure
pdf(file=outputfile, width=8, height=7)
par(mar=c(4.5,4.5,3,1))
plot(coveragedata[,3], coveragedata[,4], type='p', xlim=c(0,covmax), ylim=c(20,80), xlab="Mean coverage of mapped reads", ylab="GC%", pch=16, frame.plot=FALSE, col=pointcolor, cex.axis=1.5, cex=pchsize, main=inputfile, cex.lab=1.4)
legend(750,80, legend=legendlabels, pch=16, col=c("#39bc6799","#386edc99","#d51ea499"), pt.cex=legendpch, cex=1.1, title="Contig size (bp)")

text(1000,20,paste(length(names),"total contigs"), cex=1.2, pos=2)
dev.off()




#