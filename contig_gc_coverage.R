# contig gc coverage  last modified 2018-08-29
# output should be tab-delimited as as:
# contigname  contignumber  length  coverage  GC  gaps

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/tethya_wilhelma-genome/twilhelma_v2_dovetail/twilhelma_2014_dna_vs_dovetail_coverage.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",inputfile,perl=TRUE)

coveragedata = read.table(inputfile, header=TRUE, sep='\t')

names = coveragedata[,1]

if (length(args)>1) {
covmax = as.numeric(args[2])
} else {
covmax = 1000
}

gcmax = 80

contiglengths = coveragedata[["length"]]
magnituderange = range(log(contiglengths,base=10))
# size of each point
pchsize = log(contiglengths,base=10)-magnituderange[1]
# sizes of three reference points in the legend
legendsizes = c(0.25,0.5,0.75) * (magnituderange[2]-magnituderange[1]) + (magnituderange[1])
legendlabels = round(10^legendsizes)
legendpch = legendsizes - magnituderange[1]

totalsize = sum(contiglengths)

gc = coveragedata[["GC"]]
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
plot(coveragedata[["coverage"]], coveragedata[["GC"]], type='p', xlim=c(0,covmax), ylim=c(20,gcmax), xlab="Mean coverage of mapped reads", ylab="GC%", pch=16, frame.plot=FALSE, col=pointcolor, cex.axis=1.5, cex=pchsize, main=inputfile, cex.lab=1.4)
legend(covmax,gcmax, legend=legendlabels, pch=16, col=c("#39bc6799","#386edc99","#d51ea499"), pt.cex=legendpch, cex=1.1, title="Contig size (bp)", xjust=1)
text(covmax,23,paste(round(totalsize/1000000,digits=1),"Mb"), cex=1.2, pos=2)
text(covmax,20,paste(length(names),"total contigs"), cex=1.2, pos=2)
dev.off()




#