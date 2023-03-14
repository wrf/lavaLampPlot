# contig gc coverage
# last modified 2022-12-13

# output should be tab-delimited as as:
# contigname  contignumber  length  coverage  GC  gaps
# this is the output from hits_to_coverage.py and spadescontigstocovgc.py

args = commandArgs(trailingOnly=TRUE)

inputfile = args[1]
#inputfile = "~/genomes/ephydatia_muelleri/ASM_HIC_394/Emuelleri_lib001_final_assembly.coverage_gc.tab"
outputfile = gsub("([\\w/]+)\\....$","\\1.pdf",gsub(".gz","",inputfile),perl=TRUE)

if (inputfile==outputfile) { stop("cannot parse input file to generate output file name, add a unique 3-letter suffix") }

print(paste0("Reading ",inputfile,", writing to ",outputfile))
print("expecting input as:  contigname  contignumber  length  coverage  GC  gaps")

coveragedata = read.table(inputfile, header=TRUE, sep='\t')
# reverse row order, since it usually expects biggest contigs first
# this ends up with largest contigs plotted last, meaning top layer
coveragedata = coveragedata[rev(1:nrow(coveragedata)),]
# does not work for app, as table is displayed in reverse order

names = coveragedata[,1]

if (length(args)>1) {
covmax = as.numeric(args[2])
} else {
covmax = 1000
}

gcmax = 80
gcmin = 20

contiglengths = coveragedata[["length"]]
magnituderange = range(log(contiglengths,base=10))
# size of each point
pchsize = log(contiglengths,base=10)-magnituderange[1]
# sizes of three reference points in the legend
legendsizes = c(0.25,0.5,0.75) * (magnituderange[2]-magnituderange[1]) + (magnituderange[1])
legendlabels = round(10^legendsizes)
legendpch = legendsizes - magnituderange[1]

totalsize = sum(contiglengths)

# print GC stats to user
gc = coveragedata[["GC"]]
highgc = gc > 50
print( paste(sum(highgc), "scaffolds of GC > 50%, for", sum(contiglengths[highgc]),"bases") )
print( paste(  sum(contiglengths[!highgc]),"bases of GC <= 50%") )
lowgc = gc < 10
print( paste(sum(lowgc), "scaffolds of GC < 10%, for", sum(contiglengths[lowgc]),"bases") )


# adjust min and max of GC scale, in case scaffolds go beyond
if (min(gc) < 20) {
	print( paste("lowest GC is", min(gc), ": setting minimum on graph") )
	gcmin = min(gc)
}

if (max(gc) > 80) {
	print( paste("highest GC is", max(gc), ": setting maximum on graph") )
	gcmax = max(gc)
}


# default is green
pointcolor = rep("#39bc6744", length(names))
longcontigs = contiglengths > legendlabels[2]
massivecontigs = contiglengths > legendlabels[3]
# midsize is blue
pointcolor[longcontigs] = "#386edc66"
# longest contigs are magenta
pointcolor[massivecontigs] = "#d51ea477"

# make GC range
gcspan = gcmax - gcmin

# generate figure
pdf(file=outputfile, width=8, height=7)
par(mar=c(4.5,4.5,3,1))
plot(coveragedata[["coverage"]], coveragedata[["GC"]], type='p', 
     xlim=c(0,covmax), ylim=c(gcmin,gcmax), frame.plot=FALSE, 
     xlab="Mean coverage of mapped reads", ylab="GC%", main=inputfile,
     pch=16, col=pointcolor, cex=pchsize, 
     cex.axis=1.5, cex.lab=1.4 )
legend(covmax,gcmax, legend=legendlabels, pch=16, pt.cex=legendpch,
       col=c("#39bc6799","#386edc99","#d51ea499"), cex=1.1, 
       title="Contig size (bp)", xjust=1)
text(covmax, gcmin + 0.05*gcspan, paste(round(totalsize/1000000,digits=1),"Mb"), cex=1.2, pos=2)
text(covmax, gcmin, paste(length(names),"total contigs"), cex=1.2, pos=2)
dev.off()




#