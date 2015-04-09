# plot jellyfish kmer coverage vs gc content in lava lamp form
# v2 last modified 2015-04-02

# this script was generated using R version 3.1.3
# generate csv matrix of outputs using fastqdumps2histo.py

# run in the command line with:
# $ Rscript jellyfish_gc_coverage_blob_plot_v2.R

### CHANGE THIS FILE TO THE CSV FROM THE PREVIOUS STEP
f1 = "~/genomes/hydra_vulgaris/hydra_vulgaris_SRR1032106.k31.gc_vs_cov.csv"

### RENAME THIS TO THE DESIRED OUTPUT PDF
# change this for the output pdf, note that the pdf extension may be required
outfilename = "~/genomes/hydra_vulgaris/hydra_vulgaris_k31_vs_gc_with_line.pdf"

### CHANGE THIS VALUE BASED ON JELLYFISH -U OPTION
# or to subset or make the matrix smaller, change ymax, say to view from 0 to 500
ymax = 1000
# change zmax to correct color resolution, in case the real max is lower than zmax
zmax = 1000000

### THIS BLOCK BELOW SHOULD REMAIN MOSTLY UNCHANGED
# read in data and create log of matrix for color coding
# ymax here is also used to adjust the heatmap size
d1 = read.table(f1, header=FALSE, sep=",", colClasses='integer')[,1:ymax]

# at least one value is needed to make the matrix, in this case columns
# this alternatively could be nrow of kmer
m1 = matrix(unlist(d1), ncol=ymax, byrow=FALSE)
xmax = dim(m1)[1]

# instead calculate kmer as being number of rows - 1
kmer = xmax - 1

# this is needed to correct formatting of the GC percentage
roundedxdim = 5 * round(xmax/5)
gcpositions = pretty(c(0,roundedxdim))

# since GC is just a count, not a percentage, convert it to rounded percentage
gclabels = round(gcpositions / kmer * 100)

# must calculate sums before normalizing
d2 = colSums(m1)

# chops all values greater than zmax to zmax, which is 1 million by default
m1[m1>zmax] = zmax
# log all values to get better color resolution
m2 = log(m1, base=2)

# use this for log base 2
# first function generates color scale of 1000 colors from white to purple, then purple to red
colorcount = 1000
# 10 percent of the colors are for white to pink, the rest through the rainbow
white2pink = colorcount*0.1
pink2red = colorcount*0.9
colorset = c(colorRampPalette(c("#FFFFFF","#FF00FF"),alpha=0.9)(white2pink), rainbow(pink2red, s=0.95, e=0.93, v=0.75, alpha=0.9) )

# for coverage heatmap
### MUST UNCOMMENT FOR PDF OUTPUT
pdf(file=outfilename, width=10, heigh=8)

par(mar=c(4,4,4,4))
# this becomes the main chart label
mainlab = paste(kmer,"mer coverage vs GC",sep="")
# image() actually generates the heatmap
image(x=1:ncol(m2), y=1:nrow(m2), z=t(m2), col=colorset, xlab="Coverage", ylab="GC%", axes=FALSE, main=mainlab )
axis(1, at=pretty(c(0,ymax)), labels=pretty(c(0,ymax)) )
# because 0 is a position, all other positions have to be corrected
axis(2, at=gcpositions+1, labels=gclabels )

# for text overlay
textcol="white"
textsize = 1.3
# add one text command for each comment, such as identity of each blob
# x axis positions center at the given value, so must be increased for left alignment
# y axis positions correspond to GC counts, not percentage
### UNCOMMENT TO DISPLAY TEXT OVERLAY OR MODIFY
#text(160,9,"Hydra genomic", col=textcol, cex=textsize)
#text(125,24,"Curvibacter", col=textcol, cex=textsize)

# for line overlay
### COMMENT OUT LINES FROM PAR THROUGH PLOT TO REMOVE LINE OVERLAY
par(new=TRUE, mar=c(4,2,2,2))
# parameters for line generation
xmax2 = length(d2)
ymax2 = max(d2)
yscale = ceiling(log(ymax2,base=10))
plot(d2, type='l', xlab="Coverage", ylab="Unique kmers", xlim=c(0,xmax2), ylim=c(10^2,10^yscale), log='y', frame.plot=FALSE, axes=FALSE)

### MUST UNCOMMENT FOR PDF OUTPUT
dev.off()
print("PDF created!")