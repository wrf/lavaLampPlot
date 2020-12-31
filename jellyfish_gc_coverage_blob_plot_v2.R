# plot jellyfish kmer coverage vs gc content in lava lamp form
# v2 last modified 2020-12-29

# this script was generated using R version 3.1.3
# generate csv matrix of outputs using fastqdumps2histo.py

# run in the command line with:
# $ Rscript jellyfish_gc_coverage_blob_plot_v2.R inputfile.csv output_graph.pdf [y-max]

### FILE NAMES FOR INPUT AND OUTPUT ARE TAKEN FROM COMMAND LINE ARGUMENTS
args = commandArgs(trailingOnly=TRUE)
f1 = args[1]
outfilename = args[2]

### SET THESE PARAMETERS BASED ON DESIRED OUTPUT
# withscale draws a scalebar to the right, withline means draw the log line overlay
withscale  = TRUE
withline   = FALSE

refinecolors = TRUE
# if refinecolors is true, it will assume that the value set as zmax below (default to 1M) should be purple
# and the max in the graph should adjust to that; if that max is below zmax, the highest color changes
# this means that colors would be comparable across different samples, rather than just the scale values

# if refinecolors is false, purple will always be the highest value, regardless of the set zmax

### CHANGE THIS VALUE BASED ON JELLYFISH -U OPTION
# or to subset or make the matrix smaller, change ymax, say to view from 0 to 500
if (length(args)>2) {
ymax = as.numeric(args[3])
} else {
ymax = 1000
}
print(ymax)
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
# 16 percent of the colors are for white to pink, the rest through the rainbow
white2pink = colorcount*0.166
pink2red = colorcount*0.834

# example monochromatic colorset
#colorset = c(colorRampPalette(c("#DDDDFF","#5893db"))(white2pink), colorRampPalette(c("#5893db","#000000"))(pink2red) )
# standard rainbow colorset
colorset = c(colorRampPalette(c("#FFDDDD","#CC006E"),alpha=0.9)(white2pink), rainbow(pink2red, start=0.91, end=0.89, v=0.8, alpha=0.9) )

# zmax should therefore be either the max of m1 or the defined zmax, whichever is lower
truezmax = min(max(m1),zmax)
if (refinecolors==TRUE){
# colorset is then bound by the truezmax, so the highest color value will be constant across 
	colorset = colorset[1:as.integer(log(truezmax, base=2)/log(zmax, base=2)*1000)]
}
# must recalculate zmax if zmax is much greater than highest value in m1, otherwise colors are wrong
zmax = truezmax

# for coverage heatmap
### MUST UNCOMMENT FOR PDF OUTPUT
pdf(file=outfilename, width=10, heigh=8)

# this becomes the main chart label
mainlab = paste0(kmer,"mer coverage vs GC in ",f1)

if(withscale==TRUE){
	# this generates a 2 column figure with a narrow scalebar on the right side
	layout(matrix(c(1,2),1,2), widths=c(9,1))
	par(mar=c(4,4,4,2))	
}else{
	# this explicit sets up the single panel figure
	layout(matrix(c(1,1),1,1), widths=c(1))
	par(mar=c(4,4,4,4))
}
# image() actually generates the heatmap
image(x=1:ncol(m2), y=1:nrow(m2), z=t(m2), col=colorset, xlab="", ylab="", axes=FALSE, main=mainlab )
mtext("Coverage", side=1, line=2.5, cex=1.4)
mtext("GC%", side=2, line=2.5, cex=1.4)
axis(1, at=pretty(c(0,ymax)), labels=pretty(c(0,ymax)), cex.axis=1.4 )
# because 0 is a position, all other positions have to be corrected
axis(2, at=gcpositions+1, labels=gclabels, cex.axis=1.4 )

# for text overlay
textcol="white"
textsize = 1.3
# add one text command for each comment, such as identity of each blob
# x axis positions center at the given value, so must be increased for left alignment
# y axis positions correspond to GC counts, not percentage
# so they can be converted to values with percent GC x kmer
### UNCOMMENT TO DISPLAY TEXT OVERLAY OR MODIFY
#text(160,0.29*kmer,"Hydra genomic", col=textcol, cex=textsize)
#text(125,0.77*kmer,"Curvibacter", col=textcol, cex=textsize)

# for line overlay
if (withline==TRUE){
	par(new=TRUE, mar=c(4,2,2,2))
	# parameters for line generation
	xmax2 = length(d2)
	ymax2 = max(d2)
	yscale = ceiling(log(ymax2,base=10))
	plot(d2, type='l', xlab="Coverage", ylab="Unique kmers", xlim=c(0,xmax2), ylim=c(10^2,10^yscale), log='y', frame.plot=FALSE, axes=FALSE)
}
# for log scalebar
if (withscale==TRUE){
	par(new=FALSE, mar=c(4,2,2,1))
	l2zmax = log(zmax, base=2)
	notch = log(10^seq(0,10), base=2)
	mnotch = notch[notch<=l2zmax]
	mnotchpos = mnotch/l2zmax
	colorbar = as.matrix(seq(0,l2zmax,ceiling(l2zmax)/colorcount),1,colorcount+1)
	image(z=t(colorbar),col=colorset, axes=FALSE, xlab="Counts")
	axis(2,at=mnotchpos, labels=c(round(2^mnotch)), las=1)
}

### MUST UNCOMMENT FOR PDF OUTPUT
dev.off()
print("PDF created!")