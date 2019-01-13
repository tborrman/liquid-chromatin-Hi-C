#!/bin/env/Rscript
library(RColorBrewer)
library(argparse)

parser <- ArgumentParser(description="Make smooth scatter plot for input matrices")
parser$add_argument("-x", help="Input matrix for x-axis", type="character", dest="x", required=TRUE)
parser$add_argument("-y", help="Input matrix for y-axis", type="character", dest="y", required=TRUE)
parser$add_argument("-o", help="outfile", type="character", dest="o", required=TRUE)
args <- parser$parse_args()

x_df <- read.table(gzfile(args$x), header=TRUE, row.names=1)
y_df <- read.table(gzfile(args$y), header=TRUE, row.names=1)

x <- as.vector(t(as.matrix(x_df)))
y <- as.vector(t(as.matrix(y_df)))

png(args$o, height=2500, width=2500, res=300)
smoothScatter(x, y,
              xlab=args$x, ylab=args$y, col='white', cex=1.25,
              colramp=colorRampPalette(rev(brewer.pal(11,'Spectral'))))
dev.off()

png(paste("log_", args$o, sep=""), height=2500, width=2500, res=300)
smoothScatter(log2(x), log2(y),
              xlab=args$x, ylab=args$y, col='white', cex=1.25,
              colramp=colorRampPalette(rev(brewer.pal(11,'Spectral'))))
dev.off()
