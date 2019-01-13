#!/usr/bin/env Rscript
library(RColorBrewer)
library(argparse)

parser <- ArgumentParser(description="Make smooth scatter plot for input replicate matrices")
parser$add_argument("-x", help="Input matrix for Rep1 on x-axis", type="character", dest="x", required=TRUE)
parser$add_argument("-y", help="Input matrix for Rep2 on y-axis", type="character", dest="y", required=TRUE)
parser$add_argument("-o", help="outfile", type="character", dest="o", required=TRUE)
parser$add_argument("-t", help="title", type="character", dest="t", required=TRUE)
args <- parser$parse_args()

x_df <- read.table(gzfile(args$x), header=TRUE, row.names=1)
y_df <- read.table(gzfile(args$y), header=TRUE, row.names=1)

x <- as.vector(t(as.matrix(x_df)))
y <- as.vector(t(as.matrix(y_df)))

logx <- log2(x)
logy <- log2(y)

logx[!is.finite(logx)] <- NA
logy[!is.finite(logy)] <- NA

pdf(args$o, height=10, width=10)
par(mar=c(6, 5, 4, 2) + 0.1)
smoothScatter(logx, logy,
              xlab="Replicate 1 log2 interactions", ylab="Replicate 2 log2 interactions", 
              col='white', cex=1.25, main=args$t, cex.axis=1.5, cex.lab=1.5, cex.main=1.5,
              colramp=colorRampPalette(rev(brewer.pal(11,'Spectral'))))
text(4,15, 
	sprintf("r = %.3f", round(cor(x, y, method="pearson", use="complete.obs"),3)), 
		col="white", cex=2)
dev.off()
