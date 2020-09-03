#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Histogram of distances from read start to DpnII sites")
parser$add_argument("-f", help="distance file (ex. dpnII_distances)", type="character", dest="f", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

distances <- scan(paste(wkdir, args$f, sep="/"))

png(paste(wkdir, "/", args$f, ".png", sep=""), width=4000, height=2000, res=300)
par(mar=c(5, 5, 4, 2) + 0.1)
hist(distances, breaks=2000, xlim=c(-1250, 1250), main="HBCRACKHiC-K562-DN-TD-R1",
	ylab="# of mapped reads", xlab="distance to DpnII site (GATC)", cex.lab=1.5, cex.main=1.5)
abline(v=0, col="blue")
abline(v=96, col="red")
axis(1, at=0, col.axis="blue")
axis(1, at=96, col.axis="red")
box(bty="l")
dev.off()