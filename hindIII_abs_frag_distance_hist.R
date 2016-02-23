#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Histogram of distances from fragment start/end to nearest HindIII cut site")
parser$add_argument("-f", help="distance file (ex. hindIII_abs_frag_distances)", type="character", dest="f", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

distances <- scan(paste(wkdir, args$f, sep="/"))

png(paste(wkdir, "/", args$f, ".png", sep=""), width=4000, height=2000, res=300)
par(mar=c(5, 5, 4, 2) + 0.1)
hist(distances, breaks=seq(0,max(distances)), xlim=c(0, 100), main="HBCRACKHiC-K562-HN-TD-R1",
	ylab="# of mapped fragments", xlab="bp distance to nearest HindIII cut site", cex.lab=1.5, cex.main=1.5)
box(bty="l")
dev.off()