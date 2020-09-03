#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="Barplot of distances from fragment start to HindIII cut sites")
parser$add_argument("-f", help="distance file (ex. hindIII_abs_frag_distances)", type="character", dest="f", required=TRUE)

args <- parser$parse_args()
wkdir <- getwd()

distances <- scan(paste(wkdir, args$f, sep="/"))

png(paste(wkdir, "/", args$f, "_bar.png", sep=""), width=3000, height=2000, res=300)
par(mar=c(5, 5, 4, 2) + 0.1)
barplot(table(distances[distances<=10]), main="HBCRACKHiC-K562-HN-TD-R1", xlab="bp distance to nearest HindIII cut site",
	ylab="# of mapped fragments", cex.lab=1.5, cex.main=1.5)
box(bty="l")
dev.off()