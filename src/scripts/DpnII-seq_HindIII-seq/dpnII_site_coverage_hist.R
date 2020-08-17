#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description="create histogram of DpnII site coverage for 500kb bins")
args <- parser$parse_args()


dpnII_df <- read.table("dpnII_site_coverage_sorted.bed", sep="\t")
dpnII_df <- dpnII_df[1:4]
colnames(dpnII_df) <- c("chrom", "start", "end", "reads")
png("dpnII_site_coverage_hist.png", width=2500, height=2000, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(dpnII_df$reads, breaks=250, lty="blank", col="black", xlab= "# of DpnII sites overlapping 500kb bin", main="Histogram of DpnII site coverage per bin")
dev.off()