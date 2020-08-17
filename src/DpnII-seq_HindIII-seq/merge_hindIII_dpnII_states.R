#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description="merge hindIII and dpnII copy number states into single file")
parser$add_argument("-hn", help="hindIII GMM file (ex. hindIII_coverage_sorted_500kb_GMM.txt)", 
	type="character", required=TRUE)
parser$add_argument("-dn", help="dpnII GMM file (ex.dpnII_coverage_sorted_500kb_R1_GMM.txt)",
	type="character", required=TRUE)
args <- parser$parse_args()

hn <- read.table(args$hn, sep="\t", header=TRUE)
dn <- read.table(args$dn, sep="\t", header=TRUE)


hn_copy <- hn$copy
dn_copy <- dn$copy
copy <- hn_copy

not_eq <- hn_copy != dn_copy
nas <- is.na(hn_copy) | is.na(dn_copy)

copy[not_eq] <- NaN
copy[nas] <- NaN

df <- data.frame(hn[1:3], copy)

write.table(df, "K562_copynumber_500kb.txt", sep="\t", row.names=FALSE, quote=FALSE)
