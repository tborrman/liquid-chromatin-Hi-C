#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description = 'bin and average output from build_interactome.R')
parser$add_argument("-i", help = "output file from build_interactome.R", type="character", required=TRUE)
parser$add_argument("-o", help="name for table output file", type="character", required=TRUE)
args <- parser$parse_args()

H3K9me3_df_full <- read.table(args$i, header=TRUE, sep="\t")

# drop positions
H3K9me3_df <- H3K9me3_df_full[,c(4,5,6)]
#hist(H3K9me3_df$H3K9me3_I, breaks= 200)

# bin the data
H3K9me3_df$H3K9me3_I <- round(H3K9me3_df$H3K9me3_I, -3)
H3K9me3_df$H3K9me3_J <- round(H3K9me3_df$H3K9me3_J, -3)

# average interaction score for data at the same 
ave_H3K9me3_df <- aggregate(H3K9me3_df, by=list(H3K9me3_df$H3K9me3_I, H3K9me3_df$H3K9me3_J), FUN=mean)

write.table(ave_H3K9me3_df, args$o, sep="\t", quote=FALSE, row.names=FALSE) 
