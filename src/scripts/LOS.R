#!/usr/bin/env Rscript

library(argparse)
options(scipen=100)
parser <- ArgumentParser(description="input cis percent timepoint and return LOS")
parser$add_argument("-i", help="input cis percent timepoint file (ex. *_cispercent.bedGraph)", type="character",
                    required=TRUE)
parser$add_argument("-m", help="Mock nuclei control cispercent.bedGraph", type="character", required=TRUE)
args <- parser$parse_args()

M_cis_df <- read.table(args$m, sep="\t", header=FALSE)
D_cis_df <- read.table(args$i, sep="\t", header=FALSE)
# Fix bp indexing
M_cis_df[2] <- M_cis_df[2] - 1
D_cis_df[2] <- D_cis_df[2] - 1


colnames(M_cis_df) <- c("chrom", "start", "end", "M_cis")
colnames(D_cis_df) <- c("chrom", "start", "end", "D_cis")


# Merge M_cis and D_cis dfs
cis_df <- merge(M_cis_df, D_cis_df, c("chrom","start"))
# Compute structure metric
LOS <- round((cis_df$M_cis - cis_df$D_cis)/cis_df$M_cis,5)
cis_df <- cbind(cis_df, LOS)

# Output
out_df <- cis_df[c(1:3,7)]
colnames(out_df)[3] <- "end"

# Remove chrM
out_df <- out_df[out_df$chrom != "chrM",]

out_df$chrom <- factor(out_df$chrom, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"), ordered = TRUE)

# Order
out_df <- out_df[order(out_df$chrom, out_df$start),]

write.table(out_df, paste(substr(args$i, 1, nchar(args$i) - 20), '_LOS.bedGraph', sep="")
            , row.names=FALSE, sep="\t", quote=FALSE)



