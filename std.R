#!/usr/bin/env Rscript

library(argparse)
options(scipen=100)
parser <- ArgumentParser(description="input std timepoint and return delta std")
parser$add_argument("-i", help="input cis percent timepoint file (ex. *_std.bedGraph)", type="character",
                    required=TRUE)
parser$add_argument("-m", help="Mock nuclei control std.bedGraph", type="character", required=TRUE)
args <- parser$parse_args()

M_std_df <- read.table(args$m, sep="\t", header=FALSE)
D_std_df <- read.table(args$i, sep="\t", header=FALSE)
# Fix bp indexing
M_std_df[2] <- M_std_df[2] - 1
D_std_df[2] <- D_std_df[2] - 1


colnames(M_std_df) <- c("chrom", "start", "end", "M_std")
colnames(D_std_df) <- c("chrom", "start", "end", "D_std")


# Merge M_std and D_std dfs
std_df <- merge(M_std_df, D_std_df, c("chrom","start"))
# Compute structure metric
delta_std <- round((std_df$M_std - std_df$D_std)/std_df$M_std,5)
std_df <- cbind(std_df, delta_std)

# Output
out_df <- std_df[c(1:3,7)]
colnames(out_df)[3] <- "end"

# Remove chrM
out_df <- out_df[out_df$chrom != "chrM",]

out_df$chrom <- factor(out_df$chrom, levels=c(paste("chr", 1:22, sep=""), "chrX"), ordered = TRUE)

# Order
out_df <- out_df[order(out_df$chrom, out_df$start),]

write.table(out_df, paste(substr(args$i, 1, nchar(args$i) - 13), '_delta_std.bedGraph', sep="")
            , row.names=FALSE, sep="\t", quote=FALSE)



