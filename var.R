#!/usr/bin/env Rscript

library(argparse)
options(scipen=100)
parser <- ArgumentParser(description="input var timepoint and return delta var")
parser$add_argument("-i", help="input cis percent timepoint file (ex. *_var.bedGraph)", type="character",
                    required=TRUE)
parser$add_argument("-m", help="Mock nuclei control var.bedGraph", type="character", required=TRUE)
args <- parser$parse_args()

M_var_df <- read.table(args$m, sep="\t", header=FALSE)
D_var_df <- read.table(args$i, sep="\t", header=FALSE)
# Fix bp indexing
M_var_df[2] <- M_var_df[2] - 1
D_var_df[2] <- D_var_df[2] - 1


colnames(M_var_df) <- c("chrom", "start", "end", "M_var")
colnames(D_var_df) <- c("chrom", "start", "end", "D_var")


# Merge M_var and D_var dfs
var_df <- merge(M_var_df, D_var_df, c("chrom","start"))
# Compute structure metric
delta_var <- round((var_df$M_var - var_df$D_var)/var_df$M_var,5)
var_df <- cbind(var_df, delta_var)

# Output
out_df <- var_df[c(1:3,7)]
colnames(out_df)[3] <- "end"

# Remove chrM
out_df <- out_df[out_df$chrom != "chrM",]

out_df$chrom <- factor(out_df$chrom, levels=c(paste("chr", 1:22, sep=""), "chrX"), ordered = TRUE)

# Order
out_df <- out_df[order(out_df$chrom, out_df$start),]

write.table(out_df, paste(substr(args$i, 1, nchar(args$i) - 13), '_delta_var.bedGraph', sep="")
            , row.names=FALSE, sep="\t", quote=FALSE)



