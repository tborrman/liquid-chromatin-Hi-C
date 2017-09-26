#!/usr/bin/env Rscript

library(argparse)
options(scipen=100)
parser <- ArgumentParser(description="input cis percent timepoint and return LOS")
parser$add_argument('-i', help='input cis percent timepoint file (ex. *.cistransratio.bedGraph', type="character",
                    required=TRUE)
args <- parser$parse_args()
wkdir <- getwd()

M_cis_df <- read.table(paste(wkdir, "HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45.cistransratio.bedGraph", sep="/"), 
                       sep="\t", skip=1)
D_cis_df <- read.table(paste(wkdir, args$i, sep="/"), 
                       sep="\t", skip=1)
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

write.table(out_df, paste(substr(args$i, 1, nchar(args$i) - 23), '.LOS.bedGraph', sep="")
            , row.names=FALSE, sep="\t", quote=FALSE)



