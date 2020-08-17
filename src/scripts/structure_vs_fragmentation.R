#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description="create scatterplot of CTR metric vs DpnII signal")
args <- parser$parse_args()

wkdir <- getwd()

dpnII_df <- read.table(paste(wkdir, "dpnII_coverage_sorted.bed", sep="/"), sep="\t")
M_ctr_df <- read.table(paste(wkdir, "HBCRACKHiC-K562-MN-R1__hg19__genome__C-500000-iced.cistransratio.bedGraph", sep="/"), sep="\t")
D_ctr_df <- read.table(paste(wkdir, "HBCRACKHiC-K562-DN-R1__hg19__genome__C-500000-iced.cistransratio.bedGraph", sep="/"), sep="\t")

# Fix bp indexing
M_ctr_df[2] <- M_ctr_df[2] - 1
D_ctr_df[2] <- D_ctr_df[2] - 1

dpnII_df <- dpnII_df[1:4]
colnames(M_ctr_df) <- c("chrom", "start", "end", "M_ctr")
colnames(D_ctr_df) <- c("chrom", "start", "end", "D_ctr")
colnames(dpnII_df) <- c("chrom", "start", "end", "dpnII")

# Merge M_ctr and D_ctr dfs
ctr_df <- merge(M_ctr_df, D_ctr_df, c("chrom","start"))
# Compute structure metric
structure <- (ctr_df$M_ctr - ctr_df$D_ctr)/ctr_df$M_ctr
ctr_df <- cbind(ctr_df, structure)

# Merge dpnII data ctr dfs
full_df <- merge(ctr_df, dpnII_df, c("chrom", "start"))

png("structure_vs_fragmentation_full.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(full_df$dpnII, full_df$structure, xlab="# of DpnII reads overlapping 500kb bin",
     ylab="[(mCTR - dCTR) / mCTR] for 500kb bin")
text(10000, -0.5, paste("r = ",round((cor(full_df$dpnII, full_df$structure)),2),sep=""))
dev.off()

limit_df <- full_df[full_df$structure > 0 & full_df$dpnII < 5000,]


png("structure_vs_fragmentation_limit.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(limit_df$dpnII, limit_df$structure, xlab="# of DpnII reads overlapping 500kb bin",
     ylab="[(mCTR - dCTR) / mCTR] for 500kb bin")
text(4000, 0.2, paste("r = ",round((cor(limit_df$dpnII, limit_df$structure)),2),sep=""))
dev.off()




