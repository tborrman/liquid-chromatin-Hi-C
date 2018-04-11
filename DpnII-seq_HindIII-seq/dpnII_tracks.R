#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description="Plot dpnII signal and estimated GMM means and standard deviations")
args <- parser$parse_args()

# DpnII signal in 500kb bins
dpnII_file <- "C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/maxins_10kb/dpnII_R1_dist_3/dpnII_3_distance_coverage_sorted_500kb_R1.bed"
dpnII_df <- read.table(dpnII_file, sep="\t")
dpnII_df <- dpnII_df[1:4]
colnames(dpnII_df) <- c("chrom", "start", "end", "reads")
# Copy number
copy_file <- "C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/maxins_10kb/dpnII_R1_dist_3/K562_copynumber_500kb.bed"
copy_df <- read.table(copy_file, sep="\t")
colnames(copy_df) <- c("chrom", "start", "end", "copy")

options(scipen=100)
chroms <- c(paste(rep("chr", 22), 1:22, sep=""),"chrX")
#Plot
for (i in seq(1,23)) {
  dpnII_chr <- dpnII_df[which(dpnII_df$chrom == chroms[i]),]
  copy_chr <- copy_df[which(copy_df$chrom == chroms[i]),]
  dpnII_gen <- dpnII_chr$start + ((dpnII_chr$end - dpnII_chr$start)/ 2)
  file_prefix <- "C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/maxins_10kb/dpnII_R1_dist_3/dpnII_tracks/dpnII_R1_dist_3_chr"
  pdf(paste(file_prefix, i, ".pdf",sep=""), width=18, height=5)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnII_gen/1000000, dpnII_chr$reads, type='l',col="black", xlim=c(0,max(dpnII_gen)/1000000), 
       ylim=c(0,1600), xlab= "Genomic position (Mb)", ylab="Raw DpnII-seq signal", cex.lab = 1.5,
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  copy <- copy_chr$copy
  copy_col <- rep(1, length(copy))
  copy_col[is.na(copy)] <- 'black'
  copy_col[copy == 2] <- 'blue'
  copy_col[copy == 3] <- 'red'
  copy_col[copy == 4] <- 'green'
  points(dpnII_gen/1000000, dpnII_chr$reads, pch=20, col=copy_col)
  dev.off()
}
# Copy number corrected
dpnIIc_file <- "C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/maxins_10kb/dpnII_R1_dist_3/dpnII_3_distance_copy_correct_500kb_R1.bed"
dpnIIc_df <- read.table(dpnIIc_file, sep="\t")
dpnIIc_df <- dpnIIc_df[1:4]
colnames(dpnIIc_df) <- c("chrom", "start", "end", "reads")
#Plot
for (i in seq(1,23)) {
  dpnIIc_chr <- dpnIIc_df[which(dpnIIc_df$chrom == chroms[i]),]
  dpnIIc_gen <- dpnIIc_chr$start + ((dpnIIc_chr$end - dpnIIc_chr$start)/ 2)
  file_prefix <- "C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/maxins_10kb/dpnII_R1_dist_3/dpnII_tracks/dpnII_R1_dist_3_copy_correct_chr"
  pdf(paste(file_prefix, i, ".pdf",sep=""), width=18, height=5)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnIIc_gen/1000000, dpnIIc_chr$reads, type='l',col="black", xlim=c(0,max(dpnIIc_gen)/1000000), 
       ylim=c(0,1600), xlab= "Genomic position (Mb)", ylab="Corrected DpnII-seq signal", cex.lab = 1.5,
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  points(dpnIIc_gen/1000000, dpnIIc_chr$reads, pch=20, col="blue")
  dev.off()
}
