#!/usr/bin/env Rscript

library(argparse)
parser <- ArgumentParser(description="Plot dpnII signal and estimated GMM means and standard deviations")
parser$add_argument("-i", help="coverage bed file", type="character", required=TRUE)
args <- parser$parse_args()


# DpnII signal in 500kb bins
dpnII_df <- read.table(args$i, sep="\t", header=FALSE)
dpnII_df <- dpnII_df[1:4]
colnames(dpnII_df) <- c("chrom", "start", "end", "reads")

# Copy number
copy_file <- "C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/HiC_copynumber/COSMIC/COSMIC_K562_copy_number_40kb.bedGraph"
copy_df <- read.table(copy_file, sep="\t", header=FALSE)
colnames(copy_df) <- c("chrom", "start", "end", "copy")

options(scipen=100)
chroms <- c(paste(rep("chr", 22), 1:22, sep=""),"chrX")
#Plot
for (i in seq(1,23)) {
  dpnII_chr <- dpnII_df[which(dpnII_df$chrom == chroms[i]),]
  copy_chr <- copy_df[which(copy_df$chrom == chroms[i]),]
  dpnII_gen <- dpnII_chr$start + ((dpnII_chr$end - dpnII_chr$start)/ 2)
  args_prefix <- gsub(".*coverage/40kb/", "", args$i)
  file_prefix <- paste("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/plots/dpnII_tracks/", gsub(".bed", "_", args_prefix), sep="")
  png(paste(file_prefix, i, ".png",sep=""), width=5500, height=1500, res=300)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnII_gen/1000000, dpnII_chr$reads, type='l',col="black", xlim=c(0,max(dpnII_gen)/1000000), 
       xlab= "Genomic position (Mb)", ylab="Raw DpnII-seq signal", cex.lab = 1.5, ylim=c(0,4500),
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  copy <- copy_chr$copy
  copy_col <- rep(1, length(copy))
  copy_col[is.na(copy)] <- 'black'
  copy_col[copy > 4 | copy < 2] <- 'black'
  copy_col[copy == 2] <- 'blue'
  copy_col[copy == 3] <- 'red'
  copy_col[copy == 4] <- 'green'
  points(dpnII_gen/1000000, dpnII_chr$reads, pch=20, col=copy_col)
  dev.off()
  pdf(paste(file_prefix, i, ".pdf",sep=""), width=12, height=3.5)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnII_gen/1000000, dpnII_chr$reads, type='l',col="black", xlim=c(0,max(dpnII_gen)/1000000), 
       xlab= "Genomic position (Mb)", ylab="Raw DpnII-seq signal", cex.lab = 1.5, ylim=c(0,4500),
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  points(dpnII_gen/1000000, dpnII_chr$reads, pch=20, col=copy_col, cex=0.5)
  dev.off()
}

# Copy number corrected
copy_arg <- gsub("coverage", "copy_correct_coverage", args$i)
dpnIIc_df <- read.table(copy_arg, sep="\t")
dpnIIc_df <- dpnIIc_df[1:4]
colnames(dpnIIc_df) <- c("chrom", "start", "end", "reads")
#Plot
for (i in seq(1,23)) {
  dpnIIc_chr <- dpnIIc_df[which(dpnIIc_df$chrom == chroms[i]),]
  dpnIIc_gen <- dpnIIc_chr$start + ((dpnIIc_chr$end - dpnIIc_chr$start)/ 2)
  args_prefix <- gsub(".*coverage/40kb/", "", copy_arg)
  file_prefix <- paste("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/plots/dpnII_tracks/", gsub(".bed", "_", args_prefix), sep="")
  png(paste(file_prefix, i, ".png",sep=""), width=5500, height=1500, res=300)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnIIc_gen/1000000, dpnIIc_chr$reads, type='l',col="black", xlim=c(0,max(dpnIIc_gen)/1000000), 
       xlab= "Genomic position (Mb)", ylab="Corrected DpnII-seq signal", cex.lab = 1.5, ylim=c(0,4500),
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  points(dpnIIc_gen/1000000, dpnIIc_chr$reads, pch=20, col="blue")
  dev.off()
  pdf(paste(file_prefix, i, ".pdf",sep=""), width=12, height=3.5)
  par(mar=c(5,6,4,2) + 0.1)
  plot(dpnIIc_gen/1000000, dpnIIc_chr$reads, type='l',col="black", xlim=c(0,max(dpnIIc_gen)/1000000), 
       xlab= "Genomic position (Mb)", ylab="Corrected DpnII-seq signal", cex.lab = 1.5, ylim=c(0,4500),
       main = paste("Chromosome ", i, sep=""), cex.main=1.5)
  points(dpnIIc_gen/1000000, dpnIIc_chr$reads, pch=20, col="blue", cex=0.5)
  dev.off()
  
}

