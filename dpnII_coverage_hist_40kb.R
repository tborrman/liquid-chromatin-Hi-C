#!/usr/bin/env Rscript

library(argparse)

parser <- ArgumentParser(description="create histogram of DpnII read coverage ")
parser$add_argument("-i", help="output file of GMM_copy_number.py (ex. *_GMM.txt)", type="character",
	required=TRUE)
args <- parser$parse_args()


df <- read.table(args$i, sep="\t", header=TRUE)
print(head(df))
max_x <- mean(df$reads) + (sd(df$reads)*2.5)

xvals <- seq(0,max(df$reads), 5)
# In order of ploidy
weights <- c(0.34, 0.40, 0.26)


png(gsub(".{4}$", "_hist.png", args$i), width=2500, height=2000, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(df$reads, breaks=seq(0, max(df$reads) + 2, 2), lty="blank", col="gray60",
 xlim=c(0,max_x), xlab= "Number of reads overlapping 40kb bin",
 main = "", freq=FALSE, ylim=c(0,0.008))

# 1st Gaussian
mu = 172.3
sigma = sqrt(2306.97)
abline(v= mu, col='blue')
abline(v= mu + sigma, col='blue', lty=2)
abline(v= mu - sigma, col='blue', lty=2)
lines(xvals, dnorm(xvals, mean= mu, sd= sigma)*weights[1], col="blue", lw=1.5)

# 2nd Gaussian
mu = 219.37
sigma = sqrt(2620.27)
abline(v= mu, col='red')
abline(v= mu + sigma, col='red', lty=2)
abline(v= mu - sigma, col='red', lty=2)
lines(xvals, dnorm(xvals, mean= mu, sd= sigma)*weights[2], col="red", lw=1.5)

# 3rd Gaussian
mu = 276.64
sigma = sqrt(3118.7)
abline(v= mu, col='green')
abline(v= mu + sigma, col='green', lty=2)
abline(v= mu - sigma, col='green', lty=2)
lines(xvals, dnorm(xvals, mean= mu, sd= sigma)*weights[3], col="green", lw=1.5)
dev.off()


# Per chromosome
# chroms <- c(paste(rep("chr", 21), 1:22, sep=""), "chrX", "chrY")
# png(paste("dpnII_coverage_hist_chrom.png", sep=""), width=4000, height=2000, res=300)
# par(mfrow=c(4,6), mar=c(2,1,1,1) + 0.5,oma=c(4,4,4,2) + 0.1 )
    
# for (chrom in chroms) {
#     chrom_df <- dpnII_df[which(dpnII_df$chrom == chrom),]
#     hist(chrom_df$reads, breaks=seq(0,max(chrom_df$reads) + 40, 40), lty="blank", col="black",xlim=c(0,4500),  xlab= "", ylab="", main=chrom)
    
# }
# dev.off()