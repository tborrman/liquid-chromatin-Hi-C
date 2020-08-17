df500 <- read.table("HBDpSeqK562-DN16hR1_S11_L003_coverage_500kb.bed", sep="\t", header=FALSE)
df500 <- df500[1:4]
colnames(df500) <- c("chrom", "start", "end", "dpnIIseq")
eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph", sep="\t", header=FALSE)
colnames(eigen) <- c("chrom", "start", "end", "PC1")
PC1 <- eigen$PC1
PC1[is.na(PC1)] <- NA


png("HBDpSeqK562-DN16hR1_S11_L003_coverage_500kb_hist.png", width=3000, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(df500$dpnIIseq, breaks=seq(0, max(df500$dpnIIseq, na.rm=TRUE) + 100, 100), lty="blank", col="gray60",
     xlab= "Number of reads overlapping 500kb bin",
     main = "", freq=FALSE, cex.lab=1.5,
     xlim=c(0,50000), ylim=c(0,0.00008)
     )
dev.off()


df40 <- read.table("HBDpSeqK562-DN16hR1_S11_L003_coverage_40kb.bed", sep="\t", header=FALSE)
df40 <- df40[1:4]
colnames(df40) <- c("chrom", "start", "end", "dpnIIseq")
png("HBDpSeqK562-DN16hR1_S11_L003_coverage_40kb_hist.png", width=3000, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(df40$dpnIIseq, breaks=seq(0, max(df40$dpnIIseq, na.rm=TRUE) + 100, 10), lty="blank", col="gray60",
     xlab= "Number of reads overlapping 40kb bin",
     main = "", freq=FALSE, cex.lab=1.5,
     xlim=c(0,8000), ylim=c(0, 0.001)
)
dev.off()


get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}

df40e <- cbind(df40, PC1)
df40_clean <- na.omit(df40e)

dpnII_residuals <- get_residuals(df40_clean$dpnIIseq, df40_clean$PC1)


png("HBDpSeqK562-DN16hR1_S11_L003_coverage_40kb_residuals_from_PC1_hist.png", width=3000, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(dpnII_residuals, breaks=8000, lty="blank", col="gray60",
     xlab= "16h DpnII-seq residuals",
     main = "", freq=FALSE, cex.lab=1.5,
     xlim=c(-3000,3000)
)
dev.off()

png("HBDpSeqK562-DN16hR1_S11_L003_coverage_vs_PC1_40kb.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(df40_clean$PC1, df40_clean$dpnIIseq,
    ylim=c(0,6000), pch=20, col=rgb(1,0,0, alpha=0.01), 
    xlab= "PC1", ylab="16h DpnII-seq",
    main = "", cex.lab=1.5
)
dev.off()

# 500kb residuals
eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/Houda_Ctrl_DpnII_K562.500000_all.zScore.eigen1_noheader_nochrY_sorted_bpfix.bedGraph", sep="\t", header=FALSE)
colnames(eigen) <- c("chrom", "start", "end", "PC1")
PC1 <- eigen$PC1
PC1[is.na(PC1)] <- NA

df500e <- cbind(df500, PC1)
df500_clean <- na.omit(df500e)

dpnII_residuals <- get_residuals(df500_clean$dpnIIseq, df500_clean$PC1)

png("HBDpSeqK562-DN16hR1_S11_L003_coverage_500kb_residuals_from_PC1_hist.png", width=3000, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
hist(dpnII_residuals, breaks=2000, lty="blank", col="gray60",
     xlab= "16h DpnII-seq residuals",
     main = "", freq=FALSE, cex.lab=1.5,
     xlim=c(-30000,30000)
)
dev.off()

png("HBDpSeqK562-DN16hR1_S11_L003_coverage_vs_PC1_500kb.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(df500_clean$PC1, df500_clean$dpnIIseq,
     pch=20, col=rgb(1,0,0, alpha=0.1), 
     xlab= "PC1", ylab="16h DpnII-seq",
     main = "", cex.lab=1.5,
     ylim=c(0,50000)
)
dev.off()

