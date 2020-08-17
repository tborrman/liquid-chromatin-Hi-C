eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))
dseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))
min5 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-5m-DPnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)
hour1 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-1h-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)
hour2 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-2h-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)
hour3 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-3h-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)
hour4 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-4h-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)
ON <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-ON-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS.bedGraph", sep="\t", header=TRUE)

colnames(min5) <- c("chrom", "start", "end", "LOS5min")
colnames(hour1) <- c("chrom", "start", "end", "LOS1hour")
colnames(hour2) <- c("chrom", "start", "end", "LOS2hour")
colnames(hour3) <- c("chrom", "start", "end", "LOS3hour")
colnames(hour4) <- c("chrom", "start", "end", "LOS4hour")
colnames(ON) <- c("chrom", "start", "end", "LOS16hour")

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
thalf_chr2 <- thalf[thalf$chrom == "chr2",]
dseq_chr2 <- dseq[dseq$chrom == "chr2",]
min5_chr2 <- min5[min5$chrom == "chr2",]
hour1_chr2 <- hour1[hour1$chrom == "chr2",]
hour2_chr2 <- hour2[hour2$chrom == "chr2",]
hour3_chr2 <- hour3[hour3$chrom == "chr2",]
hour4_chr2 <- hour4[hour4$chrom == "chr2",]
ON_chr2 <- ON[ON$chrom == "chr2",]

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}


d <- data.frame(eigen_chr2[1:4], thalf_chr2["hl"], dseq_chr2["DpnIIseq"],
                min5_chr2["LOS5min"], hour1_chr2["LOS1hour"], hour2_chr2["LOS2hour"],
                hour3_chr2["LOS3hour"], hour4_chr2["LOS4hour"], ON_chr2["LOS16hour"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$hl[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$DpnIIseq[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS5min[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS1hour[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS2hour[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS3hour[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS4hour[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS16hour[which(is.na(d), arr.ind = TRUE)[,1]] <- NA

d_clean <- na.omit(d)
hl_residuals <- get_residuals(d_clean$hl, d_clean$DpnIIseq)
resid_df <- cbind(d_clean[1:3], hl_residuals)


d <- merge(d, resid_df, by=c("chrom", "start", "end"), all = TRUE)
d <- d[order(d$start),]


pdf("fig4tracks.pdf", height = 6, width = 6)
par(mfrow=c(4,1), mar=c(2, 4, 0, 2) + 0.1)
# PC1
plot(d$start/1000000, d$PC1, type="n", xlim=c(120,240),
      xlab = "", ylab= "PC1", axes=FALSE)
A <- d$PC1
A[is.na(A)] <- 0
A[A<0] <- 0
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(A, rep(0, length(A))), col="red",
        border=NA)
B <- d$PC1
B[is.na(B)] <- 0
B[B>0] <- 0
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(B, rep(0, length(B))), col="blue",
        border=NA)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
# LOS
plot(d$start/1000000, d$LOS5min, type="l", col="orange", xlim=c(120,240),
     xlab = "", ylab= "LOS", axes=FALSE, lwd=0.5, ylim=c(-0.4, 1))
lines(d$start/1000000, d$LOS1hour, type="l", col="springgreen4",
      xlab = "", lwd=0.5)
lines(d$start/1000000, d$LOS2hour, type="l", col="pink1",
      xlab = "", lwd=0.5)
lines(d$start/1000000, d$LOS3hour, type="l", col="purple",
      xlab = "", lwd=0.5)
lines(d$start/1000000, d$LOS4hour, type="l", col="turquoise",
      xlab = "", lwd=0.5)
lines(d$start/1000000, d$LOS16hour, type="l", col="red2",
      xlab = "", lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
# Half-life
plot(d$start/1000000, d$hl, type="l", col="dodgerblue", xlim=c(120,240),
     ylim=c(40,110), xlab = "", ylab= bquote("t"[1/2]), axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
# Half-life residuals
plot(d$start/1000000, d$hl_residuals, type="l", col="black", xlim=c(120,240),
     ylim=c(-30,30), xlab = "", ylab= bquote("t"[1/2]~" residuals"), axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()


