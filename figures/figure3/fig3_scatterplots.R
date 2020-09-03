eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
LOS_DpnII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse2/LOS6Mb/HBHiCK562DN10-3hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_d"))
LOS_HindIII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse2/LOS6Mb/HBHiCK562HN50-4hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS_removed_outliers.bedGraph",
                          sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_h"))
DpnIIseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN4hR1_S8_L002_copy_correct_coverage_40kb.bed",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "signal"))

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
LOS_DpnII_chr2 <- LOS_DpnII[LOS_DpnII$chrom == "chr2",]
LOS_HindIII_chr2 <- LOS_HindIII[LOS_HindIII$chrom == "chr2",]
DpnIIseq_chr2 <- DpnIIseq[DpnIIseq$chrom == "chr2",]

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}

d <- data.frame(eigen_chr2[1:4], LOS_DpnII_chr2["LOS_d"],
                LOS_HindIII_chr2["LOS_h"], DpnIIseq_chr2["signal"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_d[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_h[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA

pdf("scatterplot.pdf", height=9, width=5)
par(mfrow=c(3,2), mar=c(5, 4, 0, 2) + 0.1)
d_clean <- na.omit(d)

LOShPC1m <- lm(d_clean$LOS_h~d_clean$PC1)
plot(d$PC1, d$LOS_h, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS HindIII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_h, method="spearman", use="complete.obs"), 3))
text(0.01, -0.10, bquote(rho == .(rho) ), cex=2)
abline(LOShPC1m$coefficients[1], LOShPC1m$coefficients[2], col="grey60", lwd = 1.5)

LOSdPC1m <- lm(d_clean$LOS_d~d_clean$PC1)
plot(d$PC1, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.78, bquote(rho == .(rho) ), cex=2)
abline(LOSdPC1m$coefficients[1], LOSdPC1m$coefficients[2], col="grey60", lwd=1.5)

LOSdDpnIIseqm <- lm(d_clean$LOS_d~d_clean$signal)
plot(d$signal, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$signal, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(1700, 0.78, bquote(rho == .(rho) ), cex=2)
abline(LOSdDpnIIseqm$coefficients[1], LOSdDpnIIseqm$coefficients[2], col="grey60", lwd = 1.5)

d_subset <- d[d$signal > 1000 & d$signal < 1100,]
d_subset_clean <- na.omit(d_subset)
subset_m <- lm(d_subset_clean$LOS_d~d_subset_clean$PC1)
plot(d_subset$PC1, d_subset$LOS_d, col=ifelse(d_subset$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=1, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d_subset$PC1, d_subset$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.81, bquote(rho == .(rho) ), cex=2)
abline(subset_m$coefficients[1], subset_m$coefficients[2], col="grey60", lwd = 1.5)

LOS_residuals <- get_residuals(d_clean$LOS_d, d_clean$signal)
PC1_residuals <- get_residuals(d_clean$PC1, d_clean$signal)
resid_LOS_PC1_m <-lm(LOS_residuals ~ PC1_residuals)
plot(PC1_residuals, LOS_residuals, col=ifelse(d_clean$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab="LOS DpnII (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(PC1_residuals, LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0, -0.075, bquote(rho == .(rho) ), cex=2)
abline(resid_LOS_PC1_m$coefficients[1], resid_LOS_PC1_m$coefficients[2], col="grey60", lwd=1.5)

DpnIIseq_residuals <- get_residuals(d_clean$signal, d_clean$PC1)
LOS_residuals <- get_residuals(d_clean$LOS_d, d_clean$PC1)
resid_LOS_DpnIIseq_m <- lm(LOS_residuals ~ DpnIIseq_residuals)
plot(DpnIIseq_residuals, LOS_residuals, col=ifelse(d_clean$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq residuals", ylab="LOS DpnII (%) residuals", cex=0.5,
     xlim=c(-1000,1000), cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(DpnIIseq_residuals, LOS_residuals, method="spearman", use="complete.obs"), 3))
text(400, -0.075, bquote(rho == .(rho) ), cex=2)
abline(resid_LOS_DpnIIseq_m$coefficients[1], resid_LOS_DpnIIseq_m$coefficients[2], col="grey60", lwd=1.5)
dev.off()


png("scatterplot.png", height=3500, width=1600, res=300)
par(mfrow=c(3,2), mar=c(5, 4, 2, 2) + 0.1)
d_clean <- na.omit(d)

LOShPC1m <- lm(d_clean$LOS_h~d_clean$PC1)
plot(d$PC1, d$LOS_h, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_h, method="spearman", use="complete.obs"), 3))
# text(0.01, -0.10, bquote(rho == .(rho) ), cex=2)
abline(LOShPC1m$coefficients[1], LOShPC1m$coefficients[2], col="grey60", lwd = 3)

LOSdPC1m <- lm(d_clean$LOS_d~d_clean$PC1)
plot(d$PC1, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_d, method="spearman", use="complete.obs"), 3))
# text(0.01, 0.78, bquote(rho == .(rho) ), cex=2)
abline(LOSdPC1m$coefficients[1], LOSdPC1m$coefficients[2], col="grey60", lwd= 3)

LOSdDpnIIseqm <- lm(d_clean$LOS_d~d_clean$signal)
plot(d$signal, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(d$signal, d$LOS_d, method="spearman", use="complete.obs"), 3))
# text(1700, 0.78, bquote(rho == .(rho) ), cex=2)
abline(LOSdDpnIIseqm$coefficients[1], LOSdDpnIIseqm$coefficients[2], col="grey60", lwd = 3)

d_subset <- d[d$signal > 1000 & d$signal < 1100,]
d_subset_clean <- na.omit(d_subset)
subset_m <- lm(d_subset_clean$LOS_d~d_subset_clean$PC1)
plot(d_subset$PC1, d_subset$LOS_d, col=ifelse(d_subset$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=1, cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(d_subset$PC1, d_subset$LOS_d, method="spearman", use="complete.obs"), 3))
# text(0.01, 0.81, bquote(rho == .(rho) ), cex=2)
abline(subset_m$coefficients[1], subset_m$coefficients[2], col="grey60", lwd = 3)

LOS_residuals <- get_residuals(d_clean$LOS_d, d_clean$signal)
PC1_residuals <- get_residuals(d_clean$PC1, d_clean$signal)
resid_LOS_PC1_m <-lm(LOS_residuals ~ PC1_residuals)
plot(PC1_residuals, LOS_residuals, col=ifelse(d_clean$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(PC1_residuals, LOS_residuals, method="spearman", use="complete.obs"), 3))
# text(0, -0.075, bquote(rho == .(rho) ), cex=2)
abline(resid_LOS_PC1_m$coefficients[1], resid_LOS_PC1_m$coefficients[2], col="grey60", lwd=3)

DpnIIseq_residuals <- get_residuals(d_clean$signal, d_clean$PC1)
LOS_residuals <- get_residuals(d_clean$LOS_d, d_clean$PC1)
resid_LOS_DpnIIseq_m <- lm(LOS_residuals ~ DpnIIseq_residuals)
plot(DpnIIseq_residuals, LOS_residuals, col=ifelse(d_clean$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5,
     xlim=c(-1000,1000), cex.lab=1.5, axes=FALSE)
axis(1, labels=FALSE) 
axis(2, labels=FALSE)
box()
# rho <- sprintf("%.2f",round(cor(DpnIIseq_residuals, LOS_residuals, method="spearman", use="complete.obs"), 3))
# text(400, -0.075, bquote(rho == .(rho) ), cex=2)
abline(resid_LOS_DpnIIseq_m$coefficients[1], resid_LOS_DpnIIseq_m$coefficients[2], col="grey60", lwd=3)
dev.off()


