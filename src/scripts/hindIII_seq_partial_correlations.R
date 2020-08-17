df <- read.table("/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v5/feature_matrix_v5_40kb.txt", sep="\t", header=TRUE)
PC1 <- df$PCA_eigen1
LOS_df <- read.table("/Users/tyler/Dropbox (UMass Medical School)/digest_092718/replicates/HN/HBHiCK562HN50-4hDp2__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS_removed_outliers.bedGraph", sep="\t",header=FALSE, col.names=c("chrom", "start", "end", "LOS_4h"))
LOS <- LOS_df$LOS_4h
HindIIIseq_df <- read.table("/Users/tyler/Dropbox (UMass Medical School)/digest_092718/hindIII/HindIII-seq/HBCRACKHiC-K562-HN-TD-R1_ACAGTG_L008_copy_correct_coverage_40kb_removed_outliers.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hindIIIseq"))
HindIIIseq <- HindIIIseq_df$hindIIIseq
PC1[is.na(PC1)] <- NA
d <- data.frame(PC1, LOS, HindIIIseq)

#Chrom2
d2 <- d[df$chrom == "chr2",]
d2_clean <- na.omit(d2) 

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}

png('PC1_vs_LOS4h_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(d2_clean$PC1, d2_clean$LOS, pch=21, col="black", bg="dodgerblue",
      xlab="PC1", ylab="HindIII_LOS_4hours", cex.lab=1.5)
rho <- sprintf("%.3f", round(cor(d2_clean$PC1, d2_clean$LOS, method="spearman"), 3))
text(-0.015, -0.6, bquote(rho == .(rho) ), cex=2)
dev.off()

clean_LOS <- d2_clean$LOS
clean_HindIIIseq <- d2_clean$HindIIIseq
clean_PC1 <- d2_clean$PC1

LOS_residuals <- get_residuals(d2_clean$LOS, d2_clean$HindIIIseq)
PC1_residuals <- get_residuals(d2_clean$PC1, d2_clean$HindIIIseq)

LOSm <- lm(d2_clean$LOS~d2_clean$HindIIIseq)
PC1m <- lm(d2_clean$PC1~d2_clean$HindIIIseq)


png('HindIIIseq_vs_LOS4h_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(d2_clean$HindIIIseq, d2_clean$LOS, pch=21, col="black", bg="dodgerblue",
     xlab="HindIII-seq", ylab="HindIII_LOS_4hours", cex.lab=1.5)
rho <- sprintf("%.3f",round(cor(d2_clean$HindIIIseq, d2_clean$LOS, method="spearman"), 3))
text(20, -0.8, bquote(rho == .(rho) ), cex=2)
abline(LOSm$coefficients[1], LOSm$coefficients[2], col="red")
dev.off()

png('HindIIIseq_vs_PC1_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(d2_clean$HindIIIseq, d2_clean$PC1, pch=21, col="black", bg="dodgerblue",
     xlab="HindIII-seq", ylab="PC1", cex.lab=1.5)
rho <- sprintf("%.3f",round(cor(d2_clean$HindIIIseq, d2_clean$PC1, method="spearman"), 3))
text(20, 0.025,bquote(rho == .(rho) ), cex=2)
abline(PC1m$coefficients[1], PC1m$coefficients[2], col="red")
dev.off()

png('partial_correlation_PC1_vs_LOS4h_controlled_by_HindIII-seq.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(PC1_residuals, LOS_residuals, pch=21, col="black", bg="dodgerblue",
     xlab="PC1_residuals", ylab="HindIII_LOS_4hours_residuals", cex.lab=1.5)
rho <- sprintf("%.3f", round(cor(PC1_residuals, LOS_residuals, method="spearman"), 3))
text(-0.018,-0.6, bquote(rho == .(rho) ), cex=2)
dev.off()

pdf('partial_correlation_PC1_vs_LOS4h_controlled_by_HindIII-seq_chr2.pdf', width=8, height=9)
par(mar=c(c(5, 5, 4, 2) + 0.1))
A_eigenR <- PC1_residuals[clean_PC1 > 0]
B_eigenR <- PC1_residuals[clean_PC1 < 0]
A_LOSR <- LOS_residuals[clean_PC1 > 0]
B_LOSR <- LOS_residuals[clean_PC1 < 0]

plot(A_eigenR, A_LOSR, pch=20, col="red",
     xlab="PC1 residuals", ylab="HindIII LOS 4hours residuals", cex.lab=1.5, cex=0.5,
     ylim=c(-0.25,0.25), xlim=c(-0.025, 0.03)
)
points(B_eigenR, B_LOSR, pch=20, col="blue", cex=0.5)
rho <- sprintf("%.3f", round(cor(PC1_residuals, LOS_residuals, method="spearman"), 3))
text(0.02,0.2, bquote(rho == .(rho) ), cex=2)

dev.off()



out_df <- data.frame(clean_LOS, clean_HindIIIseq, clean_PC1, LOS_residuals, PC1_residuals)
colnames(out_df) <- c("HindIII_4h_LOS", "HindIIIseq", "PC1", "HindIII_4h_LOS_residuals", "PC1_residuals")
write.table(out_df, "HindIII-seq_controlled_partial_correlations_R2_HN50-4hDp2.txt", sep="\t", quote=FALSE, row.names=FALSE)

