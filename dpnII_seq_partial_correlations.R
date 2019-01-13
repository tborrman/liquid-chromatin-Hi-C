df <- read.table("/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v5/feature_matrix_v5_40kb.txt", sep="\t", header=TRUE)
PC1 <- df$PCA_eigen1
LOS_df <- read.table("/cygwin64/home/Tyler/Research/digest/LOS_6Mb/C-40000/rangeLOS/HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_range6Mb_LOS_removed_outliers.bedGraph", sep="\t",header=FALSE, col.names=c("chrom", "start", "end", "LOS_4h"))
LOS <- LOS_df$LOS_4h
DpnIIseq_df <- read.table("/cygwin64/home/Tyler/Research/digest/DpnII-seq/HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_copy_correct_coverage_40kb_removed_outliers.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "dpnIIseq"))
DpnIIseq <- DpnIIseq_df$dpnIIseq
PC1[is.na(PC1)] <- NA
d <- data.frame(PC1, LOS, DpnIIseq)

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
      xlab="PC1", ylab="DpnII_LOS_4hours", cex.lab=1.5)
rho <- sprintf("%.3f", round(cor(d2_clean$PC1, d2_clean$LOS, method="spearman"), 3))
text(-0.015, 0.8, bquote(rho == .(rho) ), cex=2)
dev.off()

clean_LOS <- d2_clean$LOS
clean_DpnIIseq <- d2_clean$DpnIIseq
clean_PC1 <- d2_clean$PC1

LOS_residuals <- get_residuals(d2_clean$LOS, d2_clean$DpnIIseq)
PC1_residuals <- get_residuals(d2_clean$PC1, d2_clean$DpnIIseq)

LOSm <- lm(d2_clean$LOS~d2_clean$DpnIIseq)
PC1m <- lm(d2_clean$PC1~d2_clean$DpnIIseq)


png('DpnIIseq_vs_LOS4h_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(d2_clean$DpnIIseq, d2_clean$LOS, pch=21, col="black", bg="dodgerblue",
     xlab="DpnII-seq", ylab="DpnII_LOS_4hours", cex.lab=1.5)
rho <- sprintf("%.3f",round(cor(d2_clean$DpnIIseq, d2_clean$LOS, method="spearman"), 3))
text(20, 0.8, bquote(rho == .(rho) ), cex=2)
abline(LOSm$coefficients[1], LOSm$coefficients[2], col="red")
dev.off()

png('DpnIIseq_vs_PC1_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(d2_clean$DpnIIseq, d2_clean$PC1, pch=21, col="black", bg="dodgerblue",
     xlab="DpnII-seq", ylab="PC1", cex.lab=1.5)
rho <- sprintf("%.3f",round(cor(d2_clean$DpnIIseq, d2_clean$PC1, method="spearman"), 3))
text(20, 0.025,bquote(rho == .(rho) ), cex=2)
abline(PC1m$coefficients[1], PC1m$coefficients[2], col="red")
dev.off()

png('partial_correlation_PC1_vs_LOS4h_controlled_by_DpnII-seq.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(PC1_residuals, LOS_residuals, pch=21, col="black", bg="dodgerblue",
     xlab="PC1_residuals", ylab="DpnII_LOS_4hours_residuals", cex.lab=1.5)
rho <- sprintf("%.3f", round(cor(PC1_residuals, LOS_residuals, method="spearman"), 3))
text(-0.018,0.12, bquote(rho == .(rho) ), cex=2)
dev.off()

out_df <- data.frame(clean_LOS, clean_DpnIIseq, clean_PC1, LOS_residuals, PC1_residuals)
colnames(out_df) <- c("DpnII_4h_LOS", "DpnIIseq", "PC1", "DpnII_4h_LOS_residuals", "PC1_residuals")
write.table(out_df, "DpnII-seq_controlled_partial_correlations.txt", sep="\t", quote=FALSE, row.names=FALSE)

