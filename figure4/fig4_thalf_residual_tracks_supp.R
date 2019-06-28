library(RColorBrewer)
eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
thalf_chr2 <- thalf[thalf$chrom == "chr2",]

remove_outliers <- function(x) {
  s <- sd(x, na.rm = TRUE)
  m <- mean(x, na.rm = TRUE)
  x[x > m + (s*4)] <- NA
  x[x < m - (s*4)] <- NA
  return(x)
}

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}

bluescol <-  colorRampPalette(brewer.pal(n = 9, name ="Blues")[3:9])(11)

pdf("fig4_thalf_residual_tracks_supp.pdf", height = 7, width = 6)
par(mfrow=c(11,1), mar=c(2, 4, 0, 2) + 0.1)
# plot(eigen_chr2$start/1000000, eigen_chr2$PC1, type="n",
#       xlab = "", ylab= "PC1", axes=FALSE)
# A <- eigen_chr2$PC1
# A[is.na(A)] <- 0
# A[A<0] <- 0
# polygon(c(eigen_chr2$start/1000000, rev(eigen_chr2$start/1000000)), 
#         c(A, rep(0, length(A))), col="red",
#         border=NA)
# B <- eigen_chr2$PC1
# B[is.na(B)] <- 0
# B[B>0] <- 0
# polygon(c(eigen_chr2$start/1000000, rev(eigen_chr2$start/1000000)), 
#         c(B, rep(0, length(B))), col="blue",
#         border=NA)
# axis(1, lwd=2, cex.axis=1, labels=FALSE) 
# axis(2, lwd=2, cex.axis=1)
# box(bty="l", lwd=2)
# plot(thalf_chr2$start/1000000, thalf_chr2$hl, type="l", col="dodgerblue", 
#      xlab = "", ylab= bquote("t"[1/2]), axes=FALSE, lwd=0.5)
# axis(1, lwd=2, cex.axis=1, labels=FALSE) 
# axis(2, lwd=2, cex.axis=1)
# box(bty="l", lwd=2)

f <- list(c("HBDpSeqK562-DN5mR1_S1_L001_copy_correct_coverage_40kb.bed", "5 min"),
  c("HBDpSeqK562-DN15mR1_S2_L001_copy_correct_coverage_40kb.bed", " 15 min"),
  c("HBDpSeqK562-DN30mR1_S3_L001_copy_correct_coverage_40kb.bed", "30 min"),
  c("HBDpSeqK562-DN45mR1_S4_L001_copy_correct_coverage_40kb.bed", "45 min"),
  c("HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed", "60 min"),
  c("HBDpSeqK562-DN75mR1_S9_L003_copy_correct_coverage_40kb.bed", "75 min"),
  c("HBDpSeqK562-DN90mR1_S10_L003_copy_correct_coverage_40kb.bed", "90 min"),
  c("HBDpSeqK562-DN2hR1_S6_L002_copy_correct_coverage_40kb.bed", "120 min"),
  c("HBDpSeqK562-DN3hR1_S7_L002_copy_correct_coverage_40kb.bed", "180 min"),
  c("HBDpSeqK562-DN4hR1_S8_L002_copy_correct_coverage_40kb.bed", "240 min"),
  c("HBDpSeqK562-DN16hR1_S11_L003_copy_correct_coverage_40kb.bed", "960 min"))

for (i in 1:11) {
  d <- read.table(paste("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/",f[[i]][1], sep=""),
                  sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))
  d_chr2 <- d[d$chrom == "chr2",]
  df <- cbind(thalf_chr2, d_chr2["DpnIIseq"], eigen_chr2["PC1"])
  clean_df <- na.omit(df)
  hl_resid <- get_residuals(clean_df$hl, clean_df$DpnIIseq)
  hl_res_clean <- remove_outliers(hl_resid)
  clean_df <- cbind(clean_df, hl_res_clean)
  m <- merge(d_chr2[1:3], clean_df, by=c("chrom", "start", "end"), all=TRUE)
  om <- m[order(m$start),]
  if (i == 11) {
    plot(om$start/1000000, om$hl_res_clean, type="l", col=bluescol[i], 
        xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.5)
    axis(1, lwd=2, cex.axis=1) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
  }
  else {
    plot(om$start/1000000, om$hl_res_clean, type="l", col=bluescol[i], 
         xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.5)
    axis(1, lwd=2, cex.axis=1, labels=FALSE) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
  }

  print(cor(thalf_chr2["hl"], d_chr2["DpnIIseq"], method="spearman", use="complete.obs"))
  print(cor(clean_df["PC1"], hl_resid,  method="spearman", use="complete.obs"))
}

dev.off()

