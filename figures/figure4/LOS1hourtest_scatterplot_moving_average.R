eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA

LOS2h <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse2/LOS1hourtest/HBHiCK562DN10-1hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_LOS.bedGraph",
                  sep="\t", header=TRUE)
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse2/LOS1hourtest/half-life1hourtest/half-life_exponential_40kb_removed_outliers_range6Mb_LOS1hourtest.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))
dpn <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "DpnIIseq"))

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
LOS2h_chr2 <- LOS2h[LOS2h$chrom == "chr2",]
thalf_chr2 <- thalf[thalf$chrom == "chr2",]
dpn_chr2 <- dpn[dpn$chrom == "chr2",]

get_moving_average <- function(x, y, n, s, dec) {
  ## Args:
  #  x: unordered numeric vector to slide windows across
  #  y: unordered numeric values to mean over x window
  #  n: window size
  #  s: step size
  #  dec: number of decimal places to round for merge
  # Returns:
  # df: dataframe with columns:
  #   w: window centers
  #   mu: moving average
  if ((sum(is.na(x)) > 0) | (sum(is.na(y)) > 0)) {
    stop("Remove NAs")
  }
  if (length(x)!=length(y)) {
    stop("Unequal vector lengths")
  }
  xmin <- round(min(x, na.rm=TRUE), dec)
  xmax <- round(max(x, na.rm=TRUE), dec)
  w <- seq(xmin, xmax, by=s)
  mu <- c()
  for (i in w) {
    m <- mean(y[x > (i-(n/2)) & x < (i + (n/2))])
    mu <- c(mu, m)
  }
  df <- data.frame(w, mu)
  return(df)
}

get_ma_residuals <- function(o, df, dec) {
  # Get moving average residuals
  ## Args:
  #  o: original dataframe with cols:
  #   column1: chrom
  #   column2: start
  #   column3: end
  #   column4: numeric vector original x-axis
  #   column5:  numeric vector original y-axis
  #  df: dataframe output by get_moving_average
  #      function
  #  dec: number of decimal places to round for merge
  # Returns:
  #  r: dataframe of residuals and locations
  colnames(o) <- c("chrom", "start", "end", "w", "y")
  o["w"] <- round(o["w"], dec)
  df["w"] <- round(df["w"], dec)
  m <- merge(o, df, by="w")
  m_resid <- m$y-m$mu
  r <- cbind(m[c("chrom", "start", "end")], m_resid)
  return(r)
}


d <- data.frame(eigen_chr2[1:4], LOS2h_chr2["LOS"],
                thalf_chr2["hl"], dpn_chr2["DpnIIseq"])

d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$hl[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$DpnIIseq[which(is.na(d), arr.ind = TRUE)[,1]] <- NA

d_clean <- na.omit(d)

ma <- get_moving_average(d_clean$DpnIIseq, d_clean$LOS, 200, 1, 0)
LOS_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "DpnIIseq", "LOS")], ma, 0)
colnames(LOS_r) <- c("chrom", "start", "end", "LOS_residuals")
ma <- get_moving_average(d_clean$DpnIIseq, d_clean$PC1, 200, 1, 0)
PC1_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "DpnIIseq", "PC1")], ma, 0)
colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
ma <- get_moving_average(d_clean$DpnIIseq, d_clean$hl, 200, 1, 0)
hl_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "DpnIIseq", "hl")], ma, 0)
colnames(hl_r) <- c("chrom", "start", "end", "hl_residuals")

m <- merge(d_clean, LOS_r, by=c("chrom", "start", "end"))
m <- merge(m, PC1_r, by=c("chrom", "start", "end"))
m <- merge(m, hl_r, by=c("chrom", "start", "end"))
om <- m[order(m$start),]

pdf("LOS1hourtest_LOS2hresid_vs_PC1resid_control_DpnIIseq1h_chr2_moving_average.pdf", height=8.25, width=6)
plot(om$PC1_residuals, om$LOS_residuals, col=ifelse(om$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab="2h LOS residuals", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(om$PC1_residuals, om$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0.01, -0.2, bquote(rho == .(rho) ), cex=2)
dev.off()

pdf("LOS1hourtest_thalfresid_vs_PC1resid_control_DpnIIseq1h_chr2_moving_average.pdf", height=8.25, width=6)
plot(om$PC1_residuals, om$hl_residuals, col=ifelse(om$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab=bquote("t"[1/2]~" residuals"), cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(om$PC1_residuals, om$hl_residuals, method="spearman", use="complete.obs"), 3))
text(-0.02, -15, bquote(rho == .(rho) ), cex=2)
dev.off()


