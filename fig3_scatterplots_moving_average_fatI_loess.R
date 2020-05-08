eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
LOS_FatI <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/FatI/LOS6Mb/HB-HiC-K562-Fat4-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.44_range6Mb_LOS_removed_outliers.bedGraph",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_d"))
FatIseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/FatI/FatI-seq/HB-FatISeq-4h-K562-R1_S9_L005_copy_correct_coverage_40kb_removed_outliers.bedGraph",
                      sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "signal"))

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
LOS_FatI_chr2 <- LOS_FatI[LOS_FatI$chrom == "chr2",]
FatIseq_chr2 <- FatIseq[FatIseq$chrom == "chr2",]


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

get_loess_prediction <- function(x, y) {
  # Calculate loess predictions to smooth signal
  ## Args:
  #  x: x-axis data (predictor variable)
  #  y: y-axis data (response variable)
  # Returns:
  #  d: data.frame of x and p
  #     x: x-axis for predictions (includes NAs)
  #     p: loess predictions
  x[is.na(y)] <- NA
  lo_obj <- loess(y ~ x, span = 0.01)
  p <- predict(lo_obj, x)
  d <- data.frame(x, p)
  return(d)
}

d <- data.frame(eigen_chr2[1:4], LOS_FatI_chr2["LOS_d"], FatIseq_chr2["signal"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_d[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA

# Loess
pred_df <- get_loess_prediction(d$start, d$LOS_d)
d <- cbind(d, pred_df)
d_clean <- na.omit(d)

pdf("fig3_scatterplots_moving_average_FatI_loess.pdf", height=9, width=5)
par(mfrow=c(3,2), mar=c(5, 5, 2, 2) + 0.1)

plot(d$PC1, d$p, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS FatI (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$PC1, d$p, method="spearman", use="complete.obs"), 3))
text(0.01, 0.88, bquote(rho == .(rho) ), cex=2)

ma <- get_moving_average(d_clean$signal, d_clean$p, 50, 1, 0)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "p")], ma, 0)
plot(d$signal, d$p, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="FatI-seq signal", ylab="LOS FatI (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$signal, d$p, method="spearman", use="complete.obs"), 3))
text(300, 0.88, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

d_subset <- d[d$signal > 75 & d$signal < 100,]
plot(d_subset$PC1, d_subset$p, col=ifelse(d_subset$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS FatI (%)", cex=1, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d_subset$PC1, d_subset$p, method="spearman", use="complete.obs"), 3))
text(0.01, 0.88, bquote(rho == .(rho) ), cex=2)

ma <- get_moving_average(d_clean$signal, d_clean$PC1, 50, 1, 0)
PC1_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "PC1")], ma, 0)
ma <- get_moving_average(d_clean$PC1, d_clean$signal, .004, .00001, 5)
FatI_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "signal")], ma, 5)
ma <- get_moving_average(d_clean$PC1, d_clean$p, .004, .00001, 5)
los_r_pc1 <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "p")], ma, 5)

colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
df_residuals <- merge(d_clean, los_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, PC1_r, by=c("chrom", "start", "end"))

plot(df_residuals$PC1_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab="LOS FatI (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(df_residuals$PC1_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0.01, -0.04, bquote(rho == .(rho) ), cex=2)

colnames(los_r_pc1) <- c("chrom", "start", "end", "LOS_residuals")
colnames(FatI_r) <- c("chrom", "start", "end", "FatI_residuals")
df_residuals <- merge(d_clean, los_r_pc1, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, FatI_r, by=c("chrom", "start", "end"))

plot(df_residuals$FatI_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="FatI-seq residuals", ylab="LOS FatI (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(df_residuals$FatI_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(150, -0.03, bquote(rho == .(rho) ), cex=2)
dev.off()


