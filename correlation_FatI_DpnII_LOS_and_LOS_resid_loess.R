eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
LOS_FatI <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/FatI/LOS6Mb/HB-HiC-K562-Fat4-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.44_range6Mb_LOS_removed_outliers.bedGraph",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "FatI_LOS"))
FatIseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/FatI/FatI-seq/HB-FatISeq-4h-K562-R1_S9_L005_copy_correct_coverage_40kb_removed_outliers.bedGraph",
                      sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "FatI_signal"))
LOS_DpnII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/FatI/LOS6Mb/HB-HiC-K562-DN4-R4-filter1000__hg19__genome__C-40000-iced_scaleBy_0.45_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "DpnII_LOS"))
DpnIIseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN4hR1_S8_L002_copy_correct_coverage_40kb_removed_outliers.bedGraph",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "DpnII_signal"))


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

d <- data.frame(eigen[1:4], LOS_FatI["FatI_LOS"], FatIseq["FatI_signal"],
                LOS_DpnII["DpnII_LOS"], DpnIIseq["DpnII_signal"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$FatI_LOS[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$FatI_signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$DpnII_LOS[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$DpnII_signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA

# Loess
p_DpnII <- c()
p_FatI <- c()
for (chrom in c(paste("chr", 1:22,sep=""), "chrX")) {
  chrom_start <- d$start[d$chrom == chrom]
  chrom_DpnII_LOS <- d$DpnII_LOS[d$chrom == chrom]
  chrom_FatI_LOS <- d$FatI_LOS[d$chrom  == chrom]
  chrom_pred_df <- get_loess_prediction(chrom_start, chrom_DpnII_LOS)
  p_DpnII <- c(p_DpnII, chrom_pred_df$p)
  chrom_pred_df <- get_loess_prediction(chrom_start, chrom_FatI_LOS)
  p_FatI <- c(p_FatI, chrom_pred_df$p)
}

d <- cbind(d, p_DpnII, p_FatI)
d_clean <- na.omit(d)


png("correlation_FatI_DpnII_LOS_loess.png", height=2500, width=1500, res=300)
par(mfrow=c(3,2), mar=c(5, 4.2, 0.2, 2) + 0.1)

plot(d_clean$p_DpnII, d_clean$p_FatI, col=ifelse(d_clean$PC1 > 0, rgb(1,0,0,0.1), rgb(0,0,1,0.1)),
     pch=20, xlab="LOS DpnII (%)", ylab="LOS FatI (%)", cex=0.5, cex.lab=1.5,
     xlim=c(0.65, 1.0), ylim=c(0.7, 1.0))
rho <- sprintf("%.2f",round(cor(d_clean$p_DpnII, d_clean$p_FatI, method="spearman", use="complete.obs"), 3))
text(0.9, 0.75, bquote(rho == .(rho) ), cex=2)

plot(d_clean$DpnII_signal, d_clean$FatI_signal, col=ifelse(d_clean$PC1 > 0, rgb(1,0,0,0.1), rgb(0,0,1,0.1)),
     pch=20, xlab="DpnII-seq signal", ylab="FatI-seq signal", cex=0.5, cex.lab=1.5,
     xlim=c(0, 3500), ylim=c(0, 650))
rho <- sprintf("%.2f",round(cor(d_clean$DpnII_signal, d_clean$FatI_signal, method="spearman", use="complete.obs"), 3))
text(2800, 100, bquote(rho == .(rho) ), cex=2)

ma <- get_moving_average(d_clean$DpnII_signal, d_clean$p_DpnII, 200, 1, 0)
DpnII_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "DpnII_signal", "p_DpnII")], ma, 0)
plot(d_clean$DpnII_signal, d_clean$p_DpnII, col=ifelse(d_clean$PC1 > 0, rgb(1,0,0,0.1), rgb(0,0,1,0.1)),
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5,
     xlim=c(0,3500), ylim=c(0.65, 1.0))
rho <- sprintf("%.2f",round(cor(d_clean$DpnII_signal, d_clean$p_DpnII, method="spearman", use="complete.obs"), 3))
text(2800, 0.7, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

ma <- get_moving_average(d_clean$FatI_signal, d_clean$p_FatI, 50, 1, 0)
FatI_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "FatI_signal", "p_FatI")], ma, 0)
plot(d_clean$FatI_signal, d_clean$p_FatI, col=ifelse(d_clean$PC1 > 0, rgb(1,0,0,0.1), rgb(0,0,1,0.1)),
     pch=20, xlab="FatI-seq signal", ylab="LOS FatI (%)", cex=0.5, cex.lab=1.5,
     xlim=c(0,650), ylim=c(0.7, 1.0))
rho <- sprintf("%.2f",round(cor(d_clean$FatI_signal, d_clean$p_FatI, method="spearman", use="complete.obs"), 3))
text(500, 0.75, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

colnames(DpnII_r) <- c("chrom", "start", "end", "DpnII_resid")
colnames(FatI_r) <- c("chrom", "start", "end", "FatI_resid")

df_residuals <- merge(d_clean, DpnII_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, FatI_r, by=c("chrom", "start", "end"))

plot(df_residuals$DpnII_resid, df_residuals$FatI_resid, col=ifelse(df_residuals$PC1 > 0, rgb(1,0,0,0.1), rgb(0,0,1,0.1)),
     pch=20, xlab="LOS DpnII (%) residuals", ylab="LOS FatI (%) residuals", cex=0.5, cex.lab=1.5,
     xlim=c(-0.2, 0.12), ylim=c(-0.2, 0.12))
rho <- sprintf("%.2f",round(cor(df_residuals$DpnII_resid, df_residuals$FatI_resid, method="spearman", use="complete.obs"), 3))
text(0.0, -0.15, bquote(rho == .(rho) ), cex=2)

dev.off()


