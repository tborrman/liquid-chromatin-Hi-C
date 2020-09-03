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


d <- data.frame(eigen_chr2[1:4], LOS_DpnII_chr2["LOS_d"],
                LOS_HindIII_chr2["LOS_h"], DpnIIseq_chr2["signal"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_d[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_h[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d_clean <- na.omit(d)

pdf("fig3_scatterplots_moving_average.pdf", height=9, width=5)
par(mfrow=c(3,2), mar=c(5, 4, 0, 2) + 0.1)

plot(d$PC1, d$LOS_h, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS HindIII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_h, method="spearman", use="complete.obs"), 3))
text(0.01, -0.10, bquote(rho == .(rho) ), cex=2)

plot(d$PC1, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.78, bquote(rho == .(rho) ), cex=2)

ma <- get_moving_average(d_clean$signal, d_clean$LOS_d, 200, 1, 0)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "LOS_d")], ma, 0)
plot(d$signal, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$signal, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(1700, 0.78, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

d_subset <- d[d$signal > 1000 & d$signal < 1100,]
plot(d_subset$PC1, d_subset$LOS_d, col=ifelse(d_subset$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=1, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d_subset$PC1, d_subset$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.81, bquote(rho == .(rho) ), cex=2)

ma <- get_moving_average(d_clean$signal, d_clean$PC1, 200, 1, 0)
PC1_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "PC1")], ma, 0)
ma <- get_moving_average(d_clean$PC1, d_clean$signal, .004, .00001, 5)
DpnII_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "signal")], ma, 5)
ma <- get_moving_average(d_clean$PC1, d_clean$LOS_d, .004, .00001, 5)
los_r_pc1 <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "LOS_d")], ma, 5)

colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
df_residuals <- merge(d_clean, los_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, PC1_r, by=c("chrom", "start", "end"))

plot(df_residuals$PC1_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab="LOS DpnII (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(df_residuals$PC1_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0, -0.075, bquote(rho == .(rho) ), cex=2)

colnames(los_r_pc1) <- c("chrom", "start", "end", "LOS_residuals")
colnames(DpnII_r) <- c("chrom", "start", "end", "DpnII_residuals")
df_residuals <- merge(d_clean, los_r_pc1, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, DpnII_r, by=c("chrom", "start", "end"))

plot(df_residuals$DpnII_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq residuals", ylab="LOS DpnII (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(df_residuals$DpnII_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0, -0.075, bquote(rho == .(rho) ), cex=2)
dev.off()



png("fig3_scatterplots_moving_average.png", height=3500, width=1600, res=300)
par(mfrow=c(3,2), mar=c(5, 4, 2, 2) + 0.1)

plot(d$PC1, d$LOS_h, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_h, method="spearman", use="complete.obs"), 3))
# text(0.01, -0.10, bquote(rho == .(rho) ), cex=2)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()

plot(d$PC1, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_d, method="spearman", use="complete.obs"), 3))
# text(0.01, 0.78, bquote(rho == .(rho) ), cex=2)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()

ma <- get_moving_average(d_clean$signal, d_clean$LOS_d, 200, 1, 0)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "LOS_d")], ma, 0)
plot(d$signal, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(d$signal, d$LOS_d, method="spearman", use="complete.obs"), 3))
# text(1700, 0.78, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()

d_subset <- d[d$signal > 1000 & d$signal < 1100,]
plot(d_subset$PC1, d_subset$LOS_d, col=ifelse(d_subset$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=1, cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(d_subset$PC1, d_subset$LOS_d, method="spearman", use="complete.obs"), 3))
# text(0.01, 0.81, bquote(rho == .(rho) ), cex=2)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()

ma <- get_moving_average(d_clean$signal, d_clean$PC1, 200, 1, 0)
PC1_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "PC1")], ma, 0)
ma <- get_moving_average(d_clean$PC1, d_clean$signal, .004, .00001, 5)
DpnII_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "signal")], ma, 5)
ma <- get_moving_average(d_clean$PC1, d_clean$LOS_d, .004, .00001, 5)
los_r_pc1 <- get_ma_residuals(d_clean[c("chrom", "start", "end", "PC1", "LOS_d")], ma, 5)

colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
df_residuals <- merge(d_clean, los_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, PC1_r, by=c("chrom", "start", "end"))

plot(df_residuals$PC1_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, 
     cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(df_residuals$PC1_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
# text(0, -0.075, bquote(rho == .(rho) ), cex=2)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()

colnames(los_r_pc1) <- c("chrom", "start", "end", "LOS_residuals")
colnames(DpnII_r) <- c("chrom", "start", "end", "DpnII_residuals")
df_residuals <- merge(d_clean, los_r_pc1, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, DpnII_r, by=c("chrom", "start", "end"))

plot(df_residuals$DpnII_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="", ylab="", cex=0.5, 
     cex.lab=1.5, axes=FALSE)
# rho <- sprintf("%.2f",round(cor(df_residuals$DpnII_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
# text(0, -0.075, bquote(rho == .(rho) ), cex=2)
axis(1, labels=FALSE)
axis(2, labels=FALSE)
box()
dev.off()


