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

get_moving_average <- function(x, y, n) {
  ## Args:
  #  x: unordered numeric vector to slide windows across
  #  y: unordered numeric values to mean over x window
  #  n: window size
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
  xmin <- round(min(x, na.rm=TRUE))
  xmax <- round(max(x, na.rm=TRUE))
  w <- seq(xmin, xmax)
  mu <- c()
  for (i in w) {
    m <- mean(y[x > (i-(n/2)) & x < (i + (n/2))])
    mu <- c(mu, m)
  }
  df <- data.frame(w, mu)
  return(df)
}

get_ma_residuals <- function(o, df) {
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
  # Returns:
  #  r: dataframe of residuals and locations
  colnames(o) <- c("chrom", "start", "end", "w", "y")
  o["w"] <- round(o["w"])
  print(head(o))
  m <- merge(o, df, by="w")
  print(head(m))
  m_resid <- m$y-m$mu
  r <- cbind(m[c("chrom", "start", "end")], m_resid)
  print(head(r))
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



pdf("scatterplot_get_moving_average.pdf", height=9, width=5)
par(mfrow=c(3,2), mar=c(5, 4, 0, 2) + 0.1)


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

ma <- get_moving_average(d_clean$signal, d_clean$LOS_d, 100)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "LOS_d")], ma)

plot(d$signal, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$signal, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(1700, 0.78, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

ma <- get_moving_average(d_clean$signal, d_clean$PC1, 100)
PC1_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "PC1")], ma)
plot(d$signal, d$PC1, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="DpnII-seq signal", ylab="PC1", cex=0.5, cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(d$signal, d$PC1, method="spearman", use="complete.obs"), 3))
text(1700, -0.015, bquote(rho == .(rho) ), cex=2)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
df_residuals <- merge(d_clean, los_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, PC1_r, by=c("chrom", "start", "end"))

plot(df_residuals$PC1_residuals, df_residuals$LOS_residuals, col=ifelse(df_residuals$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1 residuals", ylab="LOS DpnII (%) residuals", cex=0.5, 
     cex.lab=1.5)
rho <- sprintf("%.2f",round(cor(df_residuals$PC1_residuals, df_residuals$LOS_residuals, method="spearman", use="complete.obs"), 3))
text(0, -0.075, bquote(rho == .(rho) ), cex=2)
#abline(resid_LOS_PC1_m$coefficients[1], resid_LOS_PC1_m$coefficients[2], col="grey60", lwd=1.5)

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


png("scatterplot_get_moving_average.png", height=3500, width=1600, res=300)
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


