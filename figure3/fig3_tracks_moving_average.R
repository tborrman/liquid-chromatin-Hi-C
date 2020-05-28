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

ma <- get_moving_average(d_clean$signal, d_clean$LOS_d, 200)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "LOS_d")], ma)
colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
d <- merge(d, los_r, by=c("chrom", "start", "end"), all=TRUE)
d <- d[order(d$start),]


pdf("fig3_tracks_moving_average.pdf", height = 6, width = 6)
par(mfrow=c(5,1), mar=c(2, 4, 0, 2) + 0.1)
plot(d$start/1000000, d$PC1, type="n",
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
plot(d$start/1000000, d$LOS_h, type="l", col="red", 
     ylim=c(-0.1, 0.2), xlab = "", ylab= "LOS HindIII (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$LOS_d, type="l", col="cyan", 
      ylim=c(0.78, 0.95), xlab = "", ylab = "LOS DpnII (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$signal, type="l", col="gray60",
      ylim=c(400,2400), xlab = "", ylab= "DpnII-seq signal", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$LOS_residuals, type="l", col="black",
      ylim=c(-0.07, 0.05), xlab = "", ylab="LOS DpnII (%) residuals", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()


