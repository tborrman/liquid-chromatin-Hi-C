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



ma <- get_moving_average(d_clean$signal, d_clean$p, 50)
los_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal", "p")], ma)
colnames(los_r) <- c("chrom", "start", "end", "LOS_residuals")
d <- merge(d, los_r, by=c("chrom", "start", "end"), all=TRUE)
d <- d[order(d$start),]


pdf("fig3_tracks_moving_average_fatI_50w_loess.pdf", height = 4, width = 6)
par(mfrow=c(4,1), mar=c(2, 6, 0, 2) + 0.1)
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
plot(d$x/1000000, d$p, type="l", col="cyan3", 
      ylim=c(0.70, 1), xlab = "", ylab = "LOS FatI (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$signal, type="l", col="gray60",
      ylim=c(0,500), xlab = "", ylab= "FatI-seq\nsignal", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$LOS_residuals, type="l", col="black",
      ylim=c(-0.09, 0.09), xlab = "", ylab="LOS FatI (%)\n residuals", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()


