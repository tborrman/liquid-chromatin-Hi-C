source("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/digest-Hi-C/my_functions.R")
eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))

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


pdf("fig4_thalf_residual_tracks_supp_moving_average.pdf", height = 7, width = 6)
par(mfrow=c(11,1), mar=c(2, 4, 0, 2) + 0.1)

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
  df <- cbind(thalf, d["DpnIIseq"], eigen["PC1"])
  clean_df <- na.omit(df)
  
  ma <- get_moving_average(clean_df$DpnIIseq, clean_df$hl, 200, 1, 0)
  hl_r <- get_ma_residuals(clean_df[c("chrom", "start", "end", "DpnIIseq", "hl")], ma, 0)
  colnames(hl_r) <- c("chrom", "start", "end", "hl_residuals")
  ma <- get_moving_average(clean_df$DpnIIseq, clean_df$PC1, 200, 1, 0)
  PC1_r <- get_ma_residuals(clean_df[c("chrom", "start", "end", "DpnIIseq", "PC1")], ma, 0)
  colnames(PC1_r) <- c("chrom", "start", "end", "PC1_residuals")
  m <- merge(clean_df, hl_r, by=c("chrom", "start", "end"))
  m <- merge(m, PC1_r, by=c("chrom", "start", "end"))
  om <- OrderChromStart(m)
  om_chr2 <- om[om$chrom == "chr2",]

    if (i == 11) {
    plot(om_chr2$start/1000000, om_chr2$hl_residuals, type="l", col="dodgerblue", 
        xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.1, ylim=c(-35,35))
    axis(1, lwd=2, cex.axis=1) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
  }
  else {
    plot(om_chr2$start/1000000, om_chr2$hl_residuals, type="l", col="dodgerblue", 
         xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.1, ylim=c(-35,35))
    axis(1, lwd=2, cex.axis=1, labels=FALSE) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
  }

  # Correlation genome wide
  print(cor(om$hl_residuals, om$PC1_residuals, method="spearman", use="complete.obs"))
}

dev.off()

