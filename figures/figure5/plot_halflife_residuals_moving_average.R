source("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/digest-Hi-C/my_functions.R")
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))
eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
dseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                   sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))


df <- cbind(thalf, eigen["PC1"], dseq["DpnIIseq"])
df_clean <- na.omit(df)

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
  print(tail(o))
  m <- merge(o, df, by="w")
  print(tail(m))
  m_resid <- m$y-m$mu
  r <- cbind(m[c("chrom", "start", "end")], m_resid)
  print(tail(r))
  return(r)
}

ma <- get_moving_average(df_clean$DpnIIseq, df_clean$hl, 200)

png("halflife_vs_DpnIIseq1h_moving_average_zoom_AB.png", width=2500, height=2500, res=300)
plot(df_clean$DpnIIseq, df_clean$hl, pch=20, col=ifelse(df_clean$PC1 > 0, "red", "blue"),
     ylab=bquote("t"[1/2]~" (minutes)"), xlab="DpnII-seq signal (1h)",
     xlim=c(200,3000), ylim=c(40,130))
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=3)
dev.off()


png("halflife_vs_DpnIIseq1h_moving_average_zoom_AB_density_alpha0.1.png", width=2500, height=2500, res=300)
A <- df_clean[df_clean$PC1 > 0,]
B <- df_clean[df_clean$PC1 < 0,]
#steps <- seq(1, length(B$chrom) +100, 100)
Adone = FALSE
Bdone = FALSE
stepsize <- 10
steps <- seq(1, length(B$chrom) + stepsize, stepsize)
plot(0, ylab=bquote("t"[1/2]~" (minutes)"), xlab="DpnII-seq signal (1h)", 
     xlim=c(200,3000), ylim=c(40,130))
for (i in 1:length(steps)) {
  if (! Adone) {
    if (i > length(A$chrom)) {
      points(A$DpnIIseq[steps[i]:length(A$DpnIIseq)], A$hl[steps[i]:length(A$DpnIIseq)], pch=20, col=rgb(1,0,0, alpha=0.1))
      Adone = TRUE
    }
    else {
      points(A$DpnIIseq[steps[i]:(steps[i]+(stepsize-1))], A$hl[steps[i]:(steps[i]+(stepsize-1))], pch=20, col=rgb(1,0,0, alpha=0.1))
    }
  }
  if (! Bdone) {
    if (i > length(B$chrom)) {
      points(B$DpnIIseq[steps[i]:length(B$DpnIIseq)], B$hl[steps[i]:length(B$DpnIIseq)], pch=20, col=rgb(0,0,1, alpha=0.1))
      Bdone = TRUE
    }
    else {
      points(B$DpnIIseq[steps[i]:(steps[i]+(stepsize-1))], B$hl[steps[i]:(steps[i]+(stepsize-1))], pch=20, col=rgb(0,0,1, alpha=0.1))
    }
  }
}
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=3)
dev.off()

cor(df_clean$DpnIIseq, df_clean$hl, method="spearman", use="complete.obs")

# hl residulas vs DpnIIseq
hl_r <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life_residuals/half-life_filter1000_timecourse1_residuals_moving_average_from_1hDpnIIseq_40kb.bedGraph",
                   sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl_residuals")) 
df <- cbind(hl_r, eigen["PC1"], dseq["DpnIIseq"])
df_clean <- na.omit(df)

png("halflife_residuals_vs_DpnIIseq1h_moving_average_zoom_AB_density_alpha0.1.png", width=2500, height=2500, res=300)
A <- df_clean[df_clean$PC1 > 0,]
B <- df_clean[df_clean$PC1 < 0,]
#steps <- seq(1, length(B$chrom) +100, 100)
Adone = FALSE
Bdone = FALSE
stepsize <- 10
steps <- seq(1, length(B$chrom) + stepsize, stepsize)
plot(0, ylab=bquote("t"[1/2]~" residuals (minutes)"), xlab="DpnII-seq signal (1h)", 
     xlim=c(200,3000), ylim=c(-50,50))
for (i in 1:length(steps)) {
  if (! Adone) {
    if (i > length(A$chrom)) {
      points(A$DpnIIseq[steps[i]:length(A$DpnIIseq)], A$hl_residuals[steps[i]:length(A$DpnIIseq)], pch=20, col=rgb(1,0,0, alpha=0.1))
      Adone = TRUE
    }
    else {
      points(A$DpnIIseq[steps[i]:(steps[i]+(stepsize-1))], A$hl_residuals[steps[i]:(steps[i]+(stepsize-1))], pch=20, col=rgb(1,0,0, alpha=0.1))
    }
  }
  if (! Bdone) {
    if (i > length(B$chrom)) {
      points(B$DpnIIseq[steps[i]:length(B$DpnIIseq)], B$hl_residuals[steps[i]:length(B$DpnIIseq)], pch=20, col=rgb(0,0,1, alpha=0.1))
      Bdone = TRUE
    }
    else {
      points(B$DpnIIseq[steps[i]:(steps[i]+(stepsize-1))], B$hl_residuals[steps[i]:(steps[i]+(stepsize-1))], pch=20, col=rgb(0,0,1, alpha=0.1))
    }
  }
}
dev.off()

cor(df_clean$DpnIIseq, df_clean$hl_residuals, method="spearman", use="complete.obs")
