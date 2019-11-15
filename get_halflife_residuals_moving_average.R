source("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/digest-Hi-C/my_functions.R")
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))
dseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                   sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))


df <- cbind(thalf, dseq["DpnIIseq"])
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
hl_r <- get_ma_residuals(df_clean[c("chrom", "start", "end", "DpnIIseq", "hl")], ma)
colnames(hl_r) <- c("chrom", "start", "end", "hl_residuals")
d <- merge(df[1:3], hl_r, by=c("chrom", "start", "end"), all=TRUE)
om <- OrderChromStart(d)

write.table(om, "half-life_filter1000_timecourse1_residuals_moving_average_from_1hDpnIIseq_40kb.bedGraph", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)


png("halflife_vs_DpnIIseq1h_moving_average_zoom.png", width=2500, height=2500, res=300)
plot(df_clean$DpnIIseq, df_clean$hl, pch=20, col=rgb(.118,.565,1, alpha=0.05),
     ylab=bquote("t"[1/2]~" (minutes)"), xlab="DpnII-seq signal (1h)",
     xlim=c(200,3000), ylim=c(40,130))
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=3)
dev.off()



