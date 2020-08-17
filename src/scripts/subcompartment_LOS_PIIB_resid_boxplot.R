
TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBTD-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.59_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_TD"))
NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBNT-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.65_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_NT"))
DpnIIseq_NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/DpnII-seq/DL-20200103-DpnIIseq-K562-PIIBNT-DpnII-R1-T1_S6_L002_copy_correct_coverage_removed_outliers.bedGraph",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "signal_NT"))
DpnIIseq_TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/DpnII-seq/DL-20200103-DpnIIseq-K562-PIIBTD-DpnII-R1-T1_S7_L002_copy_correct_coverage_removed_outliers.bedGraph",
                       sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "signal_TD"))
sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                  sep="\t", header=TRUE)

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

d <- data.frame(TD[1:4], NT["LOS_NT"], DpnIIseq_TD["signal_TD"],
                DpnIIseq_NT["signal_NT"])
# NA bins
d$LOS_TD[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_NT[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$signal_TD[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$signal_NT[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d_clean <- na.omit(d)

# Treated
ma <- get_moving_average(d_clean$signal_TD, d_clean$LOS_TD, 50, 1, 0)
TD_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal_TD", "LOS_TD")], ma, 0)
plot(d_clean$signal_TD, d_clean$LOS_TD, col="blue",
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

# Not Treated
ma <- get_moving_average(d_clean$signal_NT, d_clean$LOS_NT, 50, 1, 0)
NT_r <- get_ma_residuals(d_clean[c("chrom", "start", "end", "signal_NT", "LOS_NT")], ma, 0)
plot(d_clean$signal_NT, d_clean$LOS_NT, col="blue",
     pch=20, xlab="DpnII-seq signal", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5)
lines(ma$w, ma$mu, pch=20, col="grey60", lwd=1.5)

colnames(TD_r) <- c("chrom", "start", "end", "TD_res")
colnames(NT_r) <- c("chrom", "start", "end", "NT_res")

df_residuals <- merge(sub, TD_r, by=c("chrom", "start", "end"))
df_residuals <- merge(df_residuals, NT_r, by=c("chrom", "start", "end"))

df_residuals <- na.omit(df_residuals)

# Stratify
A1 <- df_residuals[df_residuals$sub=="A1",]
A2 <- df_residuals[df_residuals$sub=="A2",]
B1 <- df_residuals[df_residuals$sub=="B1",]
B2 <- df_residuals[df_residuals$sub=="B2",]
B3 <- df_residuals[df_residuals$sub=="B3",]

pdf('subcompartment_LOS_PIIB_resid_boxplot.pdf', width=5, height=6)
par(mar=c(6,4,1,1) + 0.1)
boxplot(A1$NT_res, A1$TD_res, A2$NT_res, A2$TD_res,
        B1$NT_res, B1$TD_res, B2$NT_res, B2$TD_res,
        B3$NT_res, B3$TD_res,
        col=c("red", "red", "hotpink", "hotpink",
              "purple", "purple", "cyan", "cyan",
              "blue", "blue"), outline=FALSE,
        names=c("Not Treated", "Pol II Block", "Not Treated", "Pol II Block",
                "Not Treated", "Pol II Block","Not Treated", "Pol II Block",
                "Not Treated", "Pol II Block"),
        ylab = "LOS DpnII (%) residuals", las=2)
# text(1,median(A1_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A1_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(2,median(A1_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A1_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(3,median(A2_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A2_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(4,median(A2_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A2_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(5,median(B1_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B1_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(6,median(B1_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B1_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(7,median(B2_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B2_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(8,median(B2_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B2_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(9,median(B3_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B3_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(10,median(B3_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B3_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
legend("bottomleft", c("A1", "A2", "B1", "B2", "B3"), 
       fill=c("red", "hotpink", "purple", "cyan", "blue"), bty="n")
dev.off()

t.test(A1$NT_res, A1$TD_res)
t.test(A2$NT_res, A2$TD_res)
t.test(B1$NT_res, B1$TD_res)
t.test(B2$NT_res, B2$TD_res)
t.test(B3$NT_res, B3$TD_res)





