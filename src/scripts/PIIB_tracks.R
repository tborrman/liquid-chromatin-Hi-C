eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
LOS_DpnII_TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBTD-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.59_range6Mb_LOS.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_TD"))
LOS_DpnII_NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBNT-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.65_range6Mb_LOS.bedGraph",
                           sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_NT"))

eigen_chr2 <- eigen[eigen$chrom == "chr2",]
LOS_DpnII_TD_chr2 <- LOS_DpnII_TD[LOS_DpnII_TD$chrom == "chr2",]
LOS_DpnII_NT_chr2 <- LOS_DpnII_NT[LOS_DpnII_NT$chrom == "chr2",]


get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}


d <- data.frame(eigen_chr2[1:4], LOS_DpnII_TD_chr2["LOS_TD"],
                LOS_DpnII_NT_chr2["LOS_NT"])
# NA bins
d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_TD[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
d$LOS_NT[which(is.na(d), arr.ind = TRUE)[,1]] <- NA


d_clean <- na.omit(d)



pdf("PIIB_tracks.pdf", height = 4, width = 6)
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
plot(d$start/1000000, d$LOS_NT, type="l", col="red", 
     ylim=c(0.65, 1), xlab = "", ylab= "Not Treated\nLOS DpnII (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$LOS_TD, type="l", col="cyan", 
      ylim=c(0.65, 1), xlab = "", ylab = "Pol II Block\nLOS DpnII (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$LOS_NT, type="l", col="red", 
     ylim=c(0.65, 1), xlab = "", ylab= "Not Treated\nLOS DpnII (%)", axes=FALSE, lwd=0.5)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()


