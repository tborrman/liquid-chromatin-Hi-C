eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
FatI <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/restriction_enzymes/FatI_sites_40kb.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "FatI"))
HindIII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/restriction_enzymes/HindIII_sites_40kb.bedGraph",
                   sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "HindIII"))
DpnII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/restriction_enzymes/dpnII_sites_40kb.bedGraph",
                   sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "DpnII"))


eigen_chr2 <- eigen[eigen$chrom == "chr2",]
FatI_chr2 <- FatI[FatI$chrom == "chr2",]
HindIII_chr2 <- HindIII[HindIII$chrom == "chr2",]
DpnII_chr2 <- DpnII[DpnII$chrom == "chr2",]

d <- data.frame(eigen_chr2[1:4], FatI_chr2["FatI"], HindIII_chr2["HindIII"], DpnII_chr2["DpnII"])

pdf("restriction_site_tracks.pdf", height = 4, width = 10)
par(mfrow=c(4,1), mar=c(0, 6, 0, 2) + 0.1)
plot(d$start/1000000, d$PC1, type="n",
      xlab = "", ylab= "", axes=FALSE)
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
plot(d$start/1000000, d$HindIII, type="n", 
     ylim=c(7,28), xlab = "", ylab = "", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$HindIII, rep(0, length(d$HindIII))), col="red",
        border=NA)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$DpnII, type="n",
     ylim=c(65,175),xlab = "", ylab= "", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$DpnII, rep(0, length(d$DpnII))), col="cyan",
        border=NA)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$FatI, type="n",
     ylim=c(165,300),xlab = "", ylab="", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$FatI, rep(0, length(d$FatI))), col="orange",
        border=NA)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()


