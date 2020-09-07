eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
slice1 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz1-K562-R1/bed/HB-Dp4h-siz1-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "size1"))
slice2 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz2-K562-R1/bed/HB-Dp4h-siz2-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                     sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "size2"))
slice3 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz3-K562-R1/bed/HB-Dp4h-siz3-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                     sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "size3"))


eigen_chr2 <- eigen[eigen$chrom == "chr2",]
slice1_chr2 <- slice1[slice1$chrom == "chr2",]
slice2_chr2 <- slice2[slice2$chrom == "chr2",]
slice3_chr2 <- slice3[slice3$chrom == "chr2",]

d <- data.frame(eigen_chr2[1:4], slice1_chr2["size1"], slice2_chr2["size2"], slice3_chr2["size3"])

d$size1[is.na(d$size1)] <- 0
d$size2[is.na(d$size2)] <- 0
d$size3[is.na(d$size3)] <- 0

pdf("slice_tracks.pdf", height = 4, width = 10)
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
plot(d$start/1000000, d$size1, type="n", 
     ylim=c(0, 3e-5), xlab = "", ylab = "", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$size1, rep(0, length(d$size1))), col="gray60",
        border=NA)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$size2, type="n",
     ylim=c(0, 3e-5), xlab = "", ylab= "", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$size2, rep(0, length(d$size2))), col="gray60",
        border=NA)
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
plot(d$start/1000000, d$size3, type="n",
     ylim=c(0, 3e-5), xlab = "", ylab="", axes=FALSE)
polygon(c(d$start/1000000, rev(d$start/1000000)), 
        c(d$size3, rep(0, length(d$size3))), col="gray60",
        border=NA)
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)
dev.off()








