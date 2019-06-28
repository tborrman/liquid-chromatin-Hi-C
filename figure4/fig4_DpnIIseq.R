library(pheatmap)
library(RColorBrewer)


DpnII_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN5mR1_S1_L001_copy_correct_coverage_40kb.bed",
                       sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "min5"))

graycol <-  colorRampPalette(brewer.pal(n = 9, name ="Greys")[3:9])(11)

pdf("fig4DpnIIseqSupp.pdf", height = 7, width = 6)
par(mfrow=c(11,1), mar=c(2, 4, 0, 2) + 0.1)
f <- list(c("HBDpSeqK562-DN5mR1_S1_L001_copy_correct_coverage_40kb.bed", "min5"),
  c("HBDpSeqK562-DN15mR1_S2_L001_copy_correct_coverage_40kb.bed", "min15"),
  c("HBDpSeqK562-DN30mR1_S3_L001_copy_correct_coverage_40kb.bed", "min30"),
  c("HBDpSeqK562-DN45mR1_S4_L001_copy_correct_coverage_40kb.bed", "min45"),
  c("HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed", "min60"),
  c("HBDpSeqK562-DN75mR1_S9_L003_copy_correct_coverage_40kb.bed", "min75"),
  c("HBDpSeqK562-DN90mR1_S10_L003_copy_correct_coverage_40kb.bed", "min90"),
  c("HBDpSeqK562-DN2hR1_S6_L002_copy_correct_coverage_40kb.bed", "min120"),
  c("HBDpSeqK562-DN3hR1_S7_L002_copy_correct_coverage_40kb.bed", "min180"),
  c("HBDpSeqK562-DN4hR1_S8_L002_copy_correct_coverage_40kb.bed", "min240"),
  c("HBDpSeqK562-DN16hR1_S11_L003_copy_correct_coverage_40kb.bed", "min960"))
for (i in 1:11) {
  d <- read.table(paste("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/",f[[i]][1], sep=""),
                  sep="\t", header=FALSE, col.names=c("chrom", "start", "end", f[[i]][2]))
  d_chr2 <- d[d$chrom == "chr2",]
  
  if (i == 11) {
    plot(d_chr2$start/1000000, d_chr2[[4]], type="l", col=graycol[i], 
     xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.5,
     ylim=c(0,2500))
    axis(1, lwd=2, cex.axis=1) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
  }
  else {
    plot(d_chr2$start/1000000, d_chr2[[4]], type="l", col=graycol[i], 
         xlab = "", ylab = paste(f[[i]][2], sep=" "), axes=FALSE, lwd=0.5,
         ylim=c(0,2500))
    axis(1, lwd=2, cex.axis=1, labels=FALSE) 
    axis(2, lwd=2, cex.axis=1)
    box(bty="l", lwd=2)
    
  }
   
   if (i > 1) {
     DpnII_df <- cbind(DpnII_df, d[f[[i]][2]])
   }
}
dev.off()

pdf("fig4DpnIIseqTracks.pdf", height = 4, width = 6)
par(mfrow=c(3,1), mar=c(2, 4, 0, 2) + 0.1)
min5 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN5mR1_S1_L001_copy_correct_coverage_40kb.bed",
                sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))
h1 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                   sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))
h2 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN2hR1_S6_L002_copy_correct_coverage_40kb.bed",
                 sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))

min5_chr2 <- min5[min5$chrom == "chr2",]
h1_chr2 <- h1[h1$chrom == "chr2",]
h2_chr2 <- h2[h2$chrom == "chr2",]

plot(min5_chr2$start/1000000, min5_chr2$DpnIIseq, type="l", col="orange", 
     xlab = "", ylab = "5 min DpnII-seq", axes=FALSE, lwd=0.5,
     ylim=c(0,3000))
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)

plot(h1_chr2$start/1000000, h1_chr2$DpnIIseq, type="l", col="springgreen4", 
     xlab = "", ylab = "1 hour DpnII-seq", axes=FALSE, lwd=0.5,
     ylim=c(0,3000))
axis(1, lwd=2, cex.axis=1, labels=FALSE) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)

plot(h2_chr2$start/1000000, h2_chr2$DpnIIseq, type="l", col="pink1", 
     xlab = "Chr2 position (Mb)", ylab = "2 hour DpnII-seq", axes=FALSE, lwd=0.5,
     ylim=c(0,3000))
axis(1, lwd=2, cex.axis=1) 
axis(2, lwd=2, cex.axis=1)
box(bty="l", lwd=2)

dev.off()

# Correlation plot
pdf('fig4DpnIIcorrMatrix.pdf',width=10, height=9.5, onefile=FALSE)
c <- cor(DpnII_df[4:ncol(DpnII_df)], use="pairwise.complete.obs", method="spearman")
pheatmap(c, color=(colorRampPalette((brewer.pal(n = 7, name ="YlOrRd")))(100)),
         breaks=seq(0.85, 1, .15/100), cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()

s <- DpnII_df[c("min5", "min60", "min120", "min180", "min240", "min960")]
pdf('fig4DpnIIcorrMatrix_subset.pdf',width=10, height=9.5, onefile=FALSE)
c <- cor(s, use="pairwise.complete.obs", method="spearman")
pheatmap(c, color=(colorRampPalette((brewer.pal(n = 7, name ="YlOrRd")))(100)),
         breaks=seq(0.85, 1, .15/100), cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


