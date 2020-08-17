eigen_df <- read.table("HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted.size.bedGraph", header=FALSE, sep="\t")
LOS_df <- read.table("half-life_exponential_500kb.bedGraph", header=FALSE, sep="\t")
colnames(eigen_df) <- c("chrom", "start", "end", "eigen", "size")
colnames(LOS_df) <- c("chrom", "start", "end", "LOS.halflife")


for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {

  chr_e <- eigen_df[eigen_df$chrom == chrom,]
  chr_LOS <- LOS_df[LOS_df$chrom == chrom,]
  
  options(scipen=100)
  png(paste("compartment_size_vs_half-life_LOS_500kb_", chrom,".png",sep=""), height=2500, width=2500, res=300)
  par(mar=c(5,6,4,2) + 0.1)
  plot(chr_e$size, chr_LOS$LOS.halflife, pch=20, col=ifelse(chr_e$eigen < 0, "blue", "red"),
       xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5, main= chrom)
  dev.off()

}

# Nice chrom 2
chr_e <- eigen_df[eigen_df$chrom == "chr2",]
chr_LOS <- LOS_df[LOS_df$chrom == "chr2",]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr2_nice.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(chr_e$size, chr_LOS$LOS.halflife, pch=20, col=ifelse(chr_e$eigen < 0, "blue", "red"),
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 2", ylim=c(95, 305), xlim=c(0,13000000))
dev.off()

A_e <- chr_e[chr_e$eigen >= 0,]
A_LOS <- chr_LOS[chr_e$eigen >=0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr2_A.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(A_e$size, A_LOS$LOS.halflife, pch=20, col="red",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 2", ylim=c(95, 305), xlim=c(0,13000000))

text(2000000, 275, 
     bquote("r"[s]~ "="~  .(round(cor(A_e$size, A_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()

B_e <- chr_e[chr_e$eigen < 0,]
B_LOS <- chr_LOS[chr_e$eigen < 0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr2_B.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(B_e$size, B_LOS$LOS.halflife, pch=20, col="blue",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 2", ylim=c(95, 305), xlim=c(0,13000000))

text(2000000, 275, 
     bquote("r"[s]~ "="~  .(round(cor(B_e$size, B_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()

# Nice chrom 6
chr_e <- eigen_df[eigen_df$chrom == "chr6",]
chr_LOS <- LOS_df[LOS_df$chrom == "chr6",]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr6_nice.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(chr_e$size, chr_LOS$LOS.halflife, pch=20, col=ifelse(chr_e$eigen < 0, "blue", "red"),
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 6", ylim=c(70,265), xlim=c(0,21500000))
dev.off()

A_e <- chr_e[chr_e$eigen >= 0,]
A_LOS <- chr_LOS[chr_e$eigen >=0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr6_A.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(A_e$size, A_LOS$LOS.halflife, pch=20, col="red",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 6", ylim=c(70,265), xlim=c(0,21500000))

text(5000000, 250, 
     bquote("r"[s]~ "="~  .(round(cor(A_e$size, A_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()

B_e <- chr_e[chr_e$eigen < 0,]
B_LOS <- chr_LOS[chr_e$eigen < 0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_chr6_B.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(B_e$size, B_LOS$LOS.halflife, pch=20, col="blue",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Chromosome 6", ylim=c(70,265), xlim=c(0,21500000))

text(5000000, 250, 
     bquote("r"[s]~ "="~  .(round(cor(B_e$size, B_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()

# Nice Genome
png(paste("compartment_size_vs_half-life_LOS_500kb_genome.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(eigen_df$size, LOS_df$LOS.halflife, pch=20, col=ifelse(eigen_df$eigen < 0, "blue", "red"),
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Genome", xlim=c(0,20000000), ylim=c(50,350))
dev.off()

A_e <- eigen_df[eigen_df$eigen >= 0,]
A_LOS <- LOS_df[eigen_df$eigen >=0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_genome_A.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(A_e$size, A_LOS$LOS.halflife, pch=20, col="red",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Genome", xlim=c(0,20000000), ylim=c(50,350))

text(5000000, 325, 
     bquote("r"[s]~ "="~  .(round(cor(A_e$size, A_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()

B_e <- eigen_df[eigen_df$eigen < 0,]
B_LOS <- LOS_df[eigen_df$eigen < 0,]

png(paste("compartment_size_vs_half-life_LOS_500kb_genome_B.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(B_e$size, B_LOS$LOS.halflife, pch=20, col="blue",
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Genome", xlim=c(0,20000000), ylim=c(50,350))

text(5000000, 325, 
     bquote("r"[s]~ "="~  .(round(cor(B_e$size, B_LOS$LOS.halflife, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()
