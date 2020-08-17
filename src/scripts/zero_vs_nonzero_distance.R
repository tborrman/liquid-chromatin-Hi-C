z_df <- read.table("dpnII_zero_distance_coverage_sorted_500kb_R1.bed", sep="\t", header=FALSE)
nz_df <- read.table("dpnII_nonzero_distance_coverage_sorted_500kb_R1.bed", sep="\t", header=FALSE)
z_df <- z_df[1:4]
nz_df <- nz_df[1:4]
colnames(z_df) <- c("chrom", "start", "end", "z_reads")
colnames(nz_df) <- c("chrom", "start", "end", "nz_reads")

# Merge the two files into one table based on the genomic position
full_table <- merge(z_df, nz_df, c("chrom","start", "end"))
options(scipen=100)
chroms <- c(paste(rep("chr", 22), 1:22, sep=""),"chrX")
for (i in seq(1,23)) {
  df_chr <- full_table[which(full_table$chrom == chroms[i]),]
  df_chr <- df_chr[order(df_chr$start),]
  df_gen <- df_chr$start + ((df_chr$end - df_chr$start)/ 2)
  png(paste("zero_vs_nonzero_chr", i, ".png",sep=""), width=5000, height=3000, res=400)
    par(mfrow=c(2,1))
    if (chroms[i] == 'chr8') {
      plot(df_gen, df_chr$z_reads, pch=20, col="dodgerblue", xlim=c(0,max(df_gen)), 
           xlab= "", ylab="DpnII reads (500kb)", type="o", ylim=c(300,900),
           main = "Zero-distance fragments")
      plot(df_gen, df_chr$nz_reads, pch=20, col="red", xlim=c(0,max(df_gen)), 
           xlab= chroms[i], ylab="DpnII reads (500kb)", type="o", ylim=c(900,2700),
           main = "Nonzero-distance fragments")
    }
    else {
    plot(df_gen, df_chr$z_reads, pch=20, col="dodgerblue", xlim=c(0,max(df_gen)), 
        xlab = "",ylab="DpnII reads (500kb)", type="o",
         main = "Zero-distance fragments")
    plot(df_gen, df_chr$nz_reads, pch=20, col="red", xlim=c(0,max(df_gen)), 
         xlab= chroms[i], ylab="DpnII reads (500kb)", type="o",
         main = "Nonzero-distance fragments")
  }
  dev.off()
}
# Full genome plot
png("zero_vs_nonzero_scatter.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(full_table$z_reads, full_table$nz_reads, 
     xlab="DpnII reads 500kb (Zero-distance)",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Genome", xlim=c(200,1250), ylim=c(600,3800))
text(1000,1000, bquote("r"[s]~ "="~  
     .(round(cor(full_table$z_reads, full_table$nz_reads, method="spearman"),2))), 
     cex=1.5)
dev.off()
# Chr 8
png("zero_vs_nonzero_scatter_chr8.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
z_reads_chr8 <- full_table$z_reads[full_table$chrom =="chr8"]
nz_reads_chr8 <- full_table$nz_reads[full_table$chrom =="chr8"]
plot(z_reads_chr8, nz_reads_chr8,
     xlab="DpnII reads 500kb (Zero-distance)",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
     #xlim=c(200,1250), ylim=c(600,3800))
text(600,500, bquote("r"[s]~ "="~  
                         .(round(cor(z_reads_chr8, nz_reads_chr8, method="spearman"),2))), 
     cex=1.5)
dev.off()
# Chr 8 zoom
png("zero_vs_nonzero_scatter_chr8_zoomxabove500yabove1600.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
table_chr8 <- full_table[full_table$chrom =="chr8",]
table_chr8 <- table_chr8[table_chr8$z_reads > 500 &  table_chr8$nz_reads > 1600,]
z_reads_chr8 <- table_chr8$z_reads
nz_reads_chr8 <- table_chr8$nz_reads

plot(z_reads_chr8, nz_reads_chr8,
     xlab="DpnII reads 500kb (Zero-distance)",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8 (x above 500, y above 1600)") 
#xlim=c(200,1250), ylim=c(600,3800))
text(800,1700, bquote("r"[s]~ "="~  
                       .(round(cor(z_reads_chr8, nz_reads_chr8, method="spearman"),2))), 
     cex=1.5)
dev.off()

# Correlation with compartments, dpnII sites, and LOS of structure
table_chr8 <- full_table[full_table$chrom =="chr8",]
eigen <- read.table("C:/cygwin64/home/Tyler/Research/digest/cis_percent/compartment/HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted.bedGraph", sep="\t", header=FALSE)
sites <- read.table("C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/dpnII_site_coverage_sorted.bed", sep="\t", header=FALSE)
sites <- sites[1:4]
LOS_hl <- read.table("C:/cygwin64/home/Tyler/Research/digest/cis_percent/timecourse/half-life_LOS/half-life_exponential_500kb.bedGraph", sep="\t", header=FALSE)
colnames(LOS_hl) <- c("chrom", "start", "end", "halflife")
colnames(sites) <- c("chrom", "start", "end", "sites")
colnames(eigen) <- c("chrom", "start", "end", "eigen")
eigen["start"] <- eigen["start"] -1
eigen_chr8 <- eigen[eigen$chrom =="chr8",]
sites_chr8 <- sites[sites$chrom == "chr8",]
hl_chr8 <- LOS_hl[LOS_hl$chrom == "chr8",]

eigen_dpnII_chr8 <- merge(table_chr8, eigen_chr8, c("chrom", "start", "end"))
eigen_dpnII_chr8 <- merge(eigen_dpnII_chr8, sites_chr8, c("chrom", "start", "end"))
eigen_dpnII_chr8 <- merge(eigen_dpnII_chr8, hl_chr8, c("chrom", "start", "end"))

png("zero_vs_LOS_halflife.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$halflife, eigen_dpnII_chr8$z_reads,
     xlab="LOS Half-life (minutes)",
     ylab="DpnII reads 500kb (Zero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(100, 200, bquote("r"[s]~ "="~  
                        .(round(cor(eigen_dpnII_chr8$halflife, eigen_dpnII_chr8$z_reads, method="spearman", 
                                    use="complete.obs"),2))), 
     cex=1.5)
dev.off()


png("nonzero_vs_LOS_halflife.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$halflife, eigen_dpnII_chr8$nz_reads,
     xlab="LOS Half-life (minutes)",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(100, 500, bquote("r"[s]~ "="~  
                        .(round(cor(eigen_dpnII_chr8$halflife, eigen_dpnII_chr8$nz_reads, method="spearman",
                                    use="complete.obs"),2))), 
     cex=1.5)
dev.off()



png("zero_vs_sites.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$sites, eigen_dpnII_chr8$z_reads,
     xlab="DpnII sites 500kb",
     ylab="DpnII reads 500kb (Zero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(250, 750, bquote("r"[s]~ "="~  
                         .(round(cor(eigen_dpnII_chr8$sites, eigen_dpnII_chr8$z_reads, method="spearman"),2))), 
     cex=1.5)
dev.off()

png("nonzero_vs_sites.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$sites, eigen_dpnII_chr8$nz_reads,
     xlab="DpnII sites 500kb",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(250, 2000, bquote("r"[s]~ "="~  
                        .(round(cor(eigen_dpnII_chr8$sites, eigen_dpnII_chr8$nz_reads, method="spearman"),2))), 
     cex=1.5)
dev.off()




png("zero_vs_eigen.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$eigen, eigen_dpnII_chr8$z_reads,
     xlab="PC1",
     ylab="DpnII reads 500kb (Zero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(-0.06,175, bquote("r"[s]~ "="~  
                        .(round(cor(eigen_dpnII_chr8$eigen, eigen_dpnII_chr8$z_reads, method="spearman"),2))), 
     cex=1.5)
dev.off()

png("nonzero_vs_eigen.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(eigen_dpnII_chr8$eigen, eigen_dpnII_chr8$nz_reads,
     xlab="PC1",
     ylab="DpnII reads 500kb (Nonzero-distance)",
     pch=21, col="black", bg="dodgerblue",
     main = "Chromosome 8") 
#xlim=c(200,1250), ylim=c(600,3800))
text(-0.06,175, bquote("r"[s]~ "="~  
                         .(round(cor(eigen_dpnII_chr8$eigen, eigen_dpnII_chr8$nz_reads, method="spearman"),2))), 
     cex=1.5)
dev.off()





