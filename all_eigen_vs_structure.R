
wkdir <- getwd()
struct_table <- read.table("structure_500kb_sorted.bed", sep="\t", header=TRUE)
eigen_table <- read.table(paste(wkdir, "HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted.bedGraph", sep="/"), sep="\t")

# Fix bp indexing
eigen_table[2] <- eigen_table[2] - 1

colnames(eigen_table) <- c("chrom", "start", "end", "eigen")
colnames(struct_table) <- c("chrom", "start", "end", "struct")

# Merge
full_table <- merge(eigen_table, struct_table, c("chrom","start"))


png("all_eigen_vs_structure.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(full_table$eigen, full_table$struct, xlab= "PC1", ylab = "Loss of structure",
     col=rgb(0,0,0.7,0.3), pch=20, cex.lab=1.5)
#text(0.4, 0, paste("r = ",sprintf("%.2f",round(cor(full_table$eigen, full_table$struct),2)), sep=""), cex=1.5)
dev.off()

png("all_eigen_vs_structure_zoom.png", width=2500, height=2500, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(full_table$eigen, full_table$struct, xlab= "PC1", ylab = "Loss of structure", 
     xlim=c(-0.25, 0.25), ylim=c(0.1,1), col=rgb(0,0,0.7,0.3), pch=20, cex.lab=1.5)
#text(0.4, 0, paste("r = ",sprintf("%.2f",round(cor(full_table$eigen, full_table$struct),2)), sep=""), cex=1.5)
dev.off()


for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  chrom_df <- full_table[full_table$chrom == chrom,]
  
  png(paste("all_eigen_vs_structure_zoom_", chrom,".png", sep=""), width=2500, height=2500, res=300)
  par(mar=c(5,5,4,2) + 0.1)
  plot(chrom_df$eigen, chrom_df$struct, xlab= "PC1", ylab = "Loss of structure", 
       xlim=c(-0.25, 0.25), ylim=c(0.1,1), col=rgb(0,0,0.7,0.3), pch=20, cex.lab=1.5, 
       cex.main=1.5, main= chrom)
  #text(0.4, 0, paste("r = ",sprintf("%.2f",round(cor(full_table$eigen, full_table$struct),2)), sep=""), cex=1.5)
  dev.off()
}

# Nice plot for chr1
png(paste("all_eigen_vs_structure_zoom_nice_chr1.png", sep=""), width=2500, height=2500, res=300)
chrom_df <- full_table[full_table$chrom == "chr1",]
par(mar=c(5,5,4,2) + 0.1)
plot(chrom_df$eigen, chrom_df$struct, xlab= "PC1", ylab = "Loss of structure", 
      ylim=c(0.15,0.57), col="black", bg="blue", pch=21, cex.lab=1.5, 
      cex.main=2, main= "Chromosome 1")

dev.off()

# Nice plot for chr1 with corrleation
png(paste("all_eigen_vs_structure_zoom_nice_corr_chr1.png", sep=""), width=2500, height=2500, res=300)
chrom_df <- full_table[full_table$chrom == "chr1",]
par(mar=c(5,5,4,2) + 0.1)
plot(chrom_df$eigen, chrom_df$struct, xlab= "PC1", ylab = "Loss of structure", 
     ylim=c(0.15,0.57), col="black", bg="blue", pch=21, cex.lab=1.5, 
     cex.main=2, main= "Chromosome 1")
text(-0.05, 0.55, paste("r = ", round(cor(chrom_df$eigen, chrom_df$struct),2), sep=""), cex=1.5)
dev.off()

# Nice plot for chr1 with corrleation also A and B colored
png(paste("all_eigen_vs_structure_zoom_nice_corr_compColor_chr1_rect2.png", sep=""), width=2200, height=2800, res=300)
chrom_df <- full_table[full_table$chrom == "chr1",]
chrom_df_B <- chrom_df[chrom_df$eigen <= 0,]
chrom_df_A <- chrom_df[chrom_df$eigen > 0, ]

par(mar=c(5,5,4,2) + 0.1)
plot(chrom_df_B$eigen, chrom_df_B$struct, xlab= "PC1", ylab = "Loss of structure", 
     xlim = c(-0.07, 0.065), ylim=c(0.15,0.57), col="black", bg="blue", pch=21, cex.lab=1.5, 
     cex.main=2, main= "Chromosome 1")
points(chrom_df_A$eigen, chrom_df_A$struct, 
        col="black", bg="red", pch=21)
text(-0.05, 0.55, paste("r = ", round(cor(chrom_df$eigen, chrom_df$struct),2), sep=""), cex=1.5)
dev.off()
