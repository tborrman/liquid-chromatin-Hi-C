library(corrplot)
library(pheatmap)
library(RColorBrewer)
df <- read.table("../feature_matrix.txt", sep="\t", header=TRUE)
hl_df <- data.frame()

z_score <- function(x) {
  z <- (x - mean(x, na.rm=TRUE))/ (sd(x, na.rm=TRUE))
  return(z)
}

for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  chrom_df <- df[df$chrom == chrom,]
  features <- chrom_df[,4:ncol(chrom_df)]
  c <- as.data.frame(cor(features, use="pairwise.complete.obs"))
  hl_chr <-c[1,2:ncol(c)]
  rownames(hl_chr) <- chrom
  hl_df <- rbind(hl_df, hl_chr)
  
  # z-score chrom matrices
  z_chrom_df <- data.frame(chrom_df[1:3], apply(chrom_df[4:ncol(chrom_df)], 2, z_score))
  # Remove LADs
  z_chrom_df <- z_chrom_df[-c(24,25)]
  # Transpose
  z_chrom_df <- t(as.matrix(z_chrom_df[4:ncol(z_chrom_df)]))
  
  png(paste("z-score_features_", chrom, ".png",sep=""), width=6000, height=2500, res=300)
  pheatmap(z_chrom_df,color=rev(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)),
           cluster_cols=FALSE, labels_col = "")
  dev.off()
  # png(paste('correlate_feature_matrix_', chrom,'.png', sep=""),width=2500, height=2500, res=300)
  # corrplot(c, method="circle",type = "upper", tl.col = "black")
  # dev.off()
  # 
  # png(paste('correlate_feature_matrix_values_', chrom,'.png', sep=""),width=2500, height=2500, res=300)
  # corrplot(c, method="number", type="upper", tl.col = "black", number.cex = 0.5)
  # dev.off()
}

png("correlate_feature_matrix_half_life_chrom.png", width=2500, height=2500, res=300)
  pheatmap(hl_df,color=rev(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)))
dev.off()
  