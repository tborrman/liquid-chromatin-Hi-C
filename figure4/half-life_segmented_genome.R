library(pheatmap)
library(RColorBrewer)



df <- read.table("/cygwin64/home/Tyler/Research/digest/figure4/feature_matrix_v3_40kb.txt", sep="\t",
                 header=TRUE, check.names=FALSE)

z_score <- function(x) {
  z <- (x - mean(x, na.rm=TRUE))/ (sd(x, na.rm=TRUE))
  return(z)
}

remove_outliers_std <- function(x) {
  sigma <- sd(x, na.rm=TRUE)
  mu <- mean(x, na.rm=TRUE)
  lower_bound <- mu - (sigma*3)
  upper_bound <- mu + (sigma*3)
  x[(x > upper_bound) | (x < lower_bound)] <- NA
  return(x)
}

# Hand picked features for figure 4
features <- c("half-life_LOS", "H3K36me3_R1", "H3K27ac_R1", "H3K4me1_R1", "H4K20me1_R1",
              "H2AFZ_R1", "H3K9me1_R1", "H3K9me3_R1", "H3K9ac_R1", 
              "H3K27me3_R1", "DNase-seq_R1", "DpnII-seq", "LAD_clone14",
              "CBX1_R1", "CBX3_Myers", "CBX5", "SUZ12", 
              "HDAC2_Snyder","G1_Repli-seq", "S1_Repli-seq", "S2_Repli-seq", 
              "S3_Repli-seq","S4_Repli-seq", "G2_Repli-seq", "WGBS_R1", 
              "NADs_IMR90", "RNA-seq_total_+_R1", "gene_density")

new_labels <- c("halflife_LOS", "H3K36me3", "H3K27ac", "H3K4me1", "H4K20me1",
                "H2AFZ", "H3K9me1", "H3K9me3", "H3K9ac", 
                "H3K27me3", "DNase-seq", "DpnII-seq", "LADs",
                "CBX1", "CBX3", "CBX5", "SUZ12", 
                "HDAC2","G1 Repli-seq", "S1 Repli-seq", "S2 Repli-seq", 
                "S3 Repli-seq","S4 Repli-seq", "G2 Repli-seq", "WGBS", 
                "NADs", "RNA-seq", "gene density")



dff <- df[features]
colnames(dff) <- new_labels

# Remove half-life NA rows
dff <- dff[!is.na(dff$halflife_LOS), ]
#test <- df[!is.na(df$`half-life_LOS`),]

# Remove outliers
df_no_out <- data.frame(apply(dff, 2, remove_outliers_std))

hl_heatmap <- data.frame()


# Make zscore df
z_df <- data.frame(apply(df_no_out, 2, z_score), check.names=FALSE)


segments <- seq(50,275,15)

for (i in 1:(length(segments)-1)) {
  seg_df <- data.frame(z_df[df_no_out$halflife_LOS >= segments[i] & df_no_out$halflife_LOS < segments[i+1],2:ncol(dff)], check.names=FALSE)
  #test_seg <- data.frame(test[dff$halflife_LOS >= segments[i] & dff$halflife_LOS < segments[i+1],1:ncol(dff)])
  r <- apply(seg_df, 2, mean, na.rm=TRUE)
  print(r)
  hl_heatmap <- rbind(hl_heatmap,r)
}
colnames(hl_heatmap) <- colnames(df_no_out)[2:ncol(df_no_out)]

t_hl_heatmap <- t(hl_heatmap)
colnames(t_hl_heatmap) <- c("50-65",
                            "65-80", "80-95", "95-110",
                            "110-125", "125-140", "140-155", 
                            "155-170", "170-185", "185-200", 
                            "200-215", "215-230", "230-245", 
                            "245-260", "260-275")

png("half-life_segmented_genome_mean.png", width=3500, height=3500, res=300)
pheatmap(t_hl_heatmap,color=rev(colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100)), 
         cluster_cols=FALSE)
dev.off()


