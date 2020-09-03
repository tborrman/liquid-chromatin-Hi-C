library(pheatmap)
library(RColorBrewer)


df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

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

features <- c("half-life_LOS", "H3K36me3_R1", "H3K27ac_R1", "H3K4me1_R1",
              "H3K4me2_R1", "H3K4me3_R1", "H4K20me1_R1", "H2AFZ_R1", 
              "H3K9me1_R1", "H3K9me3_R1", "H3K9ac_R1", "H3K27me3_R1",
              "H3K79me2_R1", "LAD_K562", "CBX1_R1", "CBX3_Myers",
              "CBX5", "EHMT2", "CBX8", "RNF2", "BMI1", "SUZ12", "RBBP5",
              "CTBP1", "KAT2B", "BRD4", "NCOR1", "KDM5B", "HDAC2_Snyder",
              "SAP30", "WHSC1", "PHF8", "REST", "KDM1A_Snyder1", "POLR2A",
              "POLR2B", "POLR2G", "WGBS_R1", "NADs_IMR90", "SON_TSA-seq",
              "pSC35_TSA-seq", "LaminAC_TSA-seq", "LaminB_TSA-seq", 
              "Pol2_TSA-seq", "PML_R1")
              
    
new_labels <- c("halflife_LOS","H3K36me3", "H3K27ac", "H3K4me1",
                "H3K4me2", "H3K4me3", "H4K20me1", "H2AFZ", 
                "H3K9me1", "H3K9me3", "H3K9ac", "H3K27me3",
                "H3K79me2", "LADs.DamID", "CBX1", "CBX3",
                "CBX5", "EHMT2", "CBX8", "RNF2", "BMI1", "SUZ12", "RBBP5",
                "CTBP1", "KAT2B", "BRD4", "NCOR1", "KDM5B", "HDAC2",
                "SAP30", "WHSC1", "PHF8", "REST", "KDM1A", "POLR2A",
                "POLR2B", "POLR2G", "Methylation.CpG", "NADS.IMR90", "SON.TSA-seq",
                "pSC35.TSA-seq", "LaminAC.TSA-seq", "LaminB.TSA-seq", 
                "Pol2.TSA-seq", "PML")


dff <- df[features]
colnames(dff) <- new_labels

# Remove half-life NA rows
dff <- dff[!is.na(dff$halflife_LOS), ]

# Remove outliers
df_no_out <- data.frame(apply(dff, 2, remove_outliers_std))
# Remove half-life NA rows from outlier detection
df_no_out <- df_no_out[!is.na(df_no_out$halflife_LOS), ]
# Correlation for ordering
c <- as.data.frame(cor(df_no_out, use="pairwise.complete.obs", method="spearman"))
cordf <- data.frame(c)
cor_hlLOS <- cordf$halflife_LOS[2:length(cordf$halflife_LOS)]

hl_heatmap <- data.frame()
# Make zscore df
z_df <- data.frame(apply(df_no_out, 2, z_score), check.names=FALSE)

segments <- seq(40,120,5)

for (i in 1:(length(segments)-1)) {
  seg_df <- data.frame(z_df[df_no_out$halflife_LOS >= segments[i] & df_no_out$halflife_LOS < segments[i+1],2:ncol(df_no_out)], check.names=FALSE)
  order_df <- seg_df[,order(cor_hlLOS)]
  r <- apply(order_df, 2, mean, na.rm=TRUE)
  print(r)
  hl_heatmap <- rbind(hl_heatmap,r)
}
colnames(hl_heatmap) <- colnames(order_df)

t_hl_heatmap <- t(hl_heatmap)
colnames(t_hl_heatmap) <- c("40-45",
                            "45-50", "50-55", "55-60",
                            "60-65", "65-70", "70-75", 
                            "75-80", "80-85", "85-90", 
                            "90-95", "95-100", "100-105", 
                            "105-110", "110-115", "115-120")

t_hl_heatmap[t_hl_heatmap > 1.5] = 1.5
t_hl_heatmap[t_hl_heatmap < -1.25] = -1.25

pdf("feature_halflife_heatmap_supp.pdf", width=8, height=11, onefile=FALSE)
pheatmap(t_hl_heatmap,color=rev(colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(100)), 
         cluster_cols=FALSE)
dev.off()
