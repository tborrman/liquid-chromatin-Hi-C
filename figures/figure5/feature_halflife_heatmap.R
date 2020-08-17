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

features <- c("half-life_LOS", "SON_TSA-seq", "pSC35_TSA-seq", "Pol2_TSA-seq",
              "H3K9ac_R1", "H3K27ac_R1","H3K4me3_R1", "H3K36me3_R1", "SUZ12",
              "RNF2", "BMI1", "CBX8", "H3K27me3_R1", "CBX3_Myers", "CBX5", "CBX1_R1",
              "H3K9me3_R1", "NADs_IMR90", "LaminAC_TSA-seq", "LAD_K562", "LaminB_TSA-seq")

new_labels <- c("halflife_LOS", "SON", "pSC35", "Pol2",
                "H3K9ac", "H3K27ac","H3K4me3", "H3K36me3", "SUZ12",
                "RNF2", "BMI1", "CBX8", "H3K27me3", "CBX3", "CBX5", "CBX1",
                "H3K9me3", "NADs", "LaminA/C", "LADs", "LaminB")



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

z_df <- z_df[]

segments <- seq(40,120,5)

for (i in 1:(length(segments)-1)) {
  seg_df <- data.frame(z_df[df_no_out$halflife_LOS >= segments[i] & df_no_out$halflife_LOS < segments[i+1],2:ncol(df_no_out)], check.names=FALSE)
  #order_df <- seg_df[,order(cor_hlLOS)]
  order_df <- seg_df
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

pdf("feature_halflife_heatmap.pdf", width=8, height=7.1, onefile=FALSE)
pheatmap(t_hl_heatmap,color=rev(colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(100)), 
         cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()