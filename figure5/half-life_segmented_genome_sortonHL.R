library(pheatmap)
library(RColorBrewer)



df <- read.table("/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt", sep="\t",
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
              "H3K27me3_R1", "DNase-seq_R1", "LAD_K562",
              "CBX1_R1", "CBX3_Myers", "CBX5", "SUZ12", "HDAC2_Snyder", 
              "WGBS_R1", "NADs_IMR90", "gene_density", "SON_TSA-seq", "PML_R1",
              "H3K4me3_R1", "CTCF_R1", "SMC3_R1", "TAD_insulation")

new_labels <- c("halflife_LOS", "H3K36me3", "H3K27ac", "H3K4me1", "H4K20me1",
                "H2AFZ", "H3K9me1", "H3K9me3", "H3K9ac", 
                "H3K27me3", "DNase-seq", "LADs",
                "CBX1", "CBX3", "CBX5", "SUZ12", "HDAC2", 
                "WGBS", "NADs", "gene density", "TSA-seq", "PML",
                "H3K4me3", "CTCF", "SMC3", "Insulation")


dff <- df[features]
colnames(dff) <- new_labels

# Remove half-life NA rows
dff <- dff[!is.na(dff$halflife_LOS), ]
#test <- df[!is.na(df$`half-life_LOS`),]

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
# Order by halflife
myorder_df <- data.frame()
for (i in 1:4) {
  myseg_df <- data.frame(z_df[df_no_out$halflife_LOS >= segments[i] & df_no_out$halflife_LOS < segments[i+1],2:ncol(df_no_out)], check.names=FALSE)
  mym <- apply(myseg_df, 2, mean, na.rm=TRUE)
  myorder_df <- rbind(myorder_df, mym)
}
colnames(myorder_df) <- names(mym)
mym <- apply(myorder_df, 2, mean, na.rm=TRUE)
myorder <- rev(order(mym))


for (i in 1:(length(segments)-1)) {
  seg_df <- data.frame(z_df[df_no_out$halflife_LOS >= segments[i] & df_no_out$halflife_LOS < segments[i+1],2:ncol(df_no_out)], check.names=FALSE)
  #order_df <- seg_df[,order(cor_hlLOS)]
  r <- apply(seg_df, 2, mean, na.rm=TRUE)
  print(r)
  order_r <- r[myorder]
  hl_heatmap <- rbind(hl_heatmap,order_r)
}
colnames(hl_heatmap) <- names(order_r)

t_hl_heatmap <- t(hl_heatmap)
colnames(t_hl_heatmap) <- c("40-45",
                            "45-50", "50-55", "55-60",
                            "60-65", "65-70", "70-75",
                            "75-80", "80-85", "85-90",
                            "90-95", "95-100", "100-105",
                            "105-110", "110-115", "115-120")

# SATURATION
t_hl_heatmap[t_hl_heatmap > 1.5] = 1.5

pdf("half-life_segmented_genome_mean_saturated_sortonHL.pdf", width=8, height=8, onefile=FALSE)
pheatmap(t_hl_heatmap,color=colorRampPalette(rev(brewer.pal(n = 7, name ="RdYlBu")))(100), 
         cluster_cols=FALSE, cluster_rows=FALSE
         , breaks=seq(-1.25,1.5,2.75/100)
)
dev.off()

# Histogram
pdf(paste("half-life_hist_genome.pdf", sep=""), width=8, height=3)
hist(df_no_out$halflife_LOS, breaks=seq(5, 500, 5), xlim=c(40, 120), xlab="",
     main="", col="gray", freq=FALSE)
dev.off()



