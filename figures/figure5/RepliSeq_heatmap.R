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

features <- c("half-life_LOS", "G1_Repli-seq", "S1_Repli-seq", "S2_Repli-seq", 
              "S3_Repli-seq","S4_Repli-seq", "G2_Repli-seq")

new_labels <- c("halflife_LOS","G1 Repli-seq", "S1 Repli-seq", "S2 Repli-seq", 
                "S3 Repli-seq","S4 Repli-seq", "G2 Repli-seq")

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

pdf("RepliSeq_halflife_heatmap.pdf", width=8, height=2.5, onefile=FALSE)
pheatmap(t_hl_heatmap,color=rev(colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(100)), 
         cluster_cols=FALSE, cluster_rows=FALSE)
dev.off()

# Histogram
pdf(paste("half-life_hist_genome.pdf", sep=""), width=8, height=3)
hist(df_no_out$halflife_LOS, breaks=seq(20, 150, 5), xlim=c(40, 120), xlab="",
     main="", col="gray", probability=TRUE)
dev.off()

