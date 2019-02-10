library(corrplot)
library(ggplot2)
library(pheatmap)
library(RColorBrewer)

df <- read.table("/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt", sep="\t",
                 header=TRUE, check.names=FALSE)

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
              "WGBS_R1", "NADs_IMR90", "gene_density", "SON_TSA-seq", "PML_R1")

new_labels <- c("half-life_LOS", "H3K36me3", "H3K27ac", "H3K4me1", "H4K20me1",
                "H2AFZ", "H3K9me1", "H3K9me3", "H3K9ac", 
                "H3K27me3", "DNase-seq", "LADs",
                "CBX1", "CBX3", "CBX5", "SUZ12", "HDAC2",
                "WGBS", "NADs", "gene density", "TSA-seq", "PML")


  sub_df <- df[features]
  colnames(sub_df) <- new_labels
  
  # Remove half-life NA rows
  dff <- sub_df[!is.na(sub_df$`half-life_LOS`), ]
  
  # Remove outliers
  df_no_out <- data.frame(apply(dff, 2, remove_outliers_std))
  # Remove half-life NA rows from outlier detection
  df_no_out <- df_no_out[!is.na(df_no_out$half.life_LOS), ]
  # Correlation
  c <- as.data.frame(cor(df_no_out, use="pairwise.complete.obs", method="spearman"))
  cordf <- data.frame(c)
  cor_hlLOS <- cordf$half.life_LOS
  # Barplot
  bar_df <- data.frame(cor_hlLOS, new_labels)[2:length(cor_hlLOS),]
  order_bar <- bar_df[order(bar_df$cor_hlLOS, decreasing = TRUE),]
  
  ordered_labels <- factor(order_bar$new_labels, levels = order_bar$new_labels, ordered=TRUE)
  order_bar$new_labels <- ordered_labels
  
  mycolors <- ifelse(order_bar$cor_hlLOS < 0, "red", "dodgerblue")
  
  pdf("cor_los_halflife_40kb_genome.pdf", height=8, width=7.5)
  print(ggplot(order_bar, aes(x=new_labels, y=cor_hlLOS, fill=new_labels)) +
    geom_bar(stat = "identity") +
    scale_fill_manual(values=mycolors) +
    coord_flip() + labs(x = "", y = expression(rho)) +
    ylim(-1, 1) + 
    theme_bw() + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          axis.line.x = element_line(colour="black"), axis.ticks.y = element_blank(),
          legend.position="none")
    )
  dev.off()


