library(corrplot)
library(ggplot2)
df <- read.table("feature_matrix.txt", sep="\t", header=TRUE)
df_chrom <- df[df$chrom == "chr2",]

# plot correlation between std and LOS half-lifes
png('std_halflife_vs_LOS_half_life_scatter.png', height=2000, width=2000, res=300)
par(mar=c(5,5,4,2) + 0.1)
plot(df_chrom$half.life_LOS, df_chrom$half.life_std, 
     xlim=c(90,300), ylim=c(20,110), pch=21, col="black",
     bg="dodgerblue", xlab=bquote("t"[1/2] ~ " Loss of Structure (minutes)"),
     ylab="t"[1/2] ~ " "~ Delta~sigma ~ " (minutes)", cex.lab=1.5,
     main = "Chromosome 2", cex.main=1.5)
text(250,35, 
     bquote("r"[s]~ "="~  .(round(cor(df_chrom$half.life_LOS,df_chrom$half.life_std,
     method="spearman", use="complete.obs"),2))), cex=1.5)
dev.off()




features <- df_chrom[,4:ncol(df_chrom)]
cor_mat <- cor(features, use="pairwise.complete.obs", method="spearman")
png('correlate_feature_matrix_chr2.png',width=2500, height=2500, res=300)
  corrplot(cor_mat, method="circle",type = "upper", tl.col = "black")
dev.off()

png('correlate_feature_matrix_values_chr2.png',width=2500, height=2500, res=300)
corrplot(cor_mat, method="number", type="upper", tl.col = "black", number.cex = 0.5)
dev.off()

# holy moses this code is ugly, clean this shit up later

cor_values <- c(cor_mat["half.life_LOS","H3K4me1_R1"],cor_mat["half.life_std","H3K4me1_R1"],
                cor_mat["half.life_LOS","H4K20me1_R1"],cor_mat["half.life_std","H4K20me1_R1"],
                cor_mat["half.life_LOS","H2AFZ_R1"],cor_mat["half.life_std","H2AFZ_R1"],
                cor_mat["half.life_LOS","CTCF_R1"],cor_mat["half.life_std","CTCF_R1"],
                cor_mat["half.life_LOS","DNase.seq_R1"],cor_mat["half.life_std","DNase.seq_R1"],
                cor_mat["half.life_LOS","PCA_eigen1"],cor_mat["half.life_std","PCA_eigen1"],
                cor_mat["half.life_LOS","H3K27ac_R1"],cor_mat["half.life_std","H3K27ac_R1"],
                cor_mat["half.life_LOS","H3K36me3_R1"],cor_mat["half.life_std","H3K36me3_R1"],
                cor_mat["half.life_LOS","loops_Rao"],cor_mat["half.life_std","loops_Rao"],
                cor_mat["half.life_LOS","H3K27me3_R1"],cor_mat["half.life_std","H3K27me3_R1"],
                cor_mat["half.life_LOS","DpnII.seq"],cor_mat["half.life_std","DpnII.seq"],
                cor_mat["half.life_LOS","H3K9me3_R1"],cor_mat["half.life_std","H3K9me3_R1"],
                cor_mat["half.life_LOS","LAD_clone14"],cor_mat["half.life_std","LAD_clone14"])

feature_labs <- factor(c("H3K4me1_R1","H3K4me1_R1","H4K20me1_R1","H4K20me1_R1",
                          "H2AFZ_R1","H2AFZ_R1","CTCF_R1","CTCF_R1",
                          "DNase.seq_R1","DNase.seq_R1","PCA_eigen1","PCA_eigen1",
                          "H3K27ac_R1","H3K27ac_R1","H3K36me3_R1","H3K36me3_R1",
                          "loops_Rao","loops_Rao","H3K27me3_R1","H3K27me3_R1",
                          "DpnII.seq","DpnII.seq","H3K9me3_R1","H3K9me3_R1",
                          "LAD_clone14","LAD_clone14"),levels=c("PCA_eigen1","CTCF_R1",
                          "H2AFZ_R1", "H3K4me1_R1", "H3K27ac_R1","H4K20me1_R1",
                          "DNase.seq_R1","loops_Rao","H3K36me3_R1","H3K27me3_R1",
                          "DpnII.seq","H3K9me3_R1","LAD_clone14"), ordered=TRUE)

cor_df <- data.frame(cor_values, feature_labs)
cor_df$metric <- rep(c("LOS", "Delta sigma"), 13)

png("compare_delta_sigma_vs_los.png", height=2500, width=4500, res=600)
ggplot(cor_df, aes(x=feature_labs, y= cor_values,fill=metric)) +
        geom_bar(stat = "identity", position="dodge") +
        scale_fill_manual(values=rep(c("darkorange1", "chartreuse4"), 13)) +
        coord_flip() + labs(x = "", y = "Spearman's Correlation") +
        theme_minimal()
  
        
dev.off()


