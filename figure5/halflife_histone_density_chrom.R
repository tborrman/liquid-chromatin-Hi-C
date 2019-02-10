library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/figure4/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  print(chrom)
  df_chrom <- df[df$chrom == chrom,]
  h_na <- df_chrom[c("half-life_LOS", "H3K9me1_R1", "H3K9ac_R1", "H4K20me1_R1", "H3K4me1_R1",
                       "H3K27ac_R1", "H2AFZ_R1", "H3K27me3_R1", "H3K36me3_R1", "H3K9me3_R1")]
  colnames(h_na) <- c("hl", "H3K9me1_R1", "H3K9ac_R1", "H4K20me1_R1", "H3K4me1_R1",
                          "H3K27ac_R1", "H2AFZ_R1", "H3K27me3_R1", "H3K36me3_R1", "H3K9me3_R1")
  h <- na.omit(h_na)
  
  # Get upper quartile data 
  upq_H3K9me1_R1 <- h[h$H3K9me1_R1 > quantile(h$H3K9me1_R1)[4],]
  upq_H3K9ac_R1 <- h[h$H3K9ac_R1 > quantile(h$H3K9ac_R1)[4],]
  upq_H4K20me1_R1 <- h[h$H4K20me1_R1 > quantile(h$H4K20me1_R1)[4],]
  upq_H3K4me1_R1 <- h[h$H3K4me1_R1 > quantile(h$H3K4me1_R1)[4],]
  upq_H3K27ac_R1 <- h[h$H3K27ac_R1 > quantile(h$H3K27ac_R1)[4],]
  upq_H2AFZ_R1 <- h[h$H2AFZ_R1 > quantile(h$H2AFZ_R1)[4],]
  upq_H3K27me3_R1 <- h[h$H3K27me3_R1 > quantile(h$H3K27me3_R1)[4],]
  upq_H3K36me3_R1 <- h[h$H3K36me3_R1 > quantile(h$H3K36me3_R1)[4],]
  upq_H3K9me3_R1 <- h[h$H3K9me3_R1 > quantile(h$H3K9me3_R1)[4],]
  
  pdf(paste("chrom_histone/", chrom, "_half-life_histone_density.pdf", sep=""),
            width=7, height=6)
        print(
        ggplot() + 
          #stat_density(data=upq_H3K9me1_R1, aes(x=hl, color="H3K9me1"), geom="line") +
          stat_density(data=upq_H3K9ac_R1, aes(x=hl, color="H3K9ac"), geom="line") +
          #stat_density(data=upq_H4K20me1_R1, aes(x=hl, color="H4K20me1"), geom="line") +
          stat_density(data=upq_H3K4me1_R1 , aes(x=hl, color="H3K4me1"), geom="line") +
          #stat_density(data=upq_H3K27ac_R1 , aes(x=hl, color="H3K27ac"), geom="line") +
          #stat_density(data=upq_H2AFZ_R1 , aes(x=hl, color="H2AFZ"), geom="line") +
          #stat_density(data=upq_H3K27me3_R1 , aes(x=hl, color="H3K27me3"), geom="line") +
          stat_density(data=upq_H3K36me3_R1 , aes(x=hl, color="H3K36me3"), geom="line") +
          stat_density(data=upq_H3K9me3_R1 , aes(x=hl, color="H3K9me3"), geom="line") +
          xlim(0,400) + ylim(0,0.020) +
          xlab(bquote("t"[1/2] ~ "(minutes)")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.key = element_rect(fill="white"), legend.title = element_blank()) +
          scale_color_discrete(breaks=c("H3K9ac", "H3K4me1", "H3K36me3", "H3K9me3"))
        )
  dev.off()
  
  pdf(paste("chrom_histone/", chrom, "_half-life_histone_density_cumulative.pdf", sep=""),
      width=7, height=5)
  print(
    ggplot() + 
      #stat_density(data=upq_H3K9me1_R1, aes(x=hl, color="H3K9me1"), geom="line") +
      stat_ecdf(data=upq_H3K9ac_R1, aes(x=hl, color="H3K9ac"), geom="line") +
      #stat_density(data=upq_H4K20me1_R1, aes(x=hl, color="H4K20me1"), geom="line") +
      stat_ecdf(data=upq_H3K4me1_R1 , aes(x=hl, color="H3K4me1"), geom="line") +
      #stat_density(data=upq_H3K27ac_R1 , aes(x=hl, color="H3K27ac"), geom="line") +
      #stat_density(data=upq_H2AFZ_R1 , aes(x=hl, color="H2AFZ"), geom="line") +
      #stat_density(data=upq_H3K27me3_R1 , aes(x=hl, color="H3K27me3"), geom="line") +
      stat_ecdf(data=upq_H3K36me3_R1 , aes(x=hl, color="H3K36me3"), geom="line") +
      stat_ecdf(data=upq_H3K9me3_R1 , aes(x=hl, color="H3K9me3"), geom="line") +
      xlim(50,350) + 
      xlab(bquote("t"[1/2] ~ "(minutes)")) +
      ylab("F(x)") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill="white"), legend.title = element_blank()) +
      scale_color_discrete(breaks=c("H3K9ac", "H3K4me1", "H3K36me3", "H3K9me3"))
  )
  dev.off()
}
