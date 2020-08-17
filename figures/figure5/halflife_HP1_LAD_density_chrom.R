library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/figure4/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  if (chrom == "chr8") {
    next
  }
  print(chrom)
  df_chrom <- df[df$chrom == chrom,]
  hp1_na <- df_chrom[c("half-life_LOS", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD_clone14")]
  colnames(hp1_na) <- c("hl", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD")
  hp1 <- na.omit(hp1_na)
  
  # Get upper quartile data 
  upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1)[4],]
  upq_CBX3 <- hp1[hp1$CBX3_Myers > quantile(hp1$CBX3_Myers)[4],]
  upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5)[4],]
  upq_LAD <- hp1[hp1$LAD > quantile(hp1$LAD)[4],]
  
  
  pdf(paste("chrom_HP1_LAD/", chrom, "_half-life_HP1_LAD_density.pdf", sep=""),
            width=7, height=6)
        print(
        ggplot() + 
          stat_density(data=upq_CBX1, aes(x=hl, color="CBX1"), geom="line") +
          stat_density(data=upq_CBX3, aes(x=hl, color="CBX3"), geom="line") +
          stat_density(data=upq_CBX5, aes(x=hl, color="CBX5"), geom="line") +
          stat_density(data=upq_LAD, aes(x=hl, color="LAD"), geom="line") +
          xlim(0,400) + ylim(0,0.020) +
          xlab(bquote("t"[1/2] ~ "(minutes)")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.key = element_rect(fill="white"), legend.title = element_blank())
        )
  dev.off()
  
  
  pdf(paste("chrom_HP1_LAD/", chrom, "_half-life_HP1_LAD_density_cumulative.pdf", sep=""),
      width=7, height=5)
  print(
    ggplot() + 
      stat_ecdf(data=upq_CBX1, aes(x=hl, color="CBX1"), geom="line") +
      stat_ecdf(data=upq_CBX3, aes(x=hl, color="CBX3"), geom="line") +
      stat_ecdf(data=upq_CBX5, aes(x=hl, color="CBX5"), geom="line") +
      stat_ecdf(data=upq_LAD, aes(x=hl, color="LAD"), geom="line") +
      xlim(50, 350) +
      xlab(bquote("t"[1/2] ~ "(minutes)")) +
      ylab("F(x)") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill="white"), legend.title = element_blank())
  )
  dev.off()
}