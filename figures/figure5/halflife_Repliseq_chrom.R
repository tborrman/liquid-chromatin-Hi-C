library(ggplot2)
library(RColorBrewer)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("chrom", "half-life_LOS", "G1_Repli-seq", "S1_Repli-seq", "S2_Repli-seq", "S3_Repli-seq",
               "S4_Repli-seq", "G2_Repli-seq")]
colnames(hp1_na) <- c("chrom", "hl", "G1", "S1", "S2", "S3", "S4", "G2")
hp1 <- na.omit(hp1_na)

for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {
  hp1_chrom <- hp1[hp1$chrom == chrom,]

  # # Get upper quartile data 
  upq_G1 <- hp1_chrom[hp1_chrom$G1 > quantile(hp1_chrom$G1)[4],]
  upq_S1 <- hp1_chrom[hp1_chrom$S1 > quantile(hp1_chrom$S1)[4],]
  upq_S2 <- hp1_chrom[hp1_chrom$S2 > quantile(hp1_chrom$S2)[4],]
  upq_S3 <- hp1_chrom[hp1_chrom$S3 > quantile(hp1_chrom$S3)[4],]
  upq_S4 <- hp1_chrom[hp1_chrom$S4 > quantile(hp1_chrom$S4)[4],]
  upq_G2 <- hp1_chrom[hp1_chrom$G2 > quantile(hp1_chrom$G2)[4],]
  
  
  wave_fun <-function(x, b, n) {
    # Scale density curve
    # x = values
    # b = bin size
    # n = number of values
    d <- density(x)
    d$y <- d$y*b*n
    ddf <- data.frame(d$x, d$y)
    colnames(ddf) <- c("x", "Frequency")
    return(ddf)
  }
  
  #######################################################################################################################
  pdf(paste("Repli-seq/chrom/", chrom, "_half-life_Repli-seq_density_genome.pdf", sep=""),
            width=7, height=6)
        print(
        ggplot() + 
          geom_line(data=wave_fun(upq_G1$hl, 1, length(upq_G1$hl)), aes(x=x, y=Frequency, color="G1")) +
          geom_line(data=wave_fun(upq_S1$hl, 1, length(upq_S1$hl)), aes(x=x, y=Frequency, color="S1")) +
          geom_line(data=wave_fun(upq_S2$hl, 1, length(upq_S2$hl)), aes(x=x, y=Frequency, color="S2")) +
          geom_line(data=wave_fun(upq_S3$hl, 1, length(upq_S3$hl)), aes(x=x, y=Frequency, color="S3")) +
          geom_line(data=wave_fun(upq_S4$hl, 1, length(upq_S4$hl)), aes(x=x, y=Frequency, color="S4")) +
          geom_line(data=wave_fun(upq_G2$hl, 1, length(upq_G2$hl)), aes(x=x, y=Frequency, color="G2")) +
          xlim(30,150) +
          xlab(bquote("t"[1/2] ~ "(minutes)")) +
          theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
                panel.background = element_blank(), axis.line = element_line(colour = "black"),
                legend.key = element_rect(fill="white"), legend.title = element_blank()) +
          scale_color_manual(values=c("G1" = "#4575B4",
                                      "S1" = "#91BFDB",
                                      "S2" = "cyan",
                                      "S3" = "#FEE090",
                                      "S4" = "#FC8D59",
                                      "G2" = "#D73027"),
                             breaks=c("G1", "S1", "S2", "S3", "S4", "G2"))
        )
  dev.off()
  
  pdf(paste("Repli-seq/chrom/", chrom, "_half-life_Repli-seq_cumulative_genome.pdf", sep=""),
      width=7, height=6)
  print(
    ggplot() + 
      stat_ecdf(data=upq_G1, aes(x=hl, color="G1"), geom="line") +
      stat_ecdf(data=upq_S1, aes(x=hl, color="S1"), geom="line") +
      stat_ecdf(data=upq_S2, aes(x=hl, color="S2"), geom="line") +
      stat_ecdf(data=upq_S3, aes(x=hl, color="S3"), geom="line") +
      stat_ecdf(data=upq_S4, aes(x=hl, color="S4"), geom="line") +
      stat_ecdf(data=upq_G2, aes(x=hl, color="G2"), geom="line") +
      xlim(30,140) +
      xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab(bquote("F"["n"] ~ "(x)")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("G1" = "#4575B4",
                                "S1" = "#91BFDB",
                                "S2" = "cyan",
                                "S3" = "#FEE090",
                                "S4" = "#FC8D59",
                                "G2" = "#D73027"),
                       breaks=c("G1", "S1", "S2", "S3", "S4", "G2"))
  )
  dev.off()

}
      

