library(ggplot2)
library(RColorBrewer)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

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

plot_RNAseq <- function(d, e) {
  # Plot RNAseq wave and cumulative plots
  # e = compartment
  # d = dataframe of half-life and RNAseq for compartment
  
  # Get quantile  data 
  upq_Q1 <- d[d$RNAseq < quantile(d$RNAseq, 1/6),]
  upq_Q2 <- d[(d$RNAseq >= quantile(d$RNAseq, 1/6)) & (d$RNAseq < quantile(d$RNAseq, 2/6)),]
  upq_Q3 <- d[(d$RNAseq >= quantile(d$RNAseq, 2/6)) & (d$RNAseq < quantile(d$RNAseq, 3/6)),]
  upq_Q4 <- d[(d$RNAseq >= quantile(d$RNAseq, 3/6)) & (d$RNAseq < quantile(d$RNAseq, 4/6)),]
  upq_Q5 <- d[(d$RNAseq >= quantile(d$RNAseq, 4/6)) & (d$RNAseq < quantile(d$RNAseq, 5/6)),]
  upq_Q6 <- d[d$RNAseq >= quantile(d$RNAseq, 5/6),]
  
  #######################################################################################################################
  pdf(paste("RNA-seq/half-life_RNA-seq_density_genome_", e, ".pdf", sep=""),
      width=7, height=6)
  print(
    ggplot() + 
      geom_line(data=wave_fun(upq_Q1$hl, 1, length(upq_Q1$hl)), aes(x=x, y=Frequency, color="Q1")) +
      geom_line(data=wave_fun(upq_Q2$hl, 1, length(upq_Q2$hl)), aes(x=x, y=Frequency, color="Q2")) +
      geom_line(data=wave_fun(upq_Q3$hl, 1, length(upq_Q3$hl)), aes(x=x, y=Frequency, color="Q3")) +
      geom_line(data=wave_fun(upq_Q4$hl, 1, length(upq_Q4$hl)), aes(x=x, y=Frequency, color="Q4")) +
      geom_line(data=wave_fun(upq_Q5$hl, 1, length(upq_Q5$hl)), aes(x=x, y=Frequency, color="Q5")) +
      geom_line(data=wave_fun(upq_Q6$hl, 1, length(upq_Q6$hl)), aes(x=x, y=Frequency, color="Q6")) +
      xlim(30,100) +
      xlab(bquote("t"[1/2] ~ "(minutes)")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill="white"), legend.title = element_blank()) +
      scale_color_manual(values=c("Q1" = "#4575B4",
                                  "Q2" = "#91BFDB",
                                  "Q3" = "cyan",
                                  "Q4" = "#FEE090",
                                  "Q5" = "#FC8D59",
                                  "Q6" = "#D73027"),
                         breaks=c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"))
  )
  dev.off()
  
  pdf(paste("RNA-seq/half-life_RNA-seq_cumulative_genome_", e, ".pdf", sep=""),
      width=7, height=6)
  print(
    ggplot() + 
      stat_ecdf(data=upq_Q1, aes(x=hl, color="Q1"), geom="line") +
      stat_ecdf(data=upq_Q2, aes(x=hl, color="Q2"), geom="line") +
      stat_ecdf(data=upq_Q3, aes(x=hl, color="Q3"), geom="line") +
      stat_ecdf(data=upq_Q4, aes(x=hl, color="Q4"), geom="line") +
      stat_ecdf(data=upq_Q5, aes(x=hl, color="Q5"), geom="line") +
      stat_ecdf(data=upq_Q6, aes(x=hl, color="Q6"), geom="line") +
      xlim(30,100) +
      xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab(bquote("F"["n"] ~ "(x)")) +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black"),
            legend.key = element_rect(fill="white"), legend.title = element_blank()) +
      scale_color_manual(values=c("Q1" = "#4575B4",
                                  "Q2" = "#91BFDB",
                                  "Q3" = "cyan",
                                  "Q4" = "#FEE090",
                                  "Q5" = "#FC8D59",
                                  "Q6" = "#D73027"),
                         breaks=c("Q1", "Q2", "Q3", "Q4", "Q5", "Q6"))
  )
  dev.off()
}



hp1_na <- df[c("half-life_LOS", "RNA-seq_total_+_R1", "PCA_eigen1")]
colnames(hp1_na) <- c("hl", "RNAseq", "PCA_eigen1")
hp1 <- na.omit(hp1_na)

# Compartments
A <- hp1[hp1$PCA_eigen1 > 0,]
B <- hp1[hp1$PCA_eigen1 <= 0,]

plot_RNAseq(A, "A")
plot_RNAseq(B, "B")

