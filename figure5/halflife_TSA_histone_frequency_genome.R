library(ggplot2)
library(RColorBrewer)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "SON_TSA-seq", "H3K9me1_R1", "H3K4me1_R1", "H3K27ac_R1","PCA_eigen1")]
colnames(hp1_na) <- c("hl", "TSAseq", "H3K9me1", "H3K4me1", "H3K27ac", "PCA_eigen1")
hp1 <- na.omit(hp1_na)

# # Get upper quartile data 
upq_TSAseq <- hp1[hp1$TSAseq > quantile(hp1$TSAseq)[4],]
upq_H3K9me1 <- hp1[hp1$H3K9me1 > quantile(hp1$H3K9me1)[4],]
upq_H3K4me1 <- hp1[hp1$H3K4me1 > quantile(hp1$H3K4me1)[4],]
upq_H3K27ac <- hp1[hp1$H3K27ac > quantile(hp1$H3K27ac)[4],]


# Compartments
A_upq_TSAseq <- upq_TSAseq[upq_TSAseq$PCA_eigen1 > 0,]
B_upq_TSAseq <- upq_TSAseq[upq_TSAseq$PCA_eigen1 <= 0,]
A_upq_H3K9me1 <- upq_H3K9me1[upq_H3K9me1$PCA_eigen1 > 0,]
B_upq_H3K9me1 <- upq_H3K9me1[upq_H3K9me1$PCA_eigen1 <= 0,]
A_upq_H3K4me1 <- upq_H3K4me1[upq_H3K4me1$PCA_eigen1 > 0,]
B_upq_H3K4me1 <- upq_H3K4me1[upq_H3K4me1$PCA_eigen1 <= 0,]
A_upq_H3K27ac <- upq_H3K27ac[upq_H3K27ac$PCA_eigen1 > 0,]
B_upq_H3K27ac <- upq_H3K27ac[upq_H3K27ac$PCA_eigen1 <= 0,]


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
# B Compartments
pdf(paste("B/half-life_TSA_histone_density_genome_B.pdf", sep=""),
          width=7, height=6)
      print(
      ggplot() + 
        geom_line(data=wave_fun(B_upq_TSAseq$hl, 1, length(B_upq_TSAseq$hl)), aes(x=x, y=Frequency, color="TSA-seq")) +
        geom_line(data=wave_fun(B_upq_H3K9me1$hl, 1, length(B_upq_H3K9me1$hl)), aes(x=x, y=Frequency, color="H3K9me1")) +
        #geom_line(data=wave_fun(B_upq_H3K4me1$hl, 1, length(B_upq_H3K4me1$hl)), aes(x=x, y=Frequency, color="H3K4me1")) +
        geom_line(data=wave_fun(B_upq_H3K27ac$hl, 1, length(B_upq_H3K27ac$hl)), aes(x=x, y=Frequency, color="H3K27ac")) +
        xlim(30,150) +
        xlab(bquote("t"[1/2] ~ "(minutes)")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill="white"), legend.title = element_blank()) +
        scale_color_manual(values=c("TSA-seq" = "red3",
                                    "H3K9me1" = "pink2",
                                    #"H3K4me1" = "cyan",
                                    "H3K27ac" = "magenta"), 
                           breaks=c("TSA-seq", "H3K9me1", "H3K27ac"))
      )
dev.off()

pdf(paste("B/half-life_TSA_histone_cumulative_genome_B.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    stat_ecdf(data=B_upq_TSAseq, aes(x=hl, color="TSA-seq"), geom="line") +
    stat_ecdf(data=B_upq_H3K9me1, aes(x=hl, color="H3K9me1"), geom="line") +
    #stat_ecdf(data=B_upq_H3K4me1, aes(x=hl, color="H3K4me1"), geom="line") +
    stat_ecdf(data=B_upq_H3K27ac, aes(x=hl, color="H3K27ac"), geom="line") +
    xlim(40,130) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab(bquote("F"["n"] ~ "(x)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("TSA-seq" = "red3",
                                "H3K9me1" = "pink2",
                                #"H3K4me1" = "cyan",
                                "H3K27ac" = "magenta"), 
                       breaks=c("TSA-seq", "H3K9me1", "H3K27ac"))
)
dev.off()


####################################################################################################################
# A Compmartments     
pdf(paste("A/half-life_TSA_histone_density_genome_A.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    geom_line(data=wave_fun(A_upq_TSAseq$hl, 1, length(A_upq_TSAseq$hl)), aes(x=x, y=Frequency, color="TSA-seq")) +
    geom_line(data=wave_fun(A_upq_H3K9me1$hl, 1, length(A_upq_H3K9me1$hl)), aes(x=x, y=Frequency, color="H3K9me1")) +
    #geom_line(data=wave_fun(A_upq_H3K4me1$hl, 1, length(A_upq_H3K4me1$hl)), aes(x=x, y=Frequency, color="H3K4me1")) +
    geom_line(data=wave_fun(A_upq_H3K27ac$hl, 1, length(A_upq_H3K27ac$hl)), aes(x=x, y=Frequency, color="H3K27ac")) +
    
    xlim(30,100) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("TSA-seq" = "red3",
                                "H3K9me1" = "pink2",
                                #"H3K4me1" = "cyan",
                                "H3K27ac" = "magenta"), 
                       breaks=c("TSA-seq", "H3K9me1", "H3K27ac"))
)
dev.off()
      
pdf(paste("A/half-life_TSA_histone_cumulative_genome_A.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    stat_ecdf(data=A_upq_TSAseq, aes(x=hl, color="TSA-seq"), geom="line") +
    stat_ecdf(data=A_upq_H3K9me1, aes(x=hl, color="H3K9me1"), geom="line") +
    #stat_ecdf(data=A_upq_H3K4me1, aes(x=hl, color="H3K4me1"), geom="line") +
    stat_ecdf(data=A_upq_H3K27ac, aes(x=hl, color="H3K27ac"), geom="line") +
    xlim(35, 95) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab(bquote("F"["n"] ~ "(x)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("TSA-seq" = "red3",
                                "H3K9me1" = "pink2",
                                #"H3K4me1" = "cyan
                                "H3K27ac" = "magenta"), 
                       breaks=c("TSA-seq", "H3K9me1", "H3K27ac"))
)
dev.off()

