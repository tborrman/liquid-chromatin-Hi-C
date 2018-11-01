library(ggplot2)
library(gridExtra)
library(RColorBrewer)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "SON_TSA-seq", "H3K9me1_R1", "H3K4me1_R1", "H3K27ac_R1","PCA_eigen1",
               "CBX1_R1", "CBX5", "LAD_K562", "H3K9me3_R1")]
colnames(hp1_na) <- c("hl", "TSAseq", "H3K9me1", "H3K4me1", "H3K27ac", "PCA_eigen1",
                      "CBX1_R1", "CBX5", "LAD", "H3K9me3")
hp1 <- na.omit(hp1_na)

# # Get upper quartile data 
upq_TSAseq <- hp1[hp1$TSAseq > quantile(hp1$TSAseq)[4],]
upq_H3K9me1 <- hp1[hp1$H3K9me1 > quantile(hp1$H3K9me1)[4],]
upq_H3K4me1 <- hp1[hp1$H3K4me1 > quantile(hp1$H3K4me1)[4],]
upq_H3K27ac <- hp1[hp1$H3K27ac > quantile(hp1$H3K27ac)[4],]
# # Get upper quartile data 
upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1)[4],]
upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5)[4],]
upq_LAD <- hp1[hp1$LAD > quantile(hp1$LAD)[4],]
upq_H3K9 <- hp1[hp1$H3K9me3 > quantile(hp1$H3K9me3)[4],]


# Compartments
A_upq_TSAseq <- upq_TSAseq[upq_TSAseq$PCA_eigen1 > 0,]
B_upq_TSAseq <- upq_TSAseq[upq_TSAseq$PCA_eigen1 <= 0,]
A_upq_H3K9me1 <- upq_H3K9me1[upq_H3K9me1$PCA_eigen1 > 0,]
B_upq_H3K9me1 <- upq_H3K9me1[upq_H3K9me1$PCA_eigen1 <= 0,]
A_upq_H3K4me1 <- upq_H3K4me1[upq_H3K4me1$PCA_eigen1 > 0,]
B_upq_H3K4me1 <- upq_H3K4me1[upq_H3K4me1$PCA_eigen1 <= 0,]
A_upq_H3K27ac <- upq_H3K27ac[upq_H3K27ac$PCA_eigen1 > 0,]
B_upq_H3K27ac <- upq_H3K27ac[upq_H3K27ac$PCA_eigen1 <= 0,]
A_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 > 0,]
B_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 <= 0,]
A_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 > 0,]
B_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 <= 0,]
A_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 > 0,]
B_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 <= 0,]
A_upq_H3K9 <- upq_H3K9[upq_H3K9$PCA_eigen1 > 0,]
B_upq_H3K9 <- upq_H3K9[upq_H3K9$PCA_eigen1 <= 0,]


pdf(paste("half-life_cumulative.pdf", sep=""),
    width=8, height=6)

p1 <-  ggplot() + 
    stat_ecdf(data=A_upq_TSAseq, aes(x=hl, color="TSA-seq"), geom="line") +
    stat_ecdf(data=A_upq_H3K9me1, aes(x=hl, color="H3K9me1"), geom="line") +
    stat_ecdf(data=A_upq_H3K27ac, aes(x=hl, color="H3K27ac"), geom="line") +
    xlim(40, 125) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab("Fraction of Data") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("TSA-seq" = "red3",
                                "H3K9me1" = "pink2",
                                "H3K27ac" = "magenta"), 
                       breaks=c("TSA-seq", "H3K9me1", "H3K27ac"))
p2 <- ggplot() + 
    stat_ecdf(data=B_upq_CBX1, aes(x=hl, color="CBX1"), geom="line") +
    stat_ecdf(data=B_upq_CBX5, aes(x=hl, color="CBX5"), geom="line") +
    stat_ecdf(data=B_upq_LAD, aes(x=hl, color="LADs"), geom="line") +
    xlim(40,125) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab("Fraction of Data") +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("CBX1"= "cyan", 
                                "CBX5" = "#91BFDB", 
                                "LADs" = "blue"))
grid.arrange(p1, p2, nrow = 2)

dev.off()



