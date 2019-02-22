library(ggplot2)
library(gridExtra)
df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

B_na <- df[c("half-life_LOS", "CBX1_R1", "CBX5", "LAD_K562", 
               "H3K9me3_R1","PCA_eigen1")]
colnames(B_na) <- c("hl", "CBX1", "CBX5", "LAD", "H3K9me3", "PCA_eigen1")
B <- na.omit(B_na)

# Get upper quartile data 
upq_CBX1 <- B[B$CBX1 > quantile(B$CBX1)[4],]
upq_CBX5 <- B[B$CBX5 > quantile(B$CBX5)[4],]
upq_LAD <- B[B$LAD > quantile(B$LAD)[4],]
upq_H3K9 <- B[B$H3K9me3 > quantile(B$H3K9me3)[4],]

# B compartments
B_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 < 0,]
B_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 < 0,]
B_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 < 0,]
B_upq_H3K9 <- upq_H3K9[upq_H3K9$PCA_eigen1 < 0,]

A_na <- df[c("half-life_LOS", "H3K27ac_R1", "H3K9ac_R1", "SON_TSA-seq", 
             "PCA_eigen1")]

colnames(A_na) <- c("hl", "H3K27ac", "H3K9ac", "SON.TSA.seq", "PCA_eigen1")
A <- na.omit(A_na)

# Get upper quartile data 
upq_H3K27ac <- A[A$H3K27ac > quantile(A$H3K27ac)[4],]
upq_H3K9ac <- A[A$H3K9ac > quantile(A$H3K9ac)[4],]
upq_SON <- A[A$SON.TSA.seq > quantile(A$SON.TSA.seq)[4],]

# B compartments
A_upq_H3K27ac <- upq_H3K27ac[upq_H3K27ac$PCA_eigen1 > 0,]
A_upq_H3K9ac <- upq_H3K9ac[upq_H3K9ac$PCA_eigen1 > 0,]
A_upq_SON <- upq_SON[upq_SON$PCA_eigen1 > 0,]



pdf("B_only_A_only_cumulative.pdf", width=5, height=5)

p1 <-  ggplot() + 
    stat_ecdf(data=B_upq_CBX1, aes(x=hl, color="CBX1"), geom="line") +
    stat_ecdf(data=B_upq_CBX5, aes(x=hl, color="CBX5"), geom="line") +
    stat_ecdf(data=B_upq_LAD, aes(x=hl, color="LAD"), geom="line") +
    #stat_ecdf(data=B_upq_H3K9, aes(x=hl, color="H3K9me3"), geom="line") +
    xlim(40, 120) + ylim(0,1) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab("Fraction of Data") +  
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank(),
          legend.position = c(0.85,0.35)) +
    scale_color_manual(values=c("CBX1" = "cyan",
                                "CBX5" = "lightskyblue2",
                                "LAD" = "blue"), 
                       breaks=c("CBX1", "CBX5", "LAD"))

p2 <-  ggplot() + 
  stat_ecdf(data=A_upq_H3K27ac, aes(x=hl, color="H3K27ac"), geom="line") +
  stat_ecdf(data=A_upq_H3K9ac, aes(x=hl, color="H3K9ac"), geom="line") +
  stat_ecdf(data=A_upq_SON, aes(x=hl, color="SON.TSA.seq"), geom="line") +
  xlim(40, 120) + ylim(0,1) +
  xlab("") + ylab("Fraction of Data") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill="white"), legend.title = element_blank(),
        legend.position = c(0.85,0.35)) +
  scale_color_manual(values=c("H3K27ac" = "orange",
                              "H3K9ac" = "magenta",
                              "SON.TSA.seq" = "red3"), 
                     breaks=c("H3K27ac", "H3K9ac", "SON.TSA.seq")) 

grid.arrange(p2, p1, nrow=2)

dev.off()



