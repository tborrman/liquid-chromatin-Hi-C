library(ggplot2)
library(RColorBrewer)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "CBX1_R1", "CBX5", "LAD_K562","PCA_eigen1")]
colnames(hp1_na) <- c("hl", "CBX1_R1", "CBX5", "LAD", "PCA_eigen1")
hp1 <- na.omit(hp1_na)

# # Get upper quartile data 
upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1)[4],]
upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5)[4],]
upq_LAD <- hp1[hp1$LAD > quantile(hp1$LAD)[4],]
#upq_H3K9 <- hp1[hp1$H3K9me3 > quantile(hp1$H3K9me3)[4],]

# Get upper quartile data 
# upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1, 0.95),]
# upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5, 0.95),]
# upq_LAD <- hp1[hp1$LAD > quantile(hp1$LAD, 0.95),]
# upq_H3K9 <- hp1[hp1$H3K9me3 > quantile(hp1$H3K9me3, 0.95),]

# Compartments
A_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 > 0,]
B_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 <= 0,]
A_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 > 0,]
B_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 <= 0,]
A_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 > 0,]
B_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 <= 0,]
#A_upq_H3K9 <- upq_H3K9[upq_H3K9$PCA_eigen1 > 0,]
#B_upq_H3K9 <- upq_H3K9[upq_H3K9$PCA_eigen1 <= 0,]

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
pdf(paste("B/half-life_HP1_LAD_density_genome_B.pdf", sep=""),
          width=7, height=6)
      print(
      ggplot() + 
        geom_line(data=wave_fun(B_upq_CBX1$hl, 1, length(B_upq_CBX1$hl)), aes(x=x, y=Frequency, color="CBX1")) +
        #geom_histogram(data=B_upq_CBX3, aes(x=hl, color="CBX3"), fill="white", binwidth = 1) +
        geom_line(data=wave_fun(B_upq_CBX5$hl, 1, length(B_upq_CBX5$hl)), aes(x=x, y=Frequency, color="CBX5")) +
        #geom_line(data=wave_fun(B_upq_H3K9$hl, 1, length(B_upq_H3K9$hl)), aes(x=x, y=Frequency, color="H3K9me3")) +
        geom_line(data=wave_fun(B_upq_LAD$hl, 1, length(B_upq_LAD$hl)), aes(x=x, y=Frequency, color="LADs")) +
        xlim(40,150) + ylim(0,600) +
        xlab(bquote("t"[1/2] ~ "(minutes)")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill="white"), legend.title = element_blank()) +
        #scale_color_manual(values=rev(brewer.pal(n = 3, name ="RdYlBu")))
        scale_color_manual(values=c("dodgerblue", "mediumorchid4", "red"))
      )
dev.off()

pdf(paste("B/half-life_HP1_LAD_cumulative_genome_B.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    stat_ecdf(data=B_upq_CBX1, aes(x=hl, color="CBX1"), geom="line") +
    stat_ecdf(data=B_upq_CBX5, aes(x=hl, color="CBX5"), geom="line") +
    #stat_ecdf(data=B_upq_H3K9, aes(x=hl, color="H3K9me3"), geom="line") +
    stat_ecdf(data=B_upq_LAD, aes(x=hl, color="LADs"), geom="line") +
    xlim(50,140) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab(bquote("F"["n"] ~ "(x)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    #scale_color_manual(values=rev(brewer.pal(n = 3, name ="RdYlBu")))
    scale_color_manual(values=c("dodgerblue", "mediumorchid4", "red"))
)
dev.off()


####################################################################################################################
# A Compmartments     
pdf(paste("A/half-life_HP1_LAD_density_genome_A.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    geom_line(data=wave_fun(A_upq_CBX1$hl, 1, length(A_upq_CBX1$hl)), aes(x=x, y=Frequency, color="CBX1")) +
    #geom_histogram(data=B_upq_CBX3, aes(x=hl, color="CBX3"), fill="white", binwidth = 1) +
    geom_line(data=wave_fun(A_upq_CBX5$hl, 1, length(A_upq_CBX5$hl)), aes(x=x, y=Frequency, color="CBX5")) +
    #geom_line(data=wave_fun(A_upq_H3K9$hl, 1, length(A_upq_H3K9$hl)), aes(x=x, y=Frequency, color="H3K9me3")) +
    geom_line(data=wave_fun(A_upq_LAD$hl, 1, length(A_upq_LAD$hl)), aes(x=x, y=Frequency, color="LADs")) +
    
    xlim(30,150) + ylim(0,600) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    #scale_color_manual(values=rev(brewer.pal(n = 3, name ="RdYlBu")))
    scale_color_manual(values=c("dodgerblue", "mediumorchid4", "red"))
)
dev.off()
      

