library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD_K562", "PCA_eigen1")]
colnames(hp1_na) <- c("hl", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD", "PCA_eigen1")
hp1 <- na.omit(hp1_na)

# Get upper quartile data 
upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1)[4],]
upq_CBX3 <- hp1[hp1$CBX3_Myers > quantile(hp1$CBX3_Myers)[4],]
upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5)[4],]
upq_LAD <- hp1[hp1$LAD > quantile(hp1$LAD)[4],]

# Compartments
A_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 > 0,]
B_upq_CBX1 <- upq_CBX1[upq_CBX1$PCA_eigen1 <= 0,]
A_upq_CBX3 <- upq_CBX3[upq_CBX3$PCA_eigen1 > 0,]
B_upq_CBX3 <- upq_CBX3[upq_CBX3$PCA_eigen1 <= 0,]
A_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 > 0,]
B_upq_CBX5 <- upq_CBX5[upq_CBX5$PCA_eigen1 <= 0,]
A_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 > 0,]
B_upq_LAD <- upq_LAD[upq_LAD$PCA_eigen1 <= 0,]


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

# B Compartments
pdf(paste("half-life_HP1_LAD_density_genome_Bcomp.pdf", sep=""),
          width=7, height=6)
      print(
      ggplot() + 
        geom_line(data=wave_fun(B_upq_CBX1$hl, 1, length(B_upq_CBX1$hl)), aes(x=x, y=Frequency, color="CBX1")) +
        #geom_histogram(data=B_upq_CBX3, aes(x=hl, color="CBX3"), fill="white", binwidth = 1) +
        geom_line(data=wave_fun(B_upq_CBX3$hl, 1, length(B_upq_CBX3$hl)), aes(x=x, y=Frequency, color="CBX3")) +
        geom_line(data=wave_fun(B_upq_CBX5$hl, 1, length(B_upq_CBX5$hl)), aes(x=x, y=Frequency, color="CBX5")) +
        geom_line(data=wave_fun(B_upq_LAD$hl, 1, length(B_upq_LAD$hl)), aes(x=x, y=Frequency, color="LAD")) +
        xlim(25,165) + ylim(0,700) +
        xlab(bquote("t"[1/2] ~ "(minutes)")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
              panel.background = element_blank(), axis.line = element_line(colour = "black"),
              legend.key = element_rect(fill="white"), legend.title = element_blank())
      )
dev.off()
# A Compmartments     
pdf(paste("half-life_HP1_LAD_density_genome_Acomp.pdf", sep=""),
    width=7, height=6)
print(
  ggplot() + 
    geom_line(data=wave_fun(A_upq_CBX1$hl, 1, length(A_upq_CBX1$hl)), aes(x=x, y=Frequency, color="CBX1")) +
    #geom_histogram(data=B_upq_CBX3, aes(x=hl, color="CBX3"), fill="white", binwidth = 1) +
    geom_line(data=wave_fun(A_upq_CBX3$hl, 1, length(A_upq_CBX3$hl)), aes(x=x, y=Frequency, color="CBX3")) +
    geom_line(data=wave_fun(A_upq_CBX5$hl, 1, length(A_upq_CBX5$hl)), aes(x=x, y=Frequency, color="CBX5")) +
    geom_line(data=wave_fun(A_upq_LAD$hl, 1, length(A_upq_LAD$hl)), aes(x=x, y=Frequency, color="LAD")) +
    xlim(25,165) + ylim(0,700) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank())
)
dev.off()
      

