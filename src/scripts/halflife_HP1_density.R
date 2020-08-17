library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v3/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "CBX1_R1", "CBX3_Myers", "CBX5")]
colnames(hp1_na) <- c("hl", "CBX1_R1", "CBX3_Myers", "CBX5")
hp1 <- na.omit(hp1_na)

# Get upper quartile data 
upq_CBX1 <- hp1[hp1$CBX1_R1 > quantile(hp1$CBX1_R1)[4],]
upq_CBX3 <- hp1[hp1$CBX3_Myers > quantile(hp1$CBX3_Myers)[4],]
upq_CBX5 <- hp1[hp1$CBX5 > quantile(hp1$CBX5)[4],]

png("half-life_HP1_density.png", width=2000, height=1400, res=300)
      ggplot() + 
        geom_density(data=upq_CBX1, aes(x=hl, color="CBX1")) +
        geom_density(data=upq_CBX3, aes(x=hl, color="CBX3")) +
        geom_density(data=upq_CBX5, aes(x=hl, color="CBX5")) + 
        xlim(0,400) + ylim(0,0.02) +
        xlab(bquote("t"[1/2] ~ "(minutes)")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank())
dev.off()
