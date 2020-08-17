library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v3/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

df_na <- df[c("half-life_LOS", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD_clone14",
              "CBX8", "RNF2", "BMI1", "SUZ12")]
colnames(df_na) <- c("hl", "CBX1_R1", "CBX3_Myers", "CBX5", "LAD_clone14",
                     "CBX8", "RNF2", "BMI1", "SUZ12")
d <- na.omit(df_na)

# Get upper quartile data 
upq_CBX1 <- d[d$CBX1_R1 > quantile(d$CBX1_R1)[4],]
upq_CBX3 <- d[d$CBX3_Myers > quantile(d$CBX3_Myers)[4],]
upq_CBX5 <- d[d$CBX5 > quantile(d$CBX5)[4],]
upq_LAD <- d[d$LAD_clone14 > quantile(d$LAD_clone14)[4],]
upq_CBX8 <- d[d$CBX8 > quantile(d$CBX8)[4],]
upq_RNF2 <- d[d$RNF2 > quantile(d$RNF2)[4],]
upq_BMI1 <- d[d$BMI1 > quantile(d$BMI1)[4],]
upq_SUZ12 <- d[d$SUZ12 > quantile(d$SUZ12)[4],]
  
png("half-life_LAD_PcG_density.png", width=2000, height=1400, res=300)
      ggplot() + 
        geom_density(data=upq_LAD, aes(x=hl, color="LAD")) +
        geom_density(data=upq_CBX8, aes(x=hl, color="CBX8 (PcG)")) +
        geom_density(data=upq_RNF2, aes(x=hl, color="RNF2 (PcG)")) + 
        geom_density(data=upq_BMI1, aes(x=hl, color="BMI1 (PcG)")) +
        geom_density(data=upq_SUZ12, aes(x=hl, color="SUZ12 (PcG)")) +
        xlim(0,400) + ylim(0,0.02) +
        xlab(bquote("t"[1/2] ~ "(minutes)")) +
        scale_color_discrete(breaks=c("LAD","CBX8 (PcG)","RNF2 (PcG)",
                                      "BMI1 (PcG)","SUZ12 (PcG)")) +
        theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.title = element_blank())
dev.off()

png("half-life_LAD_PcG_HP1_density.png", width=2000, height=1400, res=300)
ggplot() + 
  geom_density(data=upq_LAD, aes(x=hl, color="LAD")) +
  geom_density(data=upq_CBX8, aes(x=hl, color="CBX8 (PcG)")) +
  geom_density(data=upq_RNF2, aes(x=hl, color="RNF2 (PcG)")) + 
  geom_density(data=upq_BMI1, aes(x=hl, color="BMI1 (PcG)")) +
  geom_density(data=upq_SUZ12, aes(x=hl, color="SUZ12 (PcG)")) +
  geom_density(data=upq_CBX1, aes(x=hl, color="CBX1 (HP1)")) + 
  geom_density(data=upq_CBX3, aes(x=hl, color="CBX3 (HP1)")) +
  geom_density(data=upq_CBX5, aes(x=hl, color="CBX5 (HP1)")) +
  xlim(0,400) + ylim(0,0.02) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) +
  scale_color_discrete(breaks=c("LAD","CBX8 (PcG)","RNF2 (PcG)",
                                "BMI1 (PcG)","SUZ12 (PcG)", "CBX1 (HP1)",
                                "CBX3 (HP1)","CBX5 (HP1)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.title = element_blank())
dev.off()
