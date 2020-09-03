library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v3/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

# Compartment density plot
##############################################################################
compartment <- ifelse(df$PCA_eigen1 > 0, 'A', 'B')
hl <- df$`half-life_LOS`
df_na<- data.frame(hl, compartment)
df_d <- df_na[!is.na(df_na$compartment),]

png("half-life_compartment_density.png", width=2000, height=1400, res=300)
ggplot(df_d, aes(x=hl, fill=compartment, color=compartment)) +
  geom_density(alpha=0.5) +
  xlim(0,400) +
  ylim(0, 0.015) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
      panel.background = element_blank(), axis.line = element_line(colour = "black")) +
  scale_fill_manual(values=c("red", "blue")) +
  scale_color_manual(values=c("red", "blue"))
dev.off()
####################################################3