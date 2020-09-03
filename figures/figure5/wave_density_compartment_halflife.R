library(ggplot2)
df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

c_na <- df[c("half-life_LOS", "PCA_eigen1")]
colnames(c_na) <- c("hl", "PCA_eigen1")

c_df <- na.omit(c_na)

A_df <- c_df[c_df$PCA_eigen1 > 0,]
B_df <- c_df[c_df$PCA_eigen1 < 0,]

pdf("wave_density_compartment_halflife.pdf", width=5, height=4)
print(
  ggplot() + 
    stat_density(data=A_df, aes(x=hl, color="A"), geom="line") +
    stat_density(data=B_df , aes(x=hl, color="B"), geom="line") +
    #xlim(0,400) + ylim(0,0.020) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) +
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"),
          legend.key = element_rect(fill="white"), legend.title = element_blank()) +
    scale_color_manual(values=c("A"="red", "B"="blue"),
                       breaks=c("A", "B"))
)
dev.off()
