library(ggplot2)
F_df <- read.table("FPKM_40kb.bedGraph",
                   sep="\t", header=FALSE)
colnames(F_df) <- c("chrom", "start", "end", "FPKM")
mat <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt",
                sep="\t", header=TRUE)

FPKM <- F_df$FPKM
F_df$FPKM <- log2(FPKM)
F_df$FPKM[!is.finite(F_df$FPKM)] <- NA


df <- merge(F_df, mat, by=c("chrom", "start", "end"))


pdf("SON_FPKM_halflife_scatter.pdf", width=5, height=4)
  ggplot(df, aes(x=half.life_LOS, y=FPKM, color=SON_TSA.seq)) + geom_point(size=0.1) +
  scale_color_gradient(low="dodgerblue", high="red") +
  ylim(-10,10) + xlim(35,130) + xlab(bquote("t"[1/2]~"(minutes)")) +
  ylab(bquote("log"[2]~"(FPKM)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
  panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

png("SON_FPKM_halflife_scatter.png", width=1300, height=1000, res=300)
ggplot(df, aes(x=half.life_LOS, y=FPKM, color=SON_TSA.seq)) + geom_point(size=0.3) +
  scale_color_gradient(low="dodgerblue", high="red") +
  ylim(-10,10) + xlim(35,130) + xlab(bquote("t"[1/2]~"(minutes)")) +
  ylab(bquote("log"[2]~"(FPKM)"))+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
         