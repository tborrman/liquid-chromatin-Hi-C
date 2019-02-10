library(ggplot2)

df <- read.table("C:/cygwin64/home/Tyler/Research/digest/figure4/feature_matrix_v3_40kb.txt",
           sep="\t", header=TRUE, check.names=FALSE)

stratify_RNAseq <- function(x) {
  g <- c()
  for (r in x) {
    if (r >=0 & r < 0.01) {
      g <- c(g, '0-0.01')
    }
    else if (r >= 0.01 & r < 0.1) {
      g <- c(g, '0.01-0.1')
    }
    else if (r >= 0.1 & r < 1) {
      g <- c(g, '0.1-1')
    }
    else if (r >= 1 & r < 10) {
      g <- c(g, '1-10')
    }
    else if (r >=10) {
      g <- c(g, '>10')
    }
    else {
      print(r)
      stop()
    }
  }
  fg <- factor(g, levels=c("0-0.01", "0.01-0.1", "0.1-1", "1-10", ">10"), ordered=TRUE)
  return(fg)
}

##############################################################################
# RNA-seq_total_+ density plot
df_NA <- df[c("half-life_LOS", "RNA-seq_total_+_R1")]
colnames(df_NA) <- c("hl", "RNAseq")
df_hl_RNA <- na.omit(df_NA)
RNAseq_signal <- stratify_RNAseq(df_hl_RNA$RNAseq)
df_hl_RNA_labs <- cbind(df_hl_RNA, RNAseq_signal)

png("genome_RNAseq/half-life_RNA-seq_total_+_R1_density.png", width=2000, height=1400, res=300)
ggplot(df_hl_RNA_labs, aes(x=hl, color=RNAseq_signal)) +
  geom_density() +
  xlim(0,400) +
  ylim(0,0.025) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) +
  ggtitle("total RNA-seq (+) strand") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
##############################################################################
####################################################################################
# A compartment only
# RNA-seq_total_+ density plot
df_NA <- df[c("half-life_LOS", "RNA-seq_total_+_R1", "PCA_eigen1")]
colnames(df_NA) <- c("hl", "RNAseq", "PC1")
df_hl_RNA <- na.omit(df_NA)
RNAseq_signal <- stratify_RNAseq(df_hl_RNA$RNAseq)
df_hl_RNA_labs <- cbind(df_hl_RNA, RNAseq_signal)
df_A <- df_hl_RNA_labs[df_hl_RNA_labs$PC1 > 0,]
df_B <- df_hl_RNA_labs[df_hl_RNA_labs$PC1 <= 0,]

png("genome_RNAseq/half-life_RNA-seq_total_+_R1_density_Acomp.png", width=2000, height=1400, res=300)
ggplot(df_A, aes(x=hl, color=RNAseq_signal)) +
  geom_density() +
  xlim(0,400) +
  ylim(0,0.025) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) +
  ggtitle("total RNA-seq (+) strand A compartment") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()

png("genome_RNAseq/half-life_RNA-seq_total_+_R1_density_Bcomp.png", width=2000, height=1400, res=300)
ggplot(df_B, aes(x=hl, color=RNAseq_signal)) +
  geom_density() +
  xlim(0,400) +
  ylim(0,0.025) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) +
  ggtitle("total RNA-seq (+) strand B compartment") + 
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"))
dev.off()
########################################################################################################
# percentiles
q1 <- df_hl_RNA_labs[df_hl_RNA_labs$RNAseq <= quantile(df_hl_RNA_labs$RNAseq)[2],]
q2 <- df_hl_RNA_labs[df_hl_RNA_labs$RNAseq > quantile(df_hl_RNA_labs$RNAseq)[2] & 
                       df_hl_RNA_labs$RNAseq <= quantile(df_hl_RNA_labs$RNAseq)[3],]
q3 <- df_hl_RNA_labs[df_hl_RNA_labs$RNAseq > quantile(df_hl_RNA_labs$RNAseq)[3] & 
                       df_hl_RNA_labs$RNAseq <= quantile(df_hl_RNA_labs$RNAseq)[4],]
q4 <- df_hl_RNA_labs[df_hl_RNA_labs$RNAseq > quantile(df_hl_RNA_labs$RNAseq)[4],]

pdf("genome_RNAseq/half-life_RNA-seq_total_+_R1_density_quartiles.pdf", 
    width=7, height=6)
  ggplot() +
    geom_density(data=q1, aes(x=hl, color="Q1")) +
    geom_density(data=q2, aes(x=hl, color="Q2")) +
    geom_density(data=q3, aes(x=hl, color="Q3")) +
    geom_density(data=q4, aes(x=hl, color="Q4")) +
    xlim(0,400) +
    ylim(0,0.020) +
    xlab(bquote("t"[1/2] ~ "(minutes)")) + 
    theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
          panel.background = element_blank(), axis.line = element_line(colour = "black"), 
          legend.key = element_rect(fill="white")) +
    labs(color="RNA-seq signal") + 
    scale_color_manual(values=rev(brewer.pal(n = 4, name ="RdYlBu")))
dev.off()

png("genome_RNAseq/half-life_RNA-seq_total_+_R1_density_scatterplot.png",
    width=2000, height=2000, res=300)
  plot(df_hl_RNA_labs$hl, log10(df_hl_RNA_labs$RNAseq), pch=20, col=rgb(0,0,1,0.1),
       xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="RNA-seq signal", 
       xlim=c(25,350))
dev.off()