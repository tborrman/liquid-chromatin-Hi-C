df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("half-life_LOS", "RNA-seq_total_+_R1", "PCA_eigen1")]
colnames(hp1_na) <- c("hl", "RNAseq", "PCA_eigen1")
hp1 <- na.omit(hp1_na)

# Compartments
A <- hp1[hp1$PCA_eigen1 > 0,]
B <- hp1[hp1$PCA_eigen1 <= 0,]

# Half-life vs RNA-seq
png("half-life_RNA-seq_total_+_R1_density_scatterplot.png",
    width=2000, height=2000, res=300)
plot(hp1$hl, hp1$RNAseq, pch=20, col=rgb(0,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="RNA-seq 40kb mean(RPM)",
     ylim = c(0,60), xlim= c(20,150)
)
text(120, 40, bquote(rho ~ " = "~ .(round(cor(hp1$hl, hp1$RNAseq, method="spearman"), 2))))
dev.off()

# Half-life vs log(RNA-seq)
png("half-life_log_RNA-seq_total_+_R1_density_scatterplot.png",
    width=2000, height=2000, res=300)
plot(hp1$hl, log2(hp1$RNAseq), pch=20, col=rgb(0,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="RNA-seq 40kb log2(mean(RPM))",
      ylim=c(-7,7), xlim= c(20,150)
)
text(120, 4, bquote(rho ~ " = "~ .(round(cor(hp1$hl, hp1$RNAseq, method="spearman"), 2))))
dev.off()

# A compartment half-life vs log(RNA-seq)
png("half-life_log_RNA-seq_total_+_R1_density_scatterplot_A.png",
    width=2000, height=2000, res=300)
plot(A$hl, log2(A$RNAseq), pch=20, col=rgb(1,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="RNA-seq 40kb log2(mean(RPM))",
     ylim=c(-7,7), xlim= c(20,150)
)
text(120, 4, bquote(rho ~ " = "~ .(round(cor(A$hl, A$RNAseq, method="spearman"), 2))))
dev.off()

# B compartment half-life vs log(RNA-seq)
png("half-life_log_RNA-seq_total_+_R1_density_scatterplot_B.png",
    width=2000, height=2000, res=300)
plot(B$hl, log2(B$RNAseq), pch=20, col=rgb(0,0,1,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="RNA-seq 40kb log2(mean(RPM))",
     ylim=c(-7,7), xlim= c(20,150)
)
text(120, 4, bquote(rho ~ " = "~ .(round(cor(B$hl, B$RNAseq, method="spearman"), 2))))
dev.off()

# Boxplots
highRNA_A <- A[A$RNAseq > 1,]
lowRNA_A <- A[A$RNAseq <= 1,]
highRNA_B <- B[B$RNAseq > 1,]
lowRNA_B <- B[B$RNAseq <= 1,]

png('boxplot_RNAseq.png', width=1000, height=2000, res=300)
boxplot(highRNA_A$hl, lowRNA_A$hl, highRNA_B$hl, lowRNA_B$hl, 
        col=c("red", "red", "blue", "blue"), outline=FALSE,
        names=c("RPM > 1", "RPM <= 1", "RPM > 1", "RPM <= 1"),
        ylab = bquote("t"[1/2] ~ "(minutes)"), las=2, ylim=c(30,140))
text(1,median(highRNA_A$hl, na.rm=TRUE) + 2, paste("n=", length(highRNA_A$hl), sep=""), 
     col= "white", cex = 0.5)
text(2,median(lowRNA_A$hl, na.rm=TRUE) + 2, paste("n=", length(lowRNA_A$hl), sep=""), 
     col= "white", cex = 0.5)
text(3,median(highRNA_B$hl, na.rm=TRUE) + 2, paste("n=", length(highRNA_B$hl), sep=""), 
     col= "white", cex = 0.5)
text(4,median(lowRNA_B$hl, na.rm=TRUE) + 2, paste("n=", length(lowRNA_B$hl), sep=""), 
     col= "white", cex = 0.5)
legend("topleft", c("A", "B"), fill=c("red", "blue"), bty="n")
dev.off()

t.test(highRNA_A$hl, lowRNA_A$hl)
t.test(highRNA_B$hl, lowRNA_B$hl)




