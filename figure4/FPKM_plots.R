library(ggplot2)
df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

hp1_na <- df[c("chrom", "start", "end", "half-life_LOS", "PCA_eigen1")]
colnames(hp1_na) <- c("chrom", "start", "end", "hl", "PCA_eigen1")
fpkm_df <- read.table("FPKM_40kb.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "FPKM"))

mdf <- merge(hp1_na, fpkm_df, by=c("chrom", "start", "end"))

#mdf$FPKM[is.na(mdf$FPKM)] <- 0
mdf <- na.omit(mdf)

# Compartments
A <- mdf[mdf$PCA_eigen1 > 0,]
B <- mdf[mdf$PCA_eigen1 <= 0,]

# Half-life vs RNA-seq
png("half-life_FPKM_density_scatterplot.png",
    width=2000, height=2000, res=300)
plot(mdf$hl, mdf$FPKM, pch=20, col=rgb(0,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="40kb mean(FPKM)"
     #ylim = c(-10,12), xlim= c(20,150)
)
text(200, 10, bquote(rho ~ " = "~ .(round(cor(mdf$hl, mdf$FPKM, method="spearman"), 2))))
dev.off()

# Half-life vs log(RNA-seq)
png("half-life_log_FPKM_density_scatterplot.png",
    width=2000, height=2000, res=300)
plot(mdf$hl, log2(mdf$FPKM), pch=20, col=rgb(0,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="40kb log2(mean(FPKM))",
      ylim=c(-10,12), xlim= c(20,150)
)
#abline(h=-4, lty="dashed", col="red")
text(120, 10, bquote(rho ~ " = "~ .(round(cor(mdf$hl, log2(mdf$FPKM), method="spearman"), 2))))
dev.off()

# A compartment half-life vs log(RNA-seq)
png("half-life_log_FPKM_density_scatterplot_A.png",
    width=2000, height=2000, res=300)
plot(A$hl, log2(A$FPKM), pch=20, col=rgb(1,0,0,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="40kb log2(mean(FPKM))",
     ylim=c(-10,12), xlim= c(20,150)
)
text(120, 10, bquote(rho ~ " = "~ .(round(cor(A$hl, log2(A$FPKM), method="spearman"), 2))))
dev.off()

# B compartment half-life vs log(RNA-seq)
png("half-life_log_FPKM_density_scatterplot_B.png",
    width=2000, height=2000, res=300)
plot(B$hl, log2(B$FPKM), pch=20, col=rgb(0,0,1,0.1),
     xlab=bquote("t"[1/2] ~ "(minutes)"), ylab="40kb log2(mean(FPKM))",
     ylim=c(-10,12), xlim= c(20,150)
)
text(120, 10, bquote(rho ~ " = "~ .(round(cor(B$hl, log2(B$FPKM), method="spearman"), 2))))
dev.off()

# Boxplots
highRNA_A <- A[A$FPKM > 1,]
lowRNA_A <- A[A$FPKM <= 1,]
highRNA_B <- B[B$FPKM > 1,]
lowRNA_B <- B[B$FPKM <= 1,]

png('boxplot_FPKM.png', width=1000, height=2000, res=300)
boxplot(highRNA_A$hl, lowRNA_A$hl, highRNA_B$hl, lowRNA_B$hl, 
        col=c("red", "red", "blue", "blue"), outline=FALSE,
        names=c("FPKM > 1", "FPKM <= 1", "FPKM > 1", "FPKM <= 1"),
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


pdf('boxplot_FPKM.pdf', width=3.5, height=8)
boxplot(highRNA_A$hl, lowRNA_A$hl, highRNA_B$hl, lowRNA_B$hl, 
        col=c("red", "red", "blue", "blue"), outline=FALSE,
        names=c("FPKM > 1", "FPKM <= 1", "FPKM > 1", "FPKM <= 1"),
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

# Stacked barplot for highRNA_A

hl <- factor(c(rep("40-45",4),
        rep("45-50",4), rep("50-55",4), rep("55-60",4),
        rep("60-65",4), rep("65-70",4), rep("70-75",4),
        rep("75-80",4), rep("80-85",4), rep("85-90",4),
        rep("90-95",4), rep("95-100",4),rep("100-105",4),
        rep("105-110",4), rep("110-115",4), rep("115-120",4)),
        levels=c("40-45",
                   "45-50", "50-55", "55-60",
                   "60-65", "65-70", "70-75",
                   "75-80", "80-85", "85-90",
                   "90-95", "95-100", "100-105",
                   "105-110", "110-115", "115-120"), 
        ordered=TRUE)

FPKM <- factor(rep(c("1-10", "10-50", "50-100", "100-1000"),16),
               levels=c("1-10", "10-50", "50-100", "100-1000"),
               ordered=TRUE)

hl_vec <- seq(40, 120,5)
FPKM_vec <- c(1,10,50,100,1000)

counts <- c()

for (i in 1:(length(hl_vec)-1)) {
  for (j in 1:(length(FPKM_vec)-1)) {
    counts <- c(counts, sum(highRNA_A$hl >= hl_vec[i] & highRNA_A$hl < hl_vec[i+1] &
                   highRNA_A$FPKM >= FPKM_vec[j] & highRNA_A$FPKM < FPKM_vec[j+1]))
  }
}

stack_df <- data.frame(hl, FPKM, counts)

pdf("FPKM_greaterthan_1_A_compartment_stacked_barplot.pdf",
    height=7,width=8)

ggplot(stack_df, aes(fill=FPKM, y=counts, x=hl)) + 
  geom_bar( stat="identity", position="fill") +
  xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab("Percent of Data") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()
  

# FPKM as x-axis for stacked barplot
hl <- factor(rep(c("40-60", "60-80", "80-100", "100-120"),4),
             levels=c("40-60", "60-80", "80-100", "100-120"),
             ordered=TRUE)

FPKM <- factor(c(rep("1-10",4), rep("10-50",4),
                 rep("50-100",4), rep("100-1000",4)),
              levels=c("1-10", "10-50", "50-100", "100-1000"),
              ordered=TRUE)
              
hl_vec <- seq(40, 120,20)
FPKM_vec <- c(1,10,50,100,1000)
counts <- c()

for (i in 1:(length(FPKM_vec)-1)) {
  for (j in 1:(length(hl_vec)-1)) {
    counts <- c(counts, sum(highRNA_A$hl >= hl_vec[i] & highRNA_A$hl < hl_vec[i+1] &
                              highRNA_A$FPKM >= FPKM_vec[j] & highRNA_A$FPKM < FPKM_vec[j+1]))
  }
}

stack_df <- data.frame(FPKM, hl, counts)
pdf("FPKM_greaterthan_1_A_compartment_stacked_barplot_FPKMxaxis.pdf",
    height=7,width=5)

ggplot(stack_df, aes(fill=hl, y=counts, x=FPKM)) + 
  geom_bar( stat="identity", position="fill") +
  xlab("FPKM") + ylab("Percent of Data") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5))
dev.off()
