library(ggplot2)
eigen_df <- read.table("../../../eigen/eigen1_40kb.bedGraph", sep="\t", header=FALSE)
thalf_df <- read.table("../rangehalf-life/half-life_exponential_40kb_removed_outliers_range6Mb.bedGraph", sep="\t", header=FALSE)
m_cis_df <- read.table("HBHiC-K562-MN-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_range6Mb_cispercent.bedGraph", sep="\t", header=FALSE)

colnames(eigen_df) <- c("chrom", "start", "end", "eigen")
colnames(thalf_df) <- c("chrom", "start", "end", "thalf")
colnames(m_cis_df) <- c("chrom", "start", "end", "cis")
m_cis_df$start <- m_cis_df$start - 1 

m <- merge(eigen_df, thalf_df, by=c("chrom", "start", "end"))
m <- merge(m, m_cis_df, by=c("chrom", "start", "end"))


chr2 <- m[m$chrom == "chr2",]
chr2 <- na.omit(chr2)

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}

png('eigen_vs_half-life_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(chr2$eigen, chr2$thalf, pch=21, col="black", bg="dodgerblue",
     xlab="PC1", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     ylim=c(25,150))
rho <- sprintf("%.3f", round(cor(chr2$eigen, chr2$thalf, method="spearman"), 3))
text(0.015, 125, bquote(rho == .(rho) ), cex=2)
dev.off()

thalf <- chr2$thalf
cis <- chr2$cis
eigen <- chr2$eigen

thalf_residuals <- get_residuals(chr2$thalf, chr2$cis)
eigen_residuals <- get_residuals(chr2$eigen, chr2$cis)

thalfm <- lm(chr2$thalf~chr2$cis)
eigenm <- lm(chr2$eigen~chr2$cis)

png('mock_cis_vs_half-life_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(chr2$cis, chr2$thalf, pch=21, col="black", bg="dodgerblue",
     xlab="Mock 6Mb cis %", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     ylim=c(25,150))
rho <- sprintf("%.3f",round(cor(chr2$cis, chr2$thalf, method="spearman"), 3))
text(30, 125, bquote(rho == .(rho) ), cex=2)
abline(thalfm$coefficients[1], thalfm$coefficients[2], col="red")
dev.off()

png('mock_cis_vs_eigen_chr2.png', width=2500, height=3500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(chr2$cis, chr2$eigen, pch=21, col="black", bg="dodgerblue",
     xlab="Mock 6Mb cis %", ylab="PC1", cex.lab=1.5)
rho <- sprintf("%.3f",round(cor(chr2$cis, chr2$eigen, method="spearman"), 3))
text(30, 0.025,bquote(rho == .(rho) ), cex=2)
abline(eigenm$coefficients[1], eigenm$coefficients[2], col="red")
dev.off()

pdf('partial_correlation_eigen_vs_half-life_controlled_by_mock_cis_chr2.pdf', width=8, height=9)
par(mar=c(c(5, 5, 4, 2) + 0.1))
A_eigenR <- eigen_residuals[eigen > 0]
B_eigenR <- eigen_residuals[eigen < 0]
A_thalfR <- thalf_residuals[eigen > 0]
B_thalfR <- thalf_residuals[eigen < 0]

plot(A_eigenR, A_thalfR, pch=20, col="red",
     xlab="PC1 residuals", ylab=bquote("t"[1/2] ~ " residuals"), cex.lab=1.5, cex=0.5,
     ylim=c(-45,60), xlim=c(-0.03, 0.03)
)
points(B_eigenR, B_thalfR, pch=20, col="blue", cex=0.5)
rho <- sprintf("%.3f", round(cor(eigen_residuals, thalf_residuals, method="spearman"), 3))
text(0.02,45, bquote(rho == .(rho) ), cex=2)

dev.off()

