df <- read.table("/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/feature_matrix_40kb.txt", sep="\t", header=TRUE)
chr2 <- df[df$chrom == "chr2",]

remove_outliers_std <- function(x) {
  sigma <- sd(x, na.rm=TRUE)
  mu <- mean(x, na.rm=TRUE)
  lower_bound <- mu - (sigma*3)
  upper_bound <- mu + (sigma*3)
  x[(x > upper_bound) | (x < lower_bound)] <- NA
  return(x)
}

remove_outliers_IQR <- function(x) {
  i <- IQR(x, na.rm=TRUE)
  mu <- mean(x, na.rm=TRUE)
  lower_bound <- mu - (i*1.5)
  upper_bound <- mu + (i*1.5)
  x[(x > upper_bound) | (x < lower_bound)] <- NA
  return(x)
}

# Remove outliers
hl <- remove_outliers_IQR(chr2$half.life_LOS)
pc1 <- remove_outliers_IQR(chr2$PCA_eigen1)
lad <- remove_outliers_IQR(chr2$LAD_clone14)

# Remove bins with NA in lad, half-life or PC1
na_df <- data.frame(hl, lad, pc1)
complete_df <- na.omit(na_df)

hl <- complete_df$hl
lad <- complete_df$lad
pc1 <- complete_df$pc1

png('halflife_LOS_vs_PC1_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(pc1, hl, pch=21, col="black", bg="dodgerblue",
     xlab="PC1", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5)
dev.off()

png('LAD_vs_PC1_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(pc1, lad, pch=21, col="black", bg="dodgerblue",
     xlab="PC1", ylab="LAD signal", cex.lab=1.5)
dev.off()

png('LAD_vs_halflife_LOS_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(lad, hl, pch=21, col="black", bg="dodgerblue",
      xlab="LAD_signal", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5)
dev.off()


lm_hl <- lm(hl~pc1)
res_hl <- lm_hl$residuals

png('halflife_LOS_vs_PC1_regress_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(pc1, hl, pch=21, col="black", bg="dodgerblue",
     xlab="PC1", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5)
abline(a=lm_hl$coefficients[1], b=lm_hl$coefficients[2], col="red", lwd=2)
dev.off()


lm_lad <- lm(lad~pc1)
res_lad <- lm_lad$residuals

png('LAD_vs_PC1_regress_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(pc1, lad, pch=21, col="black", bg="dodgerblue",
     xlab="PC1", ylab="LAD signal", cex.lab=1.5)
abline(a=lm_lad$coefficients[1], b=lm_lad$coefficients[2], col="red", lwd=2)
dev.off()

png('half-life_vs_LAD_residuals_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(res_lad, res_hl, pch=21, col="black", bg="dodgerblue",
     xlab="Residuals for LAD signal", ylab=bquote("Residuals for t"[1/2]), cex.lab=1.5)
dev.off()

print("LAD vs half-life correlation:")
print(cor(lad, hl, method="spearman", use="pairwise.complete.obs"))
print("Residual LAD vs residual half-life correlation:")
print(cor(res_lad, res_hl, method="spearman", use="pairwise.complete.obs"))




