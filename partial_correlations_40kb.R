library(corrplot)
library(ggplot2)

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

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}
# Remove chrom, start, end
chr2 <- subset(chr2, select=-c(chrom, start, end))
# Remove outliers
chr2_no_outliers <- data.frame(apply(chr2, 2, remove_outliers_std))
# Remove bins with NA
chr2_clean <- na.omit(chr2_no_outliers)
pc1 <- chr2_clean$PCA_eigen1
for_regress <- subset(chr2_clean, select=-PCA_eigen1)
residual_df <- data.frame(apply(for_regress, 2, get_residuals, p=pc1))

c <- cor(residual_df, use="pairwise.complete.obs", method="spearman")
png('partial_correlate_feature_matrix_40kb_chr2.png',width=3500, height=3500, res=300)
corrplot(c, method="circle",type = "upper", tl.col = "black")
dev.off()

cordf <- data.frame(c)
# Drop first 8 rows of loss of structure metrics
cordf <- cordf[9:nrow(cordf),]
cor_hlLOS <- cordf$half.life_LOS
cor_hlstd <- cordf$half.life_std
cor_values <- c()

for (i in 1:length(cor_hlLOS)) {
  cor_values <- c(cor_values, cor_hlLOS[i], cor_hlstd[i])
}



single_labs <- rownames(cordf)
double_labs <- c()
for (i in 1:length(single_labs)) {
  double_labs <- c(double_labs, single_labs[i], single_labs[i])
}

feature_labs <- factor(double_labs, levels=rownames(cordf), ordered=TRUE)

bar_df <- data.frame(cor_values, feature_labs)
bar_df$metric <- rep(c("LOS", "Delta sigma"), 44)


png("partial_cor_delta_sigma_vs_los_chr2_40kb.png", height=3000, width=3000, res=300)
ggplot(bar_df, aes(x=feature_labs, y= cor_values,fill=metric)) +
  geom_bar(stat = "identity", position="dodge") +
  scale_fill_manual(values=rep(c("darkorange1", "chartreuse4"), 52)) +
  coord_flip() + labs(x = "", y = "Spearman's Correlation") +
  theme_minimal()


dev.off()

suz <- chr2_clean$SUZ12
hl <- chr2_clean$half.life_LOS


png('halflife_LOS_vs_SUZ12_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(suz, hl, pch=21, col="black", bg="dodgerblue",
     xlab="SUZ12 fold change over control", ylab=bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     ylim=c(80,320))
dev.off()

res_suz <- residual_df$SUZ12
res_hl <- residual_df$half.life_LOS

png('halflife_LOS_vs_SUZ12_residuals_chr2_40kb.png', width=2500, height=2500, res=300)
par(mar=c(c(5, 5, 4, 2) + 0.1))
plot(res_suz, res_hl, pch=21, col="black", bg="dodgerblue",
     xlab="Residuals for SUZ12", ylab=bquote("Residuals for t"[1/2]), cex.lab=1.5,
     ylim=c(-70, 120))
dev.off()


print("SUZ12 vs half-life correlation:")
print(cor(suz, hl, method="spearman", use="pairwise.complete.obs"))
print("Residual SUZ12 vs residual half-life correlation:")
print(cor(res_suz, res_hl, method="spearman", use="pairwise.complete.obs"))
