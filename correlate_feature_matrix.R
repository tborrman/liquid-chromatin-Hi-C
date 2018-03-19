library(corrplot)
df <- read.table("feature_matrix.txt", sep="\t", header=TRUE)
features <- df[,4:ncol(df)]
c <- cor(features, use="pairwise.complete.obs", method="spearman")
png('correlate_feature_matrix.png',width=2500, height=2500, res=300)
  corrplot(c, method="circle",type = "upper", tl.col = "black")
dev.off()

png('correlate_feature_matrix_values.png',width=2500, height=2500, res=300)
corrplot(c, method="number", type="upper", tl.col = "black", number.cex = 0.5)
dev.off()
