#!/usr/bin/env Rscript
library(argparse)

parser <- ArgumentParser(description= "create interaction matrix from average.R output and plot")
parser$add_argument("-i", help="output from average.R", type="character", required=TRUE)
parser$add_argument("-o", help="prefix for matrix file and plot file", type="character", required=TRUE)

args <- parser$parse_args()

H3K9me3_df <- read.table(args$i, sep = "\t", header=TRUE)[,c(3,4,5)]

X = data.frame(matrix(0,nrow = length(seq(0,max(H3K9me3_df$H3K9me3_I), by=1000)), ncol = length(seq(0,max(H3K9me3_df$H3K9me3_J), by=1000))))
rownames(X) <- seq(0, max(H3K9me3_df$H3K9me3_I), by=1000)
colnames(X) <- seq(0, max(H3K9me3_df$H3K9me3_J), by=1000)

# Fill matrix
for (i in 1:nrow(H3K9me3_df)) {
  if (is.na(H3K9me3_df[i,"cScore"])){
    next
  }
  else {
    X[as.character(H3K9me3_df[i,"H3K9me3_I"]), as.character(H3K9me3_df[i,"H3K9me3_J"])] <- H3K9me3_df[i, "cScore"]
  }
}

colnames(X) <- paste(seq(0, max(H3K9me3_df$H3K9me3_J), by=1000), "|hg19|chr1:1-100", sep="")
rownames(X) <- paste(seq(0, max(H3K9me3_df$H3K9me3_I), by=1000), "|hg19|chr1:1-100", sep="")

#Only plot first 80 bins to handle outlier data
OUT  = X[1:80, 1:80]

write.table(OUT, paste(args$o, ".matrix", sep=""), col.names=NA, quote=FALSE, sep="\t", row.names=TRUE)

library(gplots)
png(paste(args$o, ".png", sep=""), height=4000, width=3500, res=300)
my_palette <- colorRampPalette(c("white", "red", "black"))(n = 1000)
heatmap.2(as.matrix(OUT), trace="none", col=my_palette, Rowv=FALSE, dendrogram = "none", Colv="Rowv")
dev.off()

          
          
          
          