son1 <- read.table("GSE81553_SON_TSA-Seq_Condition1_40kb.bedGraph", sep="\t", header=FALSE,
                   col.names=c("chrom", "start", "end", "tsa1"))
son2 <- read.table("GSE81553_SON_TSA-Seq_Condition2_40kb.bedGraph", sep="\t", header=FALSE,
                   col.names=c("chrom", "start", "end", "tsa2"))
dist1 <- read.table("GSE81553_Speckle_distance_Condition1_40kb.bedGraph", sep="\t", header=FALSE,
                   col.names=c("chrom", "start", "end", "distance1"))
dist2 <- read.table("GSE81553_Speckle_distance_Condition2_40kb.bedGraph", sep="\t", header=FALSE,
                    col.names=c("chrom", "start", "end", "distance2"))
h <- read.table("/cygwin64/home/Tyler/Research/digest/LOS_6Mb/C-40000/rangehalf-life_LOS/half-life_exponential_40kb_removed_outliers_range6Mb.bedGraph",
                 sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "hl"))

df <- data.frame(h$hl, son1$tsa1, son2$tsa2, dist1$distance1, dist2$distance2)

upper.panel <- function(x, y) {
  points(x,y, pch=20, col=rgb(red=0, green=0, blue=1, alpha=0.1))
  r <- round(cor(x,y, method="spearman", use="complete.obs"), digits=2)
  txt <- paste("R = ", r)
  usr <- par("usr"); on.exit(par(usr))
  par(usr = c(0, 1, 0, 1))
  text(0.5, 0.9, txt)
}


png("TSA-seq_compare.png", width=2500, height=2500, res=300)
pairs(df,
      labels=c(bquote(atop("t"[1/2], "(minutes)")), expression("SON TSA-seq\nCondition 1"),
      expression("SON TSA-seq\nCondition 2"), expression("Speckle Distance\nCondition 1"),
      expression("Speckle Distance\nCondition 2")), 
      lower.panel = NULL, upper.panel = upper.panel, cex.labels=1.5)
dev.off()
