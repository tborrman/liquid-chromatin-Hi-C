
mylabels <- c("control", "hindIII", "dpnII")
trans_pct <- c(52, 52, 78)


pdf("trans_percent.pdf", height=1.495*5, width=0.822*5)
par(mar=c(2,2,2,2))
  barplot(trans_pct, names=mylabels, col=c("blue", "red", "cyan"), ylim=c(0,100))
dev.off()