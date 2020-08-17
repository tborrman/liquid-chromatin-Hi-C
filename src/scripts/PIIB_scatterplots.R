eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen[is.na(eigen$PC1),"PC1"] <- NA
LOS_DpnII <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/LOS6Mb/HBHiCK562DN10-4h-DpnII-R1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_d"))
LOS_DpnII_TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBTD-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.59_range6Mb_LOS.bedGraph",
                           sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_TD"))
LOS_DpnII_NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBNT-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.65_range6Mb_LOS.bedGraph",
                           sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_NT"))


eigen_chr2 <- eigen[eigen$chrom == "chr2",]
LOS_DpnII_chr2 <- LOS_DpnII[LOS_DpnII$chrom == "chr2",]
LOS_DpnII_TD_chr2 <- LOS_DpnII_TD[LOS_DpnII_TD$chrom == "chr2",]
LOS_DpnII_NT_chr2 <- LOS_DpnII_NT[LOS_DpnII_NT$chrom == "chr2",]


d <- data.frame(eigen_chr2[1:4], LOS_DpnII_chr2["LOS_d"], LOS_DpnII_TD_chr2["LOS_TD"],
                LOS_DpnII_NT_chr2["LOS_NT"])

d_genome <- data.frame(eigen[1:4], LOS_DpnII["LOS_d"], LOS_DpnII_TD["LOS_TD"], LOS_DpnII_NT["LOS_NT"])
# NA bins
# d$PC1[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
# d$LOS_d[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
# d$signal[which(is.na(d), arr.ind = TRUE)[,1]] <- NA
# d_clean <- na.omit(d)

pdf("PIIB_scatterplot_chr2.pdf", height=5, width=10)
par(mfrow=c(1,3), mar=c(5, 5, 3, 2) + 0.1)

plot(d$PC1, d$LOS_d, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Timecourse R1\n4 hour digest",
     ylim=c(0.65,1))
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.69, bquote(rho == .(rho) ), cex=2)

plot(d$PC1, d$LOS_NT, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Not Treated\n4 hour digest",
     ylim=c(0.65,1))
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_NT, method="spearman", use="complete.obs"), 3))
text(0.01, 0.69, bquote(rho == .(rho) ), cex=2)

plot(d$PC1, d$LOS_TD, col=ifelse(d$PC1 > 0, "red", "blue"),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Pol II block\n4 hour digest",
     ylim=c(0.65,1))
rho <- sprintf("%.2f",round(cor(d$PC1, d$LOS_TD, method="spearman", use="complete.obs"), 3))
text(0.01, 0.69, bquote(rho == .(rho) ), cex=2)
dev.off()

png("PIIB_scatterplot_genome.png", height=1000, width=2200, res=300)
par(mfrow=c(1,3), mar=c(5, 5, 3, 2) + 0.1)

plot(d_genome$PC1, d_genome$LOS_d, col=ifelse(d_genome$PC1 > 0, rgb(1,0,0,0.01), rgb(0,0,1,0.01)),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Timecourse R1\n4 hour digest",
     ylim=c(0.60,1))
rho <- sprintf("%.2f",round(cor(d_genome$PC1, d_genome$LOS_d, method="spearman", use="complete.obs"), 3))
text(0.01, 0.65, bquote(rho == .(rho) ), cex=2)

plot(d_genome$PC1, d_genome$LOS_NT, col=ifelse(d_genome$PC1 > 0, rgb(1,0,0,0.01), rgb(0,0,1,0.01)),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Not Treated\n4 hour digest",
     ylim=c(0.60,1))
rho <- sprintf("%.2f",round(cor(d_genome$PC1, d_genome$LOS_NT, method="spearman", use="complete.obs"), 3))
text(0.01, 0.65, bquote(rho == .(rho) ), cex=2)

plot(d_genome$PC1, d_genome$LOS_TD, col=ifelse(d_genome$PC1 > 0, rgb(1,0,0,0.01), rgb(0,0,1,0.01)),
     pch=20, xlab="PC1", ylab="LOS DpnII (%)", cex=0.5, cex.lab=1.5, main="Pol II block\n4 hour digest",
     ylim=c(0.60,1))
rho <- sprintf("%.2f",round(cor(d_genome$PC1, d_genome$LOS_TD, method="spearman", use="complete.obs"), 3))
text(0.01, 0.65, bquote(rho == .(rho) ), cex=2)
dev.off()

