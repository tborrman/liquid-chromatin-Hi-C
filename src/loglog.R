# wkdir -> "C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse2/LOS1hourtest/log-log"
options(scipen=100)
eigen <- read.table("../../../../eigen/eigen1_40kb.bedGraph", sep="\t", header=FALSE,
                    col.names=c("chrom", "start", "end", "PC1"))
#Get data
MN <- read.table("../cispercent6mb/HBHiC-K562-MN-5mDP1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_cispercent.bedGraph",
                 sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
min5 <- read.table("../cispercent6mb/HBHiCK562DN10-5mDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.41_range6Mb_cispercent.bedGraph",
                   sep="\t", header=FALSE,  col.names=c("chrom", "start", "end", "cis"))
hour1 <- read.table("../cispercent6mb/HBHiCK562DN10-24hDp1-filter100__hg19__genome__C-40000-iced_scaleBy_0.45_range6Mb_cispercent.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour2 <- read.table("../cispercent6mb/HBHiCK562DN10-1hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour3 <- read.table("../cispercent6mb/HBHiCK562DN10-2hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph", 
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour4 <- read.table("../cispercent6mb/HBHiCK562DN10-3hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour8 <- read.table("../cispercent6mb/HBHiCK562DN10-4hDp2-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour16 <- read.table("../cispercent6mb/HBHiCK562DN10-8hDp1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph", 
                     sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))
hour24 <- read.table("../cispercent6mb/HBHiCK562DN10-16hDp1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.42_range6Mb_cispercent.bedGraph", 
                     sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "cis"))

minutes <- c(0, 60, 120, 180, 240, 480, 960, 1440)

# A compartment##################################################################################################################
#chr2 170680001 170720000
A_cis <- c(MN[10500, "cis"], hour1[10500, "cis"],
         hour2[10500, "cis"], hour3[10500, "cis"], hour4[10500, "cis"],
         hour8[10500, "cis"], hour16[10500, "cis"], hour24[10500, "cis"])
png("cis_percent_A_bin.png", height=2000, width=2000, res=300)
plot(minutes, A_cis, pch=20, col="red", ylab="Cis %",
     xlab = "Minutes after DpnII digestion", type="o",
     main = "Chr2:170680001-170720000 (A bin)", ylim=c(0,40))
dev.off()
png("log-log_cis_percent_A_bin.png", height=2000, width=2000, res=300)
plot(log2(minutes), log2(A_cis), pch=20, col="red", ylab=bquote("log"[2]~"(Cis %)"),
     xlab = bquote("log"[2]~"(Minutes after DpnII digestion)"), type="o",
     main = "Chr2:170680001-170720000 (A bin)", ylim=c(0,5.5))
dev.off()

png("cis_percent_A_500bins.png", height=2000, width=2000, res=300)
plot(minutes, A_cis, col="red", ylab="Cis %",
     xlab = "Minutes after DpnII digestion", type="l",lwd=0.5,
     main = "500 A bins", ylim=c(0,40))
counter = 0
for (i in 1:10000) {
  if (counter > 500) {
    break
  }
  if (!is.na(eigen[i, "PC1"])) {
    if(eigen[i, "PC1"] > 0) {
      A_cis <- c(MN[i, "cis"], hour1[i, "cis"],
                 hour2[i, "cis"], hour3[i, "cis"], hour4[i, "cis"],
                 hour8[i, "cis"], hour16[i, "cis"], hour24[i, "cis"])
      lines(minutes, A_cis, col="red", lwd=0.5)
      counter = counter + 1
    }
  }
}
dev.off()

png("log-log_cis_percent_A_500bins.png", height=2000, width=2000, res=300)
plot(log2(minutes), log2(A_cis), col="red", ylab=bquote("log"[2]~"(Cis %)"),
     xlab = bquote("log"[2]~"(Minutes after DpnII digestion)"), type="l",lwd=0.5,
     main = "500 A bins", ylim=c(0,5.5))
counter = 0
for (i in 1:10000) {
  if (counter > 500) {
    break
  }
  if (!is.na(eigen[i, "PC1"])) {
    if(eigen[i, "PC1"] > 0) {
      A_cis <- c(MN[i, "cis"], hour1[i, "cis"],
                 hour2[i, "cis"], hour3[i, "cis"], hour4[i, "cis"],
                 hour8[i, "cis"], hour16[i, "cis"], hour24[i, "cis"])
      lines(log2(minutes), log2(A_cis), col="red", lwd=0.5)
      counter = counter + 1
    }
  }
}
dev.off()
# B compartment##################################################################################################################
# chr2:51080000-51120000 -0.01685 (bin 7510)
B_cis <- c(MN[7510, "cis"], hour1[7510, "cis"],
           hour2[7510, "cis"], hour3[7510, "cis"], hour4[7510, "cis"],
           hour8[7510, "cis"], hour16[7510, "cis"], hour24[7510, "cis"])
png("cis_percent_B_bin.png", height=2000, width=2000, res=300)
plot(minutes, B_cis, pch=20, col="blue", ylab="Cis %",
     xlab = "Minutes after DpnII digestion", type="o",
     main = "Chr2:51080000-51120000 (B bin)", ylim=c(0,40))
dev.off()
png("log-log_cis_percent_B_bin.png", height=2000, width=2000, res=300)
plot(log2(minutes), log2(B_cis), pch=20, col="blue", ylab=bquote("log"[2]~"(Cis %)"),
     xlab = bquote("log"[2]~"(Minutes after DpnII digestion)"), type="o",
     main = "Chr2:51080000-51120000 (B bin)", ylim=c(0,5.5))
dev.off()

png("cis_percent_B_500bins.png", height=2000, width=2000, res=300)
plot(minutes, B_cis, col="blue", ylab="Cis %",
     xlab = "Minutes after DpnII digestion", type="l",lwd=0.5,
     main = "500 B bins", ylim=c(0,40))
counter = 0
for (i in 1:10000) {
  if (counter > 500) {
    break
  }
  if (!is.na(eigen[i, "PC1"])) {
    if(eigen[i, "PC1"] < 0) {
      B_cis <- c(MN[i, "cis"], hour1[i, "cis"],
                 hour2[i, "cis"], hour3[i, "cis"], hour4[i, "cis"],
                 hour8[i, "cis"], hour16[i, "cis"], hour24[i, "cis"])
      lines(minutes, B_cis, col="blue", lwd=0.5)
      counter = counter + 1
    }
  }
}
dev.off()

png("log-log_cis_percent_B_500bins.png", height=2000, width=2000, res=300)
plot(log2(minutes), log2(B_cis), col="blue", ylab=bquote("log"[2]~"(Cis %)"),
     xlab = bquote("log"[2]~"(Minutes after DpnII digestion)"), type="l",lwd=0.5,
     main = "500 B bins", ylim=c(0,5.5))
counter = 0
for (i in 1:10000) {
  if (counter > 500) {
    break
  }
  if (!is.na(eigen[i, "PC1"])) {
    if(eigen[i, "PC1"] < 0) {
      B_cis <- c(MN[i, "cis"], hour1[i, "cis"],
                 hour2[i, "cis"], hour3[i, "cis"], hour4[i, "cis"],
                 hour8[i, "cis"], hour16[i, "cis"], hour24[i, "cis"])
      lines(log2(minutes), log2(B_cis), col="blue", lwd=0.5)
      counter = counter + 1
    }
  }
}
dev.off()



 
