# plot LOS metric for each chromosome at each time point
options(scipen=100)
#Get data
min5 <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour1 <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour2 <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour3 <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour4 <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
ON <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)

for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {

  
  png(paste("plot_LOS/plot_LOS_", chrom, ".png",sep=""), width=8000, height=2000, res=400)
  par(mar=c(5,6,4,8) + 0.1, lwd=2)
    #min5
    mid <- min5$start + (abs(min5$end - min5$start)/2.0)
    chrom_mid <- mid[min5$chrom == chrom]
    chrom_LOS <- min5$LOS[min5$chrom == chrom]
    plot(chrom_mid, chrom_LOS, type='l', col="blue", xlab= chrom,
         ylab="Loss of Structure", cex.lab = 2, ylim=c(-0.30, 0.85))
    #hour1
    mid <- hour1$start + (abs(hour1$end - hour1$start)/2.0)
    chrom_mid <- mid[hour1$chrom == chrom]
    chrom_LOS <- hour1$LOS[hour1$chrom == chrom]
    lines(chrom_mid, chrom_LOS, type='l', col="purple")
    #hour2
    mid <- hour2$start + (abs(hour2$end - hour2$start)/2.0)
    chrom_mid <- mid[hour2$chrom == chrom]
    chrom_LOS <- hour2$LOS[hour2$chrom == chrom]
    lines(chrom_mid, chrom_LOS, type='l', col="red")
    #hour3
    mid <- hour3$start + (abs(hour3$end - hour3$start)/2.0)
    chrom_mid <- mid[hour3$chrom == chrom]
    chrom_LOS <- hour3$LOS[hour3$chrom == chrom]
    lines(chrom_mid, chrom_LOS, type='l', col="magenta")
    #hour4
    mid <- hour4$start + (abs(hour4$end - hour4$start)/2.0)
    chrom_mid <- mid[hour4$chrom == chrom]
    chrom_LOS <- hour4$LOS[hour4$chrom == chrom]
    lines(chrom_mid, chrom_LOS, type='l', col="cyan")
    #ON
    mid <- ON$start + (abs(ON$end - ON$start)/2.0)
    chrom_mid <- mid[ON$chrom == chrom]
    chrom_LOS <- ON$LOS[ON$chrom == chrom]
    lines(chrom_mid, chrom_LOS, type='l', col="limegreen")
  par(xpd=TRUE)
  legend(max(chrom_mid) + 13000000, 0.85, legend=c("OVN", "4h", "3h", "2h", "1h", "5min"), lty=1,
         col=c("limegreen", "cyan", "magenta", "red", "purple", "blue"))
  dev.off()
  
  pdf(paste("plot_LOS/plot_LOS_", chrom, ".pdf",sep=""), width=22, height=5)
  par(mar=c(5,6,4,8) + 0.1, lwd=2)
  #min5
  mid <- min5$start + (abs(min5$end - min5$start)/2.0)
  chrom_mid <- mid[min5$chrom == chrom]
  chrom_LOS <- min5$LOS[min5$chrom == chrom]
  plot(chrom_mid, chrom_LOS, type='l', col="blue", xlab= chrom,
       ylab="Loss of Structure", cex.lab = 2, ylim=c(-0.30, 0.85))
  #hour1
  mid <- hour1$start + (abs(hour1$end - hour1$start)/2.0)
  chrom_mid <- mid[hour1$chrom == chrom]
  chrom_LOS <- hour1$LOS[hour1$chrom == chrom]
  lines(chrom_mid, chrom_LOS, type='l', col="purple")
  #hour2
  mid <- hour2$start + (abs(hour2$end - hour2$start)/2.0)
  chrom_mid <- mid[hour2$chrom == chrom]
  chrom_LOS <- hour2$LOS[hour2$chrom == chrom]
  lines(chrom_mid, chrom_LOS, type='l', col="red")
  #hour3
  mid <- hour3$start + (abs(hour3$end - hour3$start)/2.0)
  chrom_mid <- mid[hour3$chrom == chrom]
  chrom_LOS <- hour3$LOS[hour3$chrom == chrom]
  lines(chrom_mid, chrom_LOS, type='l', col="magenta")
  #hour4
  mid <- hour4$start + (abs(hour4$end - hour4$start)/2.0)
  chrom_mid <- mid[hour4$chrom == chrom]
  chrom_LOS <- hour4$LOS[hour4$chrom == chrom]
  lines(chrom_mid, chrom_LOS, type='l', col="cyan")
  #ON
  mid <- ON$start + (abs(ON$end - ON$start)/2.0)
  chrom_mid <- mid[ON$chrom == chrom]
  chrom_LOS <- ON$LOS[ON$chrom == chrom]
  lines(chrom_mid, chrom_LOS, type='l', col="limegreen")
  par(xpd=TRUE)
  legend(max(chrom_mid) + 10000000, 0.85, legend=c("OVN", "4h", "3h", "2h", "1h", "5min"), lty=1,
         col=c("limegreen", "cyan", "magenta", "red", "purple", "blue"))
  dev.off()
  
  
}


