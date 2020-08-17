# plot delta_std metric for each chromosome at each time point
options(scipen=100)
#Get data
min5 <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-500000-iced_scaleBy_56.28_delta_std.bedGraph", sep="\t", header=TRUE)
hour1 <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.39_delta_std.bedGraph", sep="\t", header=TRUE)
hour2 <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.76_delta_std.bedGraph", sep="\t", header=TRUE)
hour3 <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.28_delta_std.bedGraph", sep="\t", header=TRUE)
hour4 <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.21_delta_std.bedGraph", sep="\t", header=TRUE)
ON <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.8_delta_std.bedGraph", sep="\t", header=TRUE)


for (chrom in c(paste("chr", 1:22, sep=""), "chrX")) {

  
  png(paste("plot_delta_std/plot_delta_std_", chrom, ".png",sep=""), width=8000, height=2000, res=400)
  par(mar=c(5,6,4,8) + 0.1, lwd=2)
    #min5
    mid <- min5$start + (abs(min5$end - min5$start)/2.0)
    chrom_mid <- mid[min5$chrom == chrom]
    chrom_delta_std <- min5$delta_std[min5$chrom == chrom]
    plot(chrom_mid, chrom_delta_std, type='l', col="blue", xlab= chrom,
         ylab=expression(paste(Delta,sigma, sep="")), cex.lab = 2, ylim=c(-0.4, 0.8))
    #hour1
    mid <- hour1$start + (abs(hour1$end - hour1$start)/2.0)
    chrom_mid <- mid[hour1$chrom == chrom]
    chrom_delta_std <- hour1$delta_std[hour1$chrom == chrom]
    lines(chrom_mid, chrom_delta_std, type='l', col="purple")
    #hour2
    mid <- hour2$start + (abs(hour2$end - hour2$start)/2.0)
    chrom_mid <- mid[hour2$chrom == chrom]
    chrom_delta_std <- hour2$delta_std[hour2$chrom == chrom]
    lines(chrom_mid, chrom_delta_std, type='l', col="red")
    #hour3
    mid <- hour3$start + (abs(hour3$end - hour3$start)/2.0)
    chrom_mid <- mid[hour3$chrom == chrom]
    chrom_delta_std <- hour3$delta_std[hour3$chrom == chrom]
    lines(chrom_mid, chrom_delta_std, type='l', col="magenta")
    #hour4
    mid <- hour4$start + (abs(hour4$end - hour4$start)/2.0)
    chrom_mid <- mid[hour4$chrom == chrom]
    chrom_delta_std <- hour4$delta_std[hour4$chrom == chrom]
    lines(chrom_mid, chrom_delta_std, type='l', col="cyan")
    #ON
    mid <- ON$start + (abs(ON$end - ON$start)/2.0)
    chrom_mid <- mid[ON$chrom == chrom]
    chrom_delta_std <- ON$delta_std[ON$chrom == chrom]
    lines(chrom_mid, chrom_delta_std, type='l', col="limegreen")
  par(xpd=TRUE)
  legend(max(chrom_mid) + 13000000, 0.85, legend=c("OVN", "4h", "3h", "2h", "1h", "5min"), lty=1,
         col=c("limegreen", "cyan", "magenta", "red", "purple", "blue"))
  dev.off()
}


