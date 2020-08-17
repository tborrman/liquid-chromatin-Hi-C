mn_df <- read.table("HBHiC-K562-MN-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_57.11_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
m5_df <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-500000-iced_scaleBy_56.28_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
h1_df <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.39_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
h2_df <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.76_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
h3_df <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.28_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
h4_df <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.21_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))
ovn_df <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.8_std.bedGraph", sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "std"))

options(scipen=100)
# Chrom boxplot
plot_boxplot <- function(df, timepoint) {
png(paste("std_chrom_boxplots/",timepoint, "_std_chrom_boxplots.png", sep=""), width= 3000, height = 2000, res = 300)

par(mar=c(5, 5, 4, 2) + 0.1)
boxplot(df[df$chrom == "chr1",]$std,
        df[df$chrom == "chr2",]$std,
        df[df$chrom == "chr3",]$std,
        df[df$chrom == "chr4",]$std,
        df[df$chrom == "chr5",]$std,
        df[df$chrom == "chr6",]$std,
        df[df$chrom == "chr7",]$std,
        df[df$chrom == "chr8",]$std,
        df[df$chrom == "chr9",]$std,
        df[df$chrom == "chr10",]$std,
        df[df$chrom == "chr11",]$std,
        df[df$chrom == "chr12",]$std,
        df[df$chrom == "chr13",]$std,
        df[df$chrom == "chr14",]$std,
        df[df$chrom == "chr15",]$std,
        df[df$chrom == "chr16",]$std,
        df[df$chrom == "chr17",]$std,
        df[df$chrom == "chr18",]$std,
        df[df$chrom == "chr19",]$std,
        df[df$chrom == "chr20",]$std,
        df[df$chrom == "chr21",]$std,
        df[df$chrom == "chr22",]$std,
        df[df$chrom == "chrX",]$std,range=0, col=rainbow(23), 
        names=c(1:22, 'X'), ylab=expression(sigma), cex.lab=1.5, main=timepoint, cex.main=1.5)
dev.off()
}

plot_boxplot(mn_df, "MN")
plot_boxplot(m5_df, "5min")
plot_boxplot(h1_df, "1h")
plot_boxplot(h2_df, "2h")
plot_boxplot(h3_df, "3h")
plot_boxplot(h4_df, "4h")
plot_boxplot(ovn_df, "OVN")

# plot std across timescales for chromosome 2
plot_chr2_std <- function(df, timepoint) {
png(paste("std_chr2/", timepoint, "_std_chr2.png",sep=""), width=8000, height=2000, res=400)
  par(mar=c(5,6,4,8) + 0.1, lwd=2)
  mid <- df$start + (abs(df$end - df$start)/2.0)
  chrom_mid <- mid[df$chrom == "chr2"]
  chrom_std <- df$std[df$chrom == "chr2"]
  plot(chrom_mid, chrom_std, type='l', col="dodgerblue", xlab= "chr2",
       ylab=expression(sigma), cex.lab = 2, ylim=c(0,3750), main= timepoint, 
       cex.main=2)
dev.off()
}

plot_chr2_std(mn_df, "MN")
plot_chr2_std(m5_df, "5min")
plot_chr2_std(h1_df, "1h")
plot_chr2_std(h2_df, "2h")
plot_chr2_std(h3_df, "3h")
plot_chr2_std(h4_df, "4h")
plot_chr2_std(ovn_df, "OVN")

# All in one
png("std_chr2/plot_std_chr2.png", width=8000, height=2000, res=400)
par(mar=c(5,6,4,8) + 0.1, lwd=2)
  # MN
  mid <- mn_df$start + (abs(mn_df$end - mn_df$start)/2.0)
  chrom_mid <- mid[mn_df$chrom == "chr2"]
  chrom_std <- mn_df$std[mn_df$chrom == "chr2"]
  plot(chrom_mid, chrom_std, type='l', col="orange", xlab= "chr2",
       ylab=expression(sigma), cex.lab = 2, ylim=c(0,3200))
  # 5min
  mid <- m5_df$start + (abs(m5_df$end - m5_df$start)/2.0)
  chrom_mid <- mid[m5_df$chrom == "chr2"]
  chrom_std <- m5_df$std[m5_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="blue")
  # 1h
  mid <- h1_df$start + (abs(h1_df$end - h1_df$start)/2.0)
  chrom_mid <- mid[h1_df$chrom == "chr2"]
  chrom_std <- h1_df$std[h1_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="purple")
  # 2h
  mid <- h2_df$start + (abs(h2_df$end - h2_df$start)/2.0)
  chrom_mid <- mid[h2_df$chrom == "chr2"]
  chrom_std <- h2_df$std[h2_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="red")
  # 3h
  mid <- h3_df$start + (abs(h3_df$end - h3_df$start)/2.0)
  chrom_mid <- mid[h3_df$chrom == "chr2"]
  chrom_std <- h3_df$std[h3_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="magenta")
  # 4h
  mid <- h4_df$start + (abs(h4_df$end - h4_df$start)/2.0)
  chrom_mid <- mid[h4_df$chrom == "chr2"]
  chrom_std <- h4_df$std[h4_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="cyan")
  # OVN
  mid <- ovn_df$start + (abs(ovn_df$end - ovn_df$start)/2.0)
  chrom_mid <- mid[ovn_df$chrom == "chr2"]
  chrom_std <- ovn_df$std[ovn_df$chrom == "chr2"]
  lines(chrom_mid, chrom_std, type='l', col="limegreen")
  par(xpd=TRUE)
  legend(max(chrom_mid) + 13000000,3000, legend=c("OVN", "4h", "3h", "2h", "1h", "5min", "MN"), lty=1,
         col=c("limegreen", "cyan", "magenta", "red", "purple", "blue", "orange"))
dev.off()

# Scatter standard deviation vs LOS genome wide per timepoint
los_min5 <- read.table("../LOS/HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-500000-iced_scaleBy_56.28_LOS.bedGraph", sep="\t", header=TRUE)
los_hour1 <- read.table("../LOS/HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.39_LOS.bedGraph", sep="\t", header=TRUE)
los_hour2 <- read.table("../LOS/HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.76_LOS.bedGraph", sep="\t", header=TRUE)
los_hour3 <- read.table("../LOS/HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.28_LOS.bedGraph", sep="\t", header=TRUE)
los_hour4 <- read.table("../LOS/HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.21_LOS.bedGraph", sep="\t", header=TRUE)
los_ON <- read.table("../LOS/HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.8_LOS.bedGraph", sep="\t", header=TRUE)

# fix bp indexing
m5_df[2] <- m5_df[2] - 1
h1_df[2] <- h1_df[2] - 1
h2_df[2] <- h2_df[2] - 1
h3_df[2] <- h3_df[2] - 1
h4_df[2] <- h4_df[2] - 1
ovn_df[2] <- ovn_df[2] - 1

# 5min 
m5m <- merge(m5_df, los_min5, c("chrom", "start"))
png("std_LOS_cor/5min_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
  par(mar=c(5, 5, 4, 2) + 0.1)
  plot(m5m$std[m5m$chrom == "chr2"], m5m$LOS[m5m$chrom == "chr2"], pch=21,col="black",
       bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
       cex.lab = 1.5, ylim=c(-0.3, 0.3), xlim=c(1000,3000), main="5min", 
       cex.main = 1.5)
  text(2500, 0.2, bquote("r"[s]~ "="~ 
             .(round(cor(m5m$std[m5m$chrom == "chr2"], m5m$LOS[m5m$chrom == "chr2"], 
                 method="spearman", use="complete.obs"), 2)))
             ,cex= 1.5)
dev.off()
# 1h
h1m <- merge(h1_df, los_hour1, c("chrom", "start"))
png("std_LOS_cor/1h_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(h1m$std[h1m$chrom == "chr2"], h1m$LOS[h1m$chrom == "chr2"], pch=21,col="black",
     bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
     cex.lab = 1.5, ylim=c(-0.1, 0.5), xlim=c(500,2500), main="1h", 
     cex.main = 1.5)
text(2000, 0.4, bquote("r"[s]~ "="~ 
                         .(round(cor(h1m$std[h1m$chrom == "chr2"], h1m$LOS[h1m$chrom == "chr2"], 
                                     method="spearman", use="complete.obs"), 2)))
     ,cex= 1.5)
dev.off()
# 2h
h2m <- merge(h2_df, los_hour2, c("chrom", "start"))
png("std_LOS_cor/2h_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(h2m$std[h2m$chrom == "chr2"], h2m$LOS[h2m$chrom == "chr2"], pch=21,col="black",
     bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
     cex.lab = 1.5, ylim=c(-0.1, 0.5), xlim=c(200,2500), main="2h", 
     cex.main = 1.5)
text(2000, 0.4, bquote("r"[s]~ "="~ 
                         .(round(cor(h2m$std[h2m$chrom == "chr2"], h2m$LOS[h2m$chrom == "chr2"], 
                                     method="spearman", use="complete.obs"), 2)))
     ,cex= 1.5)
dev.off()
# 3h
h3m <- merge(h3_df, los_hour3, c("chrom", "start"))
png("std_LOS_cor/3h_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(h3m$std[h3m$chrom == "chr2"], h3m$LOS[h3m$chrom == "chr2"], pch=21,col="black",
     bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
     cex.lab = 1.5, ylim=c(0.2, 0.6), xlim=c(200,2000), main="3h", 
     cex.main = 1.5)
text(1500, 0.5, bquote("r"[s]~ "="~ 
                         .(round(cor(h3m$std[h3m$chrom == "chr2"], h3m$LOS[h3m$chrom == "chr2"], 
                                     method="spearman", use="complete.obs"), 2)))
     ,cex= 1.5)
dev.off()
# 4h
h4m <- merge(h4_df, los_hour4, c("chrom", "start"))
png("std_LOS_cor/4h_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(h4m$std[h4m$chrom == "chr2"], h4m$LOS[h4m$chrom == "chr2"], pch=21,col="black",
     bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
     cex.lab = 1.5, ylim=c(0.2, 0.6), xlim=c(200,2000), main="4h", 
     cex.main = 1.5)
text(1500, 0.5, bquote("r"[s]~ "="~ 
                         .(round(cor(h4m$std[h4m$chrom == "chr2"], h4m$LOS[h4m$chrom == "chr2"], 
                                     method="spearman", use="complete.obs"), 2)))
     ,cex= 1.5)
dev.off()
# OVN
om <- merge(ovn_df, los_ON, c("chrom", "start"))
png("std_LOS_cor/ovn_std_LOS_cor_chr2.png", width=3000, height=3000, res=350)
par(mar=c(5, 5, 4, 2) + 0.1)
plot(om$std[om$chrom == "chr2"], om$LOS[om$chrom == "chr2"], pch=21,col="black",
     bg="dodgerblue", xlab=expression(sigma),ylab="Loss of Structure", 
     cex.lab = 1.5, ylim=c(0.5, 0.85), xlim=c(200,2000),main="OVN", 
     cex.main = 1.5)
text(1500, 0.75, bquote("r"[s]~ "="~ 
                         .(round(cor(om$std[om$chrom == "chr2"], om$LOS[om$chrom == "chr2"], 
                                     method="spearman", use="complete.obs"), 2)))
     ,cex= 1.5)
dev.off()


