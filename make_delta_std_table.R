# half-lifes are fun!!!!
source("C:/cygwin64/home/Tyler/Research/digest/digest-Hi-C/my_functions.R")

options(scipen=100)
#Get data
min5 <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-500000-iced_scaleBy_56.28_delta_std.bedGraph", sep="\t", header=FALSE, skip=1, col.names=c("chrom", "start", "end", "min5_std"))
hour1 <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.39_delta_std.bedGraph", sep="\t", header=FALSE, skip=1, col.names=c("chrom", "start", "end", "h1_std"))
hour2 <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_58.76_delta_std.bedGraph", sep="\t", header=FALSE, skip=1,col.names=c("chrom", "start", "end", "h2_std"))
hour3 <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.28_delta_std.bedGraph", sep="\t", header=FALSE, skip=1,col.names=c("chrom", "start", "end", "h3_std"))
hour4 <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.21_delta_std.bedGraph", sep="\t", header=FALSE, skip=1,col.names=c("chrom", "start", "end", "h4_std"))
ON <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-500000-iced_scaleBy_59.8_delta_std.bedGraph", sep="\t", header=FALSE, skip=1,col.names=c("chrom", "start", "end", "OVN_std"))

# # merging data
# colnames(min5)[4] <- "min5_LOS"
# colnames(hour1)[4] <- "h1_LOS"
# colnames(hour2)[4] <- "h2_LOS"
# colnames(hour3)[4] <- "h3_LOS"
# colnames(hour4)[4] <- "h4_LOS"
# colnames(ON)[4] <- "OVN_LOS"

# Merging
m <- merge(min5, hour1, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour2, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour3, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour4, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, ON, by=c("chrom", "start", "end"), all=TRUE)
m <- OrderChromStart(m)


minutes <- c(5, 60, 120, 180, 240, 960)


# Example half-life
png("plot_delta_std/chr1_49.5Mb_delta_std.png", height=1800, width=2500, res=300)
  plot(minutes, m[100, 4:9], xlab = "Minutes of DpnII Digestion", 
       ylab= expression(paste(Delta,sigma,sep="")), type='o', pch=20, col="chartreuse4",
       main="Chr1:49.5Mb-50Mb", ylim=c(-0.3, 0.9))
dev.off()

# Genome wide at 500kb
b <- boxplot(m[4:9])
png("plot_delta_std/genome_std.png", height=1800, width=2500, res=300)
boxplot(m[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes, 
    xlab="Minutes of DpnII Digestion", outline=FALSE,
    ylab= expression(paste(Delta,sigma,sep="")), main= "Genome", col="chartreuse4")
lines(minutes, b$stats[3,], col="black")
dev.off()

# Write LOS data to table
write.table(m, "delta_std_500kb.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Compartments ####################################################################################
eigenPath <- "C:/cygwin64/home/Tyler/Research/digest/cis_percent/compartment/HBHiC-K562-MN-Dp-1__hg19__genome__C-500000-raw_scaleBy_2.72.balanced_scaleBy_51.45__all.zScore.eigen1.sorted.bedGraph"
eigen_table <- read.table(eigenPath, sep="\t")

# Fix bp indexing
eigen_table[2] <- eigen_table[2] - 1
colnames(eigen_table) <- c("chrom", "start", "end", "eigen")

# Merge
me <- merge(m, eigen_table, by=c("chrom", "start"), all=TRUE)

# Plot boxplots for A and B
A <- me[me$eigen > 0,]
B <- me[me$eigen <= 0,]

# Remove NA Y chrom
A <- A[rowSums(is.na(A)) != ncol(A),]
B <- B[rowSums(is.na(B)) != ncol(B),]

# Compartments at 500kb
b <- boxplot(A[4:9])
A_med <- b$stats[3,]
png("plot_delta_std/A_compartment_delta_std.png", height=1800, width=2500, res=300)
boxplot(A[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= expression(paste(Delta,sigma,sep="")), main= "A compartment", col="red")
lines(minutes, A_med, col="black")
dev.off()

b <- boxplot(B[4:9])
B_med <- b$stats[3,]
png("plot_delta_std/B_compartment_delta_std.png", height=1800, width=2500, res=300)
boxplot(B[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= expression(paste(Delta,sigma,sep="")), main= "B compartment", col="blue")
lines(minutes, B_med, col="black")
dev.off()

# Overplot
png("plot_delta_std/A_B_compartment_delta_std.png", height=1800, width=2500, res=300)
boxplot(A[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= expression(paste(Delta,sigma,sep="")), main= "A and B compartments", col="red")
lines(minutes, A_med, col="red")

boxplot(B[4:9], at=minutes, boxwex=50,
        col="blue", add=TRUE, outline=FALSE, names=FALSE)
lines(minutes, B_med, col="blue")
dev.off()






