# half-lifes are fun!!!!
source("C:/cygwin64/home/Tyler/Research/digest/digest-Hi-C/my_functions.R")

options(scipen=100)
#Get data
min5 <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-500000-raw_scaleBy_6.68.balanced_scaleBy_53.35.LOS.bedGraph", sep="\t", header=TRUE)
hour1 <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-500000-raw_scaleBy_8.9.balanced_scaleBy_55.21.LOS.bedGraph", sep="\t", header=TRUE)
hour2 <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-500000-raw_scaleBy_4.56.balanced_scaleBy_55.65.LOS.bedGraph", sep="\t", header=TRUE)
hour3 <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-500000-raw_scaleBy_4.54.balanced_scaleBy_56.08.LOS.bedGraph", sep="\t", header=TRUE)
hour4 <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-500000-raw_scaleBy_4.12.balanced_scaleBy_56.03.LOS.bedGraph", sep="\t", header=TRUE)
ON <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-500000-raw_scaleBy_3.89.balanced_scaleBy_56.69.LOS.bedGraph", sep="\t", header=TRUE)

# merging data
colnames(min5)[4] <- "min5_LOS"
colnames(hour1)[4] <- "h1_LOS"
colnames(hour2)[4] <- "h2_LOS"
colnames(hour3)[4] <- "h3_LOS"
colnames(hour4)[4] <- "h4_LOS"
colnames(ON)[4] <- "OVN_LOS"

# Remove end column
min5 <- min5[,-3]
hour1 <- hour1[,-3]
hour2 <- hour2[,-3]
hour3 <- hour3[,-3]
hour4 <- hour4[,-3]
ON <- ON[,-3]

# Merging
m <- merge(min5, hour1, by=c("chrom", "start"), all=TRUE)
m <- merge(m, hour2, by=c("chrom", "start"), all=TRUE)
m <- merge(m, hour3, by=c("chrom", "start"), all=TRUE)
m <- merge(m, hour4, by=c("chrom", "start"), all=TRUE)
m <- merge(m, ON, by=c("chrom", "start"), all=TRUE)

m <- OrderChromStart(m)

minutes <- c(5, 60, 120, 180, 240, 960)

# Example half-life
png("plot_half_life/chr1_49.5Mb_hl.png", height=1800, width=2500, res=300)
  plot(minutes, m[100, 3:8], xlab = "Minutes of DpnII Digestion", 
       ylab= "Loss of Structure", type='o', pch=20, col="chartreuse4",
       main="Chr1:49.5Mb-50Mb",ylim=c(0,1))
dev.off()

# Genome wide at 500kb
b <- boxplot(m[3:8])
png("plot_half_life/genome_hl.png", height=1800, width=2500, res=300)
boxplot(m[3:8], ylim=c(0,1), at=minutes, boxwex=50, names=minutes, 
    xlab="Minutes of DpnII Digestion", outline=FALSE,
    ylab= "Loss of Structure", main= "Genome", col="chartreuse4")
lines(minutes, b$stats[3,], col="black")
dev.off()

end <- m$start + 500000
o <- cbind(m[1:2], end, m[3:8])
# Write LOS data to table
write.table(o, "LOS_500kb.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

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
b <- boxplot(A[3:8])
A_med <- b$stats[3,]
png("plot_half_life/A_compartment_hl.png", height=1800, width=2500, res=300)
boxplot(A[3:8], ylim=c(0,1), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "A compartment", col="red")
lines(minutes, A_med, col="black")
dev.off()

b <- boxplot(B[3:8])
B_med <- b$stats[3,]
png("plot_half_life/B_compartment_hl.png", height=1800, width=2500, res=300)
boxplot(B[3:8], ylim=c(0,1), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "B compartment", col="blue")
lines(minutes, B_med, col="black")
dev.off()

# Overplot
png("plot_half_life/A_B_compartment_hl.png", height=1800, width=2500, res=300)
boxplot(A[3:8], ylim=c(0,1), at=minutes, boxwex=50, names=minutes, 
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "A and B compartments", col="red")
lines(minutes, A_med, col="red")

boxplot(B[3:8], at=minutes, boxwex=50,
        col="blue", add=TRUE, outline=FALSE, names=FALSE)
lines(minutes, B_med, col="blue")
dev.off()






