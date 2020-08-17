# half-lifes are fun!!!!
source("C:/cygwin64/home/Tyler/Research/digest/digest-Hi-C/my_functions.R")

options(scipen=100)
#Get data
min5 <- read.table("HBHiCK562DN10-5m-DPnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour1 <- read.table("HBHiCK562DN10-1h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour2 <- read.table("HBHiCK562DN10-2h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour3 <- read.table("HBHiCK562DN10-3h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
hour4 <- read.table("HBHiCK562DN10-4h-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)
ON <- read.table("HBHiCK562DN10-ON-DpnII-R1__hg19__genome__C-40000-iced_scaleBy_0.39_LOS.bedGraph", sep="\t", header=TRUE)

# merging data
colnames(min5)[4] <- "min5_LOS"
colnames(hour1)[4] <- "h1_LOS"
colnames(hour2)[4] <- "h2_LOS"
colnames(hour3)[4] <- "h3_LOS"
colnames(hour4)[4] <- "h4_LOS"
colnames(ON)[4] <- "OVN_LOS"

# Merging
m <- merge(min5, hour1, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour2, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour3, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, hour4, by=c("chrom", "start", "end"), all=TRUE)
m <- merge(m, ON, by=c("chrom", "start", "end"), all=TRUE)

m <- OrderChromStart(m)


minutes <- c(5, 60, 120, 180, 240, 960)

# Example half-life
png("plot_LOS/chr1_49.96Mb_hl.png", height=1800, width=2500, res=300)
  plot(minutes, m[1250, 4:9], xlab = "Minutes of DpnII Digestion", 
       ylab= "Loss of Structure", type='o', pch=20, col="chartreuse4",
       main="Chr1:49.96Mb-50Mb",ylim=c(-0.3,0.9))
dev.off()

# Genome wide at 40kb
b <- boxplot(m[4:9])
png("plot_LOS/genome_LOS.png", height=1800, width=2500, res=300)
boxplot(m[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes, 
    xlab="Minutes of DpnII Digestion", outline=FALSE,
    ylab= "Loss of Structure", main= "Genome", col="chartreuse4")
lines(minutes, b$stats[3,], col="black")
dev.off()

# Write LOS data to table
write.table(m, "LOS_40kb.txt", sep="\t", col.names=TRUE, row.names=FALSE, quote=FALSE)

# Compartments ####################################################################################
eigenPath <- "C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph"
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

# Compartments at 40kb
b <- boxplot(A[4:9])
A_med <- b$stats[3,]
png("plot_LOS/A_compartment_LOS.png", height=1800, width=2500, res=300)
boxplot(A[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes,
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "A compartment", col="red")
lines(minutes, A_med, col="black")
dev.off()

b <- boxplot(B[4:9])
B_med <- b$stats[3,]
png("plot_LOS/B_compartment_LOS.png", height=1800, width=2500, res=300)
boxplot(B[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes,
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "B compartment", col="blue")
lines(minutes, B_med, col="black")
dev.off()

# Overplot
png("plot_LOS/A_B_compartment_LOS.png", height=1800, width=2500, res=300)
boxplot(A[4:9], ylim=c(-0.3,0.9), at=minutes, boxwex=50, names=minutes,
        xlab="Minutes of DpnII Digestion", outline=FALSE,
        ylab= "Loss of Structure", main= "A and B compartments", col="red")
lines(minutes, A_med, col="red")

boxplot(B[4:9], at=minutes, boxwex=50,
        col="blue", add=TRUE, outline=FALSE, names=FALSE)
lines(minutes, B_med, col="blue")
dev.off()
# 





