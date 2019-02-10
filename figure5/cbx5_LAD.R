df <- read.table("C:/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt",
                 sep="\t", header=TRUE, check.names=FALSE)

dfNA <- df[!(is.na(df$`half-life_LOS`) | is.na(df$LAD_K562) | is.na(df$CBX5) | is.na(df$H3K9me3_R1)),]

# Get rid of cbx5, LADs, and t1/2
png("cbx5_hist.png", height=2000, width=2000, res=300)
hist(dfNA$CBX5, breaks=200, xlim=c(0,3))
dev.off()

png("LADs_hist.png", height=2000, width=2000, res=300)
hist(dfNA$LAD_K562, breaks=200)
dev.off()

png("H3K9me3_hist.png", height=2000, width=2000, res=300)
hist(dfNA$H3K9me3_R1, breaks=200, xlim=c(0,3))
dev.off()

mycols <- c("chrom", "start","end", "half-life_LOS", "CBX5", "LAD_K562", "H3K9me3_R1")

CBX5 <- dfNA[dfNA$CBX5 > 1.5 & dfNA$LAD_K562 < 0, mycols]

LADs <- dfNA[dfNA$LAD_K562 > 0.75 & dfNA$CBX5 < 1, mycols]

H3K9me3 <- dfNA[dfNA$H3K9me3_R1 > 2 & dfNA$LAD_K562 < 0, mycols]

LADsnoH <- dfNA[dfNA$LAD_K562 > 0.75 & dfNA$H3K9me3_R1 < 1.25, mycols]

png('cbx5_LAD_boxplot.png', width=1000, height=2500, res=300)
par(mar=c(9,4,4,2) + 0.1)
boxplot(CBX5$`half-life_LOS`, LADs$`half-life_LOS`, H3K9me3$`half-life_LOS`, LADsnoH$`half-life_LOS`,
        col=c("grey", "grey", "grey", "grey"), outline=FALSE,
        names=c("CBX5_no_LADs", "LADs_no_CBX5", "H3K9me3_no_LADs", "LADs_no_H3K9me3"),
        ylab = bquote("t"[1/2] ~ "(minutes)"), las=2, ylim=c(30,140))
text(1,median(CBX5$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(CBX5$`half-life_LOS`),sep=""), 
     col= "black", cex = 0.75)
text(2,median(LADs$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(LADs$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
text(3,median(H3K9me3$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(H3K9me3$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
text(4,median(LADsnoH$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(LADsnoH$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
dev.off()

############################################################################################
# JUST B compartment
#############################################################################################
dfB <- df[df$PCA_eigen1 < 0,]
dfNA <- dfB[!(is.na(dfB$`half-life_LOS`) | is.na(dfB$LAD_K562) | is.na(dfB$CBX5) | is.na(dfB$H3K9me3_R1)),]


mycols <- c("chrom", "start","end", "half-life_LOS", "CBX5", "LAD_K562", "H3K9me3_R1", "PCA_eigen1")

CBX5 <- dfNA[dfNA$CBX5 > 1.5 & dfNA$LAD_K562 < 0, mycols]

LADs <- dfNA[dfNA$LAD_K562 > 0.75 & dfNA$CBX5 < 1, mycols]

H3K9me3 <- dfNA[dfNA$H3K9me3_R1 > 2 & dfNA$LAD_K562 < 0, mycols]

LADsnoH <- dfNA[dfNA$LAD_K562 > 0.75 & dfNA$H3K9me3_R1 < 1.25, mycols]

png('cbx5_LAD_boxplot_Bcompartment.png', width=1000, height=2500, res=300)
par(mar=c(9,4,4,2) + 0.1)
boxplot(CBX5$`half-life_LOS`, LADs$`half-life_LOS`, H3K9me3$`half-life_LOS`, LADsnoH$`half-life_LOS`,
        col=c("dodgerblue", "dodgerblue", "dodgerblue", "dodgerblue"), outline=FALSE,
        names=c("CBX5_no_LADs", "LADs_no_CBX5", "H3K9me3_no_LADs", "LADs_no_H3K9me3"),
        ylab = bquote("t"[1/2] ~ "(minutes)"), las=2, ylim=c(30,140))
text(1,median(CBX5$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(CBX5$`half-life_LOS`),sep=""), 
     col= "black", cex = 0.75)
text(2,median(LADs$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(LADs$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
text(3,median(H3K9me3$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(H3K9me3$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
text(4,median(LADsnoH$`half-life_LOS`, na.rm=TRUE) + 2, paste("n=", length(LADsnoH$`half-life_LOS`), sep=""), 
     col= "black", cex = 0.75)
dev.off()


