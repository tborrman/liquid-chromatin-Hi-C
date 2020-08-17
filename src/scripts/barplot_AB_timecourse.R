wkdir <- getwd()
DN1_df<- read.table(paste(wkdir, "/HBHiCK562DN10-1h_interactome.txt", sep=""), sep="\t", header=TRUE)
DN2_df<- read.table(paste(wkdir, "/HBHiCK562DN10-2h_interactome.txt", sep=""), sep="\t", header=TRUE)
DN3_df<- read.table(paste(wkdir, "/HBHiCK562DN10-3h_interactome.txt", sep=""), sep="\t", header=TRUE)
DN4_df<- read.table(paste(wkdir, "/HBHiCK562DN10-4h_interactome.txt", sep=""), sep="\t", header=TRUE)
DN5_df<- read.table(paste(wkdir, "/HBHiCK562DN10-5m_interactome.txt", sep=""), sep="\t", header=TRUE)
DNON_df<- read.table(paste(wkdir, "/HBHiCK562DN10-ON_interactome.txt", sep=""), sep="\t", header=TRUE)
MN_df<- read.table(paste(wkdir, "/HBHiC-K562-MN-Dp-1_interactome.txt", sep=""), sep="\t", header=TRUE)


get_means <- function(df) {
  AA <- mean(df$cScore[df$eigen_I > 0 & df$eigen_J > 0], na.rm = TRUE)
  BB <- mean(df$cScore[df$eigen_I < 0 & df$eigen_J < 0], na.rm = TRUE)
  AB <- mean(df$cScore[(df$eigen_I < 0 & df$eigen_J > 0) | (df$eigen_I > 0 & df$eigen_J < 0)], na.rm=TRUE)
  return(c(AA, BB, AB))
}


DN1_means <- get_means(DN1_df)
DN2_means <- get_means(DN2_df)
DN3_means <- get_means(DN3_df)
DN4_means <- get_means(DN4_df)
DN5_means <- get_means(DN5_df)
DNON_means <- get_means(DNON_df)
MN_means <- get_means(MN_df)

png("barplot_timecourse_AB.png", width = 3500, height=2000, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))
compartment_table <- matrix(c(MN_means[1], DN5_means[1], DN1_means[1], DN2_means[1], DN3_means[1], DN4_means[1], DNON_means[1],
                              MN_means[2], DN5_means[2], DN1_means[2], DN2_means[2], DN3_means[2], DN4_means[2], DNON_means[2],
                              MN_means[3], DN5_means[3], DN1_means[3], DN2_means[3], DN3_means[3], DN4_means[3], DNON_means[3]),
                              ncol=7, byrow=TRUE)
colnames(compartment_table) <- c("MN", "DN-5m", "DN-1h", "DN-2h", "DN-3h", "DN-4h", "DN-ON")
rownames(compartment_table) <- c("AA", "BB", "AB")

compartment_table <- as.table(compartment_table)
barplot(compartment_table, beside=TRUE, col=c("red", "blue", "mediumorchid4"), 
        ylab = "Mean interactions", cex.names=1.5, cex.axis = 1.5, cex.lab = 1.5, ylim=c(0,800))
legend("topright", rownames(compartment_table), fill= c("red","blue", "mediumorchid4"), cex=1.5)
dev.off()


