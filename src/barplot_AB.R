wkdir <- getwd()
DN_df<- read.table(paste(wkdir, "/eigenI_eigenJ_DN_interactome.txt", sep=""), sep="\t", header=TRUE)
MN_df<- read.table(paste(wkdir, "/eigenI_eigenJ_MN_interactome.txt", sep=""), sep="\t", header=TRUE)


DN_AA <- DN_df$cScore[DN_df$eigen_I > 0 & DN_df$eigen_J > 0]
DN_BB <- DN_df$cScore[DN_df$eigen_I < 0 & DN_df$eigen_J < 0]
DN_AB <- DN_df$cScore[(DN_df$eigen_I < 0 & DN_df$eigen_J > 0) | (DN_df$eigen_I > 0 & DN_df$eigen_J < 0)]

MN_AA <- MN_df$cScore[MN_df$eigen_I > 0 & MN_df$eigen_J > 0]
MN_BB <- MN_df$cScore[MN_df$eigen_I < 0 & MN_df$eigen_J < 0]
MN_AB <- MN_df$cScore[(MN_df$eigen_I < 0 & MN_df$eigen_J > 0) | (MN_df$eigen_I > 0 & MN_df$eigen_J < 0)]

png("barplot_AB.png", width=2500, height=2500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))

compartment_table <- matrix(c(log(mean(MN_AA, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE)),
                              log(mean(MN_BB, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE)),
                              log(mean(DN_AA, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE)),
                              log(mean(DN_BB, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE))),
                              ncol=2, byrow=TRUE)
colnames(compartment_table) <- c("AA/AB", "BB/AB")
rownames(compartment_table) <- c("Non-fragmented", "Fragmented")

compartment_table <- as.table(compartment_table)

barplot(compartment_table, beside=TRUE, col = c("black", "gray"),
        ylab = "Log Ratio of enrichment over AB interactions", cex.names=2,ylim=c(0,2), cex.axis=1.5, cex.lab=1.5)
legend("topright", rownames(compartment_table), fill= c("black", "gray"), cex=1.5)
dev.off()

        
print("Fold decrease log:")
print("MN_AA --> DN_AA")
print((log(mean(MN_AA, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE)))/(log(mean(DN_AA, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE))))
print("MN_BB --> DN_BB")
print((log(mean(MN_BB, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE)))/(log(mean(DN_BB, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE))))

print("Fold decrease:")
print("MN_AA --> DN_AA")
print((mean(MN_AA, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE))/((mean(DN_AA, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE))))
print("MN_BB --> DN_BB")
print((mean(MN_BB, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE))/(mean(DN_BB, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE)))

png("barplot_AB_nonlog.png", width=2500, height=2500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))

compartment_table <- matrix(c(mean(MN_AA, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE),
                              mean(MN_BB, na.rm = TRUE) / mean(MN_AB, na.rm = TRUE),
                              mean(DN_AA, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE),
                              mean(DN_BB, na.rm = TRUE) / mean(DN_AB, na.rm = TRUE)),
                            ncol=2, byrow=TRUE)
colnames(compartment_table) <- c("AA/AB", "BB/AB")
rownames(compartment_table) <- c("Non-fragmented", "Fragmented")

compartment_table <- as.table(compartment_table)

barplot(compartment_table, beside=TRUE, col = c("black", "gray"),
        ylab = "Ratio of enrichment over AB interactions", cex.names=2,ylim=c(0,5.5), cex.axis=1.5, cex.lab=1.5)
legend("topright", rownames(compartment_table), fill= c("black", "gray"), cex=1.5)
dev.off()
