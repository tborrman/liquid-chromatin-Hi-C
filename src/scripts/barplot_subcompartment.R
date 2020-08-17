
#MN
MN_df<- read.table("HBCRACKHiC-K562-MN-R1__hg19__genome__C-100000-iced__allchr_merge_pairwise_subcompartment.txt", sep="\t", header=FALSE)
colnames(MN_df) <-  c("chrom", "I_start", "J_start", "I_sub", "J_sub", "IF")
MN_A1 <- MN_df$IF[MN_df$I_sub == "A1" & MN_df$J_sub == "A1"]
print("step 1")
MN_A2 <- MN_df$IF[MN_df$I_sub == "A2" & MN_df$J_sub == "A2"]
MN_B1 <- MN_df$IF[MN_df$I_sub == "B1" & MN_df$J_sub == "B1"]
MN_B2 <- MN_df$IF[MN_df$I_sub == "B2" & MN_df$J_sub == "B2"]
MN_B3 <- MN_df$IF[MN_df$I_sub == "B3" & MN_df$J_sub == "B3"]
MN_A1B3 <- MN_df$IF[(MN_df$I_sub == "A1" & MN_df$J_sub == "B3") | (MN_df$I_sub == "B3" & MN_df$J_sub == "A1")]


#DN
DN_df<- read.table("HBCRACKHiC-K562-DN__hg19__genome__C-100000-iced__allchr_merge_pairwise_subcompartment.txt", sep="\t", header=FALSE)
colnames(DN_df) <- c("chrom", "I_start", "J_start", "I_sub", "J_sub", "IF")
DN_A1 <- DN_df$IF[DN_df$I_sub == "A1" & DN_df$J_sub == "A1"]
print("step 2")
DN_A2 <- DN_df$IF[DN_df$I_sub == "A2" & DN_df$J_sub == "A2"]
DN_B1 <- DN_df$IF[DN_df$I_sub == "B1" & DN_df$J_sub == "B1"]
DN_B2 <- DN_df$IF[DN_df$I_sub == "B2" & DN_df$J_sub == "B2"]
DN_B3 <- DN_df$IF[DN_df$I_sub == "B3" & DN_df$J_sub == "B3"]
DN_A1B3 <- DN_df$IF[(DN_df$I_sub == "A1" & DN_df$J_sub == "B3") | (DN_df$I_sub == "B3" & DN_df$J_sub == "A1")]

png("barplot_subcompartment.png", width=4000, height=2500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))

compartment_table <- matrix(c(log(mean(MN_A1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)),
                              log(mean(MN_A2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)),
                              log(mean(MN_B1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)),
                              log(mean(MN_B2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)),
                              log(mean(MN_B3, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)),
                              log(mean(DN_A1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)),
                              log(mean(DN_A2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)),
                              log(mean(DN_B1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)),
                              log(mean(DN_B2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)),
                              log(mean(DN_B3, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))),
                              ncol=5, byrow=TRUE)



colnames(compartment_table) <- c("A1/A1B3", "A2/A1B3","B1/A1B3", "B2/A1B3", "B3/A1B3")
rownames(compartment_table) <- c("MN", "DN")

compartment_table <- as.table(compartment_table)

barplot(compartment_table, beside=TRUE, col = c("red", "lightpink1", "red", "lightpink1", "blue", "lightblue", "blue", "lightblue", "blue", "lightblue" ),
        cex.names=2, cex.axis = 1.5, cex.lab=1.5,
        ylab = "Log Ratio of enrichment over AB interactions", ylim=c(0,5.5))
legend("topright", c("Non-fragmented A", "Fragmented A", "Non-fragmented B", "Fragmented B"), fill= c("red", "lightpink1", "blue", "lightblue"), cex=1.5)
dev.off()

print("Fold decrease log: ")
print((log(mean(MN_A1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE))) / (log(mean(DN_A1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)))) 
print((log(mean(MN_A2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE))) / (log(mean(DN_A2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))))
print((log(mean(MN_B1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE))) / (log(mean(DN_B1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))))
print((log(mean(MN_B2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE))) / (log(mean(DN_B2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))))
print((log(mean(MN_B3, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE))) / (log(mean(DN_B3, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))))

print("Fold decrease: ")
print((mean(MN_A1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)) / (mean(DN_A1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE))) 
print((mean(MN_A2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)) / (mean(DN_A2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)))
print((mean(MN_B1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)) / (mean(DN_B1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)))
print((mean(MN_B2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)) / (mean(DN_B2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)))
print((mean(MN_B3, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE)) / (mean(DN_B3, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)))

png("barplot_subcompartment_nonlog.png", width=4000, height=2500, res=300)
par(lwd=2, mar=c(5.1, 5.1, 4.1, 2.1))

compartment_table <- matrix(c(mean(MN_A1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE),
                              mean(MN_A2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE),
                              mean(MN_B1, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE),
                              mean(MN_B2, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE),
                              mean(MN_B3, na.rm = TRUE) / mean(MN_A1B3, na.rm = TRUE),
                              mean(DN_A1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE),
                              mean(DN_A2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE),
                              mean(DN_B1, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE),
                              mean(DN_B2, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE),
                              mean(DN_B3, na.rm = TRUE) / mean(DN_A1B3, na.rm = TRUE)),
                            ncol=5, byrow=TRUE)



colnames(compartment_table) <- c("A1/A1B3", "A2/A1B3","B1/A1B3", "B2/A1B3", "B3/A1B3")
rownames(compartment_table) <- c("MN", "DN")

compartment_table <- as.table(compartment_table)

barplot(compartment_table, beside=TRUE, col = c("red", "lightpink1", "red", "lightpink1", "blue", "lightblue", "blue", "lightblue", "blue", "lightblue" ),
        cex.names=2, cex.axis = 1.5, cex.lab=1.5,
        ylab = "Ratio of enrichment over AB interactions")
legend("topright", c("Non-fragmented A", "Fragmented A", "Non-fragmented B", "Fragmented B"), fill= c("red", "lightpink1", "blue", "lightblue"), cex=1.5)
dev.off()
