
TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBTD-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.59_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_TD"))
NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/LOS6Mb/DL-20200103-LCHCSeq-K562-PIIBNT-DpnII-R1-T1-filter1000__hg19__genome__C-40000-iced_scaleBy_0.65_range6Mb_LOS_removed_outliers.bedGraph",
                        sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "LOS_NT"))

sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                  sep="\t", header=TRUE)
df <- merge(TD, NT, by=c("chrom", "start", "end"))
df <- merge(df, sub, by=c("chrom", "start", "end"))

df <- na.omit(df)

# Stratify
A1 <- df[df$sub=="A1",]
A2 <- df[df$sub=="A2",]
B1 <- df[df$sub=="B1",]
B2 <- df[df$sub=="B2",]
B3 <- df[df$sub=="B3",]



pdf('subcompartment_LOS_PIIB_boxplot.pdf', width=5, height=6)
par(mar=c(6,4,1,1) + 0.1)
boxplot(A1$LOS_NT, A1$LOS_TD, A2$LOS_NT, A2$LOS_TD,
        B1$LOS_NT, B1$LOS_TD, B2$LOS_NT, B2$LOS_TD,
        B3$LOS_NT, B3$LOS_TD,
        col=c("red", "red", "hotpink", "hotpink",
              "purple", "purple", "cyan", "cyan",
              "blue", "blue"), outline=FALSE,
        names=c("Not Treated", "Pol II Block", "Not Treated", "Pol II Block",
                "Not Treated", "Pol II Block","Not Treated", "Pol II Block",
                "Not Treated", "Pol II Block"),
        ylab = "LOS DpnII (%)", las=2)
# text(1,median(A1_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A1_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(2,median(A1_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A1_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(3,median(A2_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A2_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(4,median(A2_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(A2_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(5,median(B1_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B1_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(6,median(B1_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B1_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(7,median(B2_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B2_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(8,median(B2_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B2_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(9,median(B3_highFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B3_highFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
# text(10,median(B3_lowFPKM$hl, na.rm=TRUE) + 2, paste("n=", length(B3_lowFPKM$hl), sep=""), 
#      col= "white", cex = 0.5)
legend("bottomleft", c("A1", "A2", "B1", "B2", "B3"), 
       fill=c("red", "hotpink", "purple", "cyan", "blue"), bty="n")
dev.off()

t.test(A1$LOS_NT, A1$LOS_TD)
t.test(A2$LOS_NT, A2$LOS_TD)
t.test(B1$LOS_NT, B1$LOS_TD)
t.test(B2$LOS_NT, B2$LOS_TD)
t.test(B3$LOS_NT, B3$LOS_TD)



