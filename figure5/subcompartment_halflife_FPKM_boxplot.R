hl <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                 sep="\t", header=FALSE)
colnames(hl) <- c("chrom", "start", "end", "hl")
FPKM <- read.table("FPKM_40kb.bedGraph",
                   sep="\t", header=FALSE)
colnames(FPKM) <- c("chrom", "start", "end", "FPKM")
sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                  sep="\t", header=TRUE)
df <- merge(hl, FPKM, by=c("chrom", "start", "end"))
df <- merge(df, sub, by=c("chrom", "start", "end"))

df <- na.omit(df)

# Stratify
A1_highFPKM <- df[df$sub=="A1" & df$FPKM > 1,]
A1_lowFPKM <- df[df$sub=="A1" & df$FPKM <= 1,]
A2_highFPKM <- df[df$sub=="A2" & df$FPKM > 1,]
A2_lowFPKM <- df[df$sub=="A2" & df$FPKM <= 1,]
B1_highFPKM <- df[df$sub=="B1" & df$FPKM > 1,]
B1_lowFPKM <- df[df$sub=="B1" & df$FPKM <= 1,]
B2_highFPKM <- df[df$sub=="B2" & df$FPKM > 1,]
B2_lowFPKM <- df[df$sub=="B2" & df$FPKM <= 1,]
B3_highFPKM <- df[df$sub=="B3" & df$FPKM > 1,]
B3_lowFPKM <- df[df$sub=="B3" & df$FPKM <= 1,]


pdf('subcompartment_halflife_FPKM_boxplot.pdf', width=5, height=6)
boxplot(A1_highFPKM$hl, A1_lowFPKM$hl, A2_highFPKM$hl, A2_lowFPKM$hl,
        B1_highFPKM$hl, B1_lowFPKM$hl, B2_highFPKM$hl, B2_lowFPKM$hl,
        B3_highFPKM$hl, B3_lowFPKM$hl,
        col=c("red", "red", "hotpink", "hotpink",
              "purple", "purple", "cyan", "cyan",
              "blue", "blue"), outline=FALSE,
        names=c("FPKM > 1", "FPKM <= 1", "FPKM > 1", "FPKM <= 1",
                "FPKM > 1", "FPKM <= 1","FPKM > 1", "FPKM <= 1",
                "FPKM > 1", "FPKM <= 1"),
        ylab = bquote("t"[1/2] ~ "(minutes)"), las=2, ylim=c(30,130))
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
legend("topleft", c("A1", "A2", "B1", "B2", "B3"), 
       fill=c("red", "hotpink", "purple", "cyan", "blue"), bty="n")
dev.off()

t.test(A1_highFPKM$hl, A1_lowFPKM$hl)
t.test(A2_highFPKM$hl, A2_lowFPKM$hl)
t.test(B1_highFPKM$hl, B1_lowFPKM$hl)
t.test(B2_highFPKM$hl, B2_lowFPKM$hl)
t.test(B3_highFPKM$hl, B3_lowFPKM$hl)

