library(RColorBrewer)
size_df <- read.table("eigen1_40kb.size.bedGraph",
                       sep="\t", header=FALSE)
colnames(size_df) <- c("chrom", "start", "end", "eigen", "size")
hl_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph", 
                    sep="\t", header=FALSE)
colnames(hl_df) <- c("chrom", "start", "end", "hl")


png("compartment_size_vs_half-life_40kb_filter1000.png", height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(size_df$size, hl_df$hl, pch=20, 
     col=ifelse(size_df$eigen < 0, rgb(red=0, green=0, blue=1, alpha=0.1), rgb(red=1, green=0, blue=0, alpha=0.1)),
     xlab="Compartment size (bp)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, main= "Genome" 
     #, xlim=c(0,20000000), ylim=c(50,350)
     )
dev.off()

A_s <- size_df[size_df$eigen > 0,]
A_hl <- hl_df[size_df$eigen > 0,]

png("compartment_size_vs_half-life_40kb_filter1000_A.png", height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(A_s$size/1000000, A_hl$hl, pch=20, col=rgb(red=1, green=0, blue=0, alpha=0.1),
     xlab="Compartment size (Mb)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5 , ylim=c(35, 150), xlim=c(0, 22)
     )

text(18, 100,
     bquote(rho ~ "="~  .(round(cor(A_s$size, A_hl$hl, method="spearman", use="complete.obs"), 2))),
     cex=1.5)
dev.off()


B_s <- size_df[size_df$eigen < 0,]
B_hl <- hl_df[size_df$eigen < 0,]

png(paste("compartment_size_vs_half-life_40kb_filter1000_B.png",sep=""), height=2500, width=2500, res=300)
par(mar=c(5,6,4,2) + 0.1)
plot(B_s$size/1000000, B_hl$hl, pch=20, col=rgb(red=0, green=0, blue=1, alpha=0.1),
     xlab="Compartment size (Mb)", ylab = bquote("t"[1/2] ~ "(minutes)"), cex.lab=1.5,
     cex.main=1.5, ylim=c(35, 150), xlim=c(0, 22)
     )
text(10, 140,
     bquote(rho ~ "="~  .(round(cor(B_s$size, B_hl$hl, method="spearman", use="complete.obs"), 2))),
     cex=1.5)

dev.off()


