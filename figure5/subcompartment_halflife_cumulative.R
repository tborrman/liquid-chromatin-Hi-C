library(ggplot2)
hl <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                 sep="\t", header=FALSE)
colnames(hl) <- c("chrom", "start", "end", "hl")
sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                  sep="\t", header=TRUE)

dfna <- merge(hl, sub, by=c("chrom", "start", "end"))
df <- na.omit(dfna)

# Subcompartments
A1 <- df[df$sub =="A1",]
A2 <- df[df$sub =="A2",]
B1 <- df[df$sub =="B1",]
B2 <- df[df$sub =="B2",]
B3 <- df[df$sub =="B3",]

pdf("subcompartment_halflife_cumulative.pdf", width=6, height=3)

ggplot() + 
  stat_ecdf(data=A1, aes(x=hl, color="A1"), geom="line") +
  stat_ecdf(data=A2, aes(x=hl, color="A2"), geom="line") +
  stat_ecdf(data=B1, aes(x=hl, color="B1"), geom="line") +
  stat_ecdf(data=B2, aes(x=hl, color="B2"), geom="line") +
  stat_ecdf(data=B3, aes(x=hl, color="B3"), geom="line") +
  xlim(40, 135) +
  xlab(bquote("t"[1/2] ~ "(minutes)")) + ylab("Fraction of Data") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        legend.key = element_rect(fill="white"), legend.title = element_blank()) +
  scale_color_manual(values=c("A1" = "red",
                              "A2" = "hotpink",
                              "B1" = "purple",
                              "B2" = "cyan",
                              "B3" = "blue"),
                     breaks=c("A1", "A2", "B1", "B2", "B3"))
dev.off()

