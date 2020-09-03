library(ggplot2)

sub <- read.table("SNIPER_K562_subcompartments_40kb.bed", sep="\t", header=TRUE)

comp_size <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/compartment_size/filter1000/eigen1_40kb.size.bedGraph",
                        sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "eigen", "size"),
                        na.strings = "nan")

m <- cbind(comp_size, sub[4])
m <- na.omit(m)

A <- m[m$eigen > 0,]
B <- m[m$eigen < 0,]

# A compartment
A_sizes <- sort(unique(A$size))

A_df <- data.frame()
for (s in A_sizes) {
  t <- table(A$sub[A$size == s])
  A_df <- rbind(A_df, t)
}
colnames(A_df) <- rownames(t)
rownames(A_df) <- A_sizes

Asizes <- factor(rep(A_sizes, each=5),
                 levels=A_sizes,
                 ordered=TRUE)
Asubs <- factor(rep(rownames(t), 166),
                levels=rownames(t),
                ordered=TRUE)

counts <- c()
for (i in 1:length(Asizes)) {
  counts <- c(counts, A_df[as.character(Asizes[i]), Asubs[i]])
}

stack_df <- data.frame(Asizes, Asubs, counts)
png("subcompartment_size_stacked_A.png", width=2500, 1200, res=300)

ggplot(stack_df, aes(fill=Asubs, y=counts, x=Asizes)) + 
  geom_bar(stat="identity", position="fill") +
  xlab("A compartment size") + ylab("Percent of Subcompartment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=3)) +
        scale_fill_manual(values=c("A1" = "#ED2024",
                             "A2" = "#F06BA8",
                             "B1" = "#8150A0",
                             "B2" = "#6FCCDD",
                             "B3" = "#3953A4"),
                    breaks=c("A1", "A2", "B1", "B2", "B3"))
        
          
dev.off()

# B compartment
B_sizes <- sort(unique(B$size))

B_df <- data.frame()
for (s in B_sizes) {
  t <- table(B$sub[B$size == s])
  B_df <- rbind(B_df, t)
}
colnames(B_df) <- rownames(t)
rownames(B_df) <- B_sizes

Bsizes <- factor(rep(B_sizes, each=5),
                 levels=B_sizes,
                 ordered=TRUE)
Bsubs <- factor(rep(rownames(t), 150),
                levels=rownames(t),
                ordered=TRUE)

counts <- c()
for (i in 1:length(Bsizes)) {
  counts <- c(counts, B_df[as.character(Bsizes[i]), Bsubs[i]])
}

stack_df <- data.frame(Bsizes, Bsubs, counts)
png("subcompartment_size_stacked_B.png", width=2500, 1200, res=300)

ggplot(stack_df, aes(fill=Bsubs, y=counts, x=Bsizes)) + 
  geom_bar(stat="identity", position="fill") +
  xlab("B compartment size") + ylab("Percent of Subcompartment") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5, size=3)) +
  scale_fill_manual(values=c("A1" = "#ED2024",
                             "A2" = "#F06BA8",
                             "B1" = "#8150A0",
                             "B2" = "#6FCCDD",
                             "B3" = "#3953A4"),
                    breaks=c("A1", "A2", "B1", "B2", "B3"))


dev.off()

