eigen <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/eigen/eigen1_40kb.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "PC1"))
eigen$PC1[is.na(eigen$PC1)] <- NA

path <- "C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/"
files <- c("HBDpSeqK562-DN5mR1_S1_L001_copy_correct_coverage_40kb.bed",
           "HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
           "HBDpSeqK562-DN2hR1_S6_L002_copy_correct_coverage_40kb.bed",
           "HBDpSeqK562-DN3hR1_S7_L002_copy_correct_coverage_40kb.bed",
           "HBDpSeqK562-DN4hR1_S8_L002_copy_correct_coverage_40kb.bed",
           "HBDpSeqK562-DN16hR1_S11_L003_copy_correct_coverage_40kb.bed")

labels <- c("DN_5min_R3", "DN_60min_R3", "DN_120min_R3", "DN_180min_R3",
            "DN_240min_R3", "DN_960min_R3")

df <- read.table(paste(path, files[1], sep=""), sep="\t", header=FALSE)
colnames(df) <- c("chrom", "start", "end", labels[1])
for (i in 2:length(files)) {
  x <- read.table(paste(path, files[i], sep=""), sep="\t", header=FALSE)
  colnames(x) <- c("chrom", "start", "end", labels[i])
  df <- cbind(df, x[labels[i]])
} 


c_df <- cbind(eigen["PC1"], df[4:ncol(df)])
cor_matrix <- as.data.frame(cor(c_df, use="pairwise.complete.obs", method="spearman"))

pdf("fig4_PC1_corr_DpnIIseq_barplot.pdf", width=6, height=7)
par(mar=c(8, 5, 4, 2) + 0.1)
barplot(cor_matrix$PC1[2:length(cor_matrix$PC1)], 
        names=colnames(cor_matrix)[2:length(colnames(cor_matrix))],
        col="grey60", las=2, ylim=c(0, 1), 
        ylab="Genome wide Spearman correlation with PC1",
        cex.lab=1.5)
dev.off()


