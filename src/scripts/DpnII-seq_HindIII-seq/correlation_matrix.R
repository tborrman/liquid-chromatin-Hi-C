library(pheatmap)
library(RColorBrewer)

# 40kb coverage correlation matrix
path <- "../coverage/40kb/"
files <- c("HBCRACKHiC-K562-DN-TD-R1_GCCAAT_L008_coverage_40kb.bed",
           "HBCRACKHiC-K562-DN-TD-R2_ATGTCA_L008_coverage_40kb.bed",
           "HBCRACKHiC-K562-HN-TD-R1_ACAGTG_L008_coverage_40kb.bed",
           "HBDpSeqK562-DN15mR1_S2_L001_coverage_40kb.bed",
           "HBDpSeqK562-DN16hR1_S11_L003_coverage_40kb.bed",
           "HBDpSeqK562-DN1hR1_S5_L002_coverage_40kb.bed",
           "HBDpSeqK562-DN2hR1_S6_L002_coverage_40kb.bed",
           "HBDpSeqK562-DN30mR1_S3_L001_coverage_40kb.bed",
           "HBDpSeqK562-DN3hR1_S7_L002_coverage_40kb.bed",
           "HBDpSeqK562-DN45mR1_S4_L001_coverage_40kb.bed",
           "HBDpSeqK562-DN4hR1_S8_L002_coverage_40kb.bed",
           "HBDpSeqK562-DN5mR1_S1_L001_coverage_40kb.bed",
           "HBDpSeqK562-DN75mR1_S9_L003_coverage_40kb.bed",
           "HBDpSeqK562-DN90mR1_S10_L003_coverage_40kb.bed")
labels <- c("DN_240min_R1", "DN_240min_R2", "HN_240min_R1", "DN_15min_R3", "DN_960min_R3", "DN_60min_R3", "DN_120min_R3", "DN_30min_R3",
            "DN_180min_R3", "DN_45min_R3", "DN_240min_R3", "DN_5min_R3", "DN_75min_R3", "DN_90min_R3")

coverage_df <- read.table(paste(path, files[1], sep=""), sep="\t", header=FALSE)
coverage_df <- coverage_df[1:4]
colnames(coverage_df) <- c("chrom", "start", "end", labels[1])
for (i in 2:length(files)) {
  x <- read.table(paste(path, files[i], sep=""), sep="\t", header=FALSE)
  x <- x[1:4]
  colnames(x) <- c("chrom", "start", "end", labels[i])
  coverage_df <- cbind(coverage_df, x[labels[i]])
} 

cor_matrix <- as.data.frame(cor(coverage_df[4:length(coverage_df)], use="pairwise.complete.obs", method="spearman"))
pdf("correlation_matrix/coverage_correlation_matrix.pdf", width=8, height=7.5, onefile=FALSE)
pheatmap(cor_matrix, color=rev(colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(100)), 
         breaks=seq(-0.2,1,1.2/100))
dev.off()

# 40kb copy correct coverage correlation matrix
path <- "../copy_correct_coverage/40kb/"
copy_files <- c()
for (f in files) {
  cf <- gsub("coverage", "copy_correct_coverage", f)
  copy_files <- c(copy_files, cf)
}

cc_coverage_df <- read.table(paste(path, copy_files[1], sep=""), sep="\t", header=FALSE)
colnames(cc_coverage_df) <- c("chrom", "start", "end", labels[1])
for (i in 2:length(copy_files)) {
  x <- read.table(paste(path, copy_files[i], sep=""), sep="\t", header=FALSE)
  colnames(x) <- c("chrom", "start", "end", labels[i])
  cc_coverage_df <- cbind(cc_coverage_df, x[labels[i]])
}
# Order to mactch clustering of pre copy number correction
ordered_cc_coverage_df <- cc_coverage_df[c("HN_240min_R1",
                                           "DN_180min_R3",
                                           "DN_960min_R3",
                                           "DN_240min_R3",
                                           "DN_90min_R3",
                                           "DN_120min_R3",
                                           "DN_75min_R3",
                                           "DN_45min_R3",
                                           "DN_60min_R3",
                                           "DN_30min_R3",
                                           "DN_15min_R3",
                                           "DN_5min_R3",
                                           "DN_240min_R1",
                                           "DN_240min_R2")]

cc_cor_matrix <- as.data.frame(cor(ordered_cc_coverage_df, use="pairwise.complete.obs", method="spearman"))
pdf("correlation_matrix/copy_correct_coverage_correlation_matrix.pdf", width=8, height=7.5, onefile=FALSE)
pheatmap(cc_cor_matrix, color=rev(colorRampPalette((brewer.pal(n = 7, name ="RdYlBu")))(100)),
         breaks=seq(-0.2,1,1.2/100), cluster_rows = FALSE, cluster_cols = FALSE)
dev.off()


