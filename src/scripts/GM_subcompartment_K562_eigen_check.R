sub_df <- read.table("GSE63525_GM12878_subcompartments_sorted_100kb.bed", header=FALSE, sep="\t")
eigen_df <- read.table("HBCRACKHiC-K562-MN-R1__hg19__genome__C-500000_binned_to_100000-iced__all.zScore.eigen1.sorted.bedGraph", header=FALSE, sep="\t")

colnames(sub_df) <- c("chrom", "start" ,"end", "subcompartment")
colnames(eigen_df) <- c("chrom", "start" ,"end", "eigen")

merged_df <- merge(sub_df, eigen_df, c("chrom", "start"))

merged_df$compartment[merged_df$eigen < 0] <- "B"
merged_df$compartment[merged_df$eigen >= 0] <- "A"

table_df <- merged_df[c("compartment", "subcompartment")]
compartment_table <- table(table_df)
png("subcompartment_check_barplot.png", height= 2000, width= 3000, res=300)
barplot(compartment_table, beside=TRUE, col= c("red", "blue"), legend= rownames(compartment_table))
dev.off()
