#!/usr/bin/env Rscript

library(argparse)


parser <- ArgumentParser(description = "build table to associate histone mods with interactome")
parser$add_argument("-i", help="3-column tab delimited interaction file (cols = [loci_i, loci_j, interaction_score])", type="character", required=TRUE)
parser$add_argument("-s", help= "histone signal bed file, binned at same resolution as interaction file", type="character", required=TRUE)
parser$add_argument("-o", help="name for table output file", type="character", required=TRUE)


args <- parser$parse_args()

#Interactome
interact_df <- read.table(args$i, sep="\t", header=TRUE)


# Clean up
# change to only start values
colnames(interact_df) <- c("Istart", "Jstart", "cScore")
matches <- regmatches(interact_df$Istart, regexpr(":[0-9]+",interact_df$Istart))
Istart <- substr(matches, 2, nchar(matches))

matches <- regmatches(interact_df$Jstart, regexpr(":[0-9]+",interact_df$Jstart))
Jstart <- substr(matches, 2, nchar(matches))

matches <- regmatches(interact_df$Istart, regexpr("chr[0-9X]+", interact_df$Istart))
chrom  <- matches

interact_df$Istart <- Istart
interact_df$Jstart <- Jstart
interact_df <- cbind(chrom, interact_df)

#H3K9me3
H3K9me3_df <- read.table(args$s, sep="\t", header=TRUE)
H3K9me3_df$start = H3K9me3_df$start + 1

#change column names
write(c("chrom","I",	"J",	"cScore",	"signal_I",	"signal_J"), sep="\t",ncolumns=6, append=TRUE, file=args$o)
for (i in 1:nrow(interact_df)) {
  row <- interact_df[i,]
  loc1 <- as.numeric(row$Istart)
  loc2 <- as.numeric(row$Jstart)
  chrom <- as.character(row$chrom)
  
  # Look up histone values for loci
  loc1_his <- H3K9me3_df[H3K9me3_df$start == loc1 & H3K9me3_df$chr == chrom, "score"]
  loc2_his <- H3K9me3_df[H3K9me3_df$start == loc2 & H3K9me3_df$chr == chrom, "score"]
  
  
  new_row = c(chrom, loc1, loc2, as.numeric(as.character(row$cScore)), loc1_his, loc2_his)
  
  write(new_row, sep="\t",ncolumns=6, append=TRUE, file=args$o)
  
  if (i%%1000 == 0) {
    
    
    print(paste("On row ", i, sep=""))
    
  }
} 


