digest_table <- read.table("C:/cygwin64/home/Tyler/Research/digest/dpnII/hg19/GMM/GMM_post_probs.txt", sep="\t", header=TRUE)
# NOTE! The columns will change if I rerun GMM_copy_number.py!"
colnames(digest_table) <- c("chrom", "start", "end", "reads", "trip", "dip", "tetra")

for (i in 1:nrow(digest_table)) {
  # trip
  if (which.max(digest_table[i,5:7]) == 1) {
    digest_table[i, "reads"] <- digest_table[i, "reads"]/ 1.5
    
  } 
  # tetra
  else if (which.max(digest_table[i,5:7]) == 3) {
    digest_table[i, "reads"] <- digest_table[i, "reads"]/ 2.0
  }
  
}

write.table(digest_table, file="DpnII-seq_copy_correct_500kb.bed", sep="\t", quote=FALSE,
            row.names=FALSE)

