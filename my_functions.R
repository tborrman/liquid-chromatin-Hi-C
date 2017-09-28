# R functions for digest-Hi-C

OrderChromStart <- function(df) {
  # Order dataframe based on chromosome and start position
  #
  # Args:
  #  df: "chrom" column and "start" column
  #      human genome with no chrM
  #      chr1, chr2, ...., chrX, chrY
  #
  # Returns:
  # dataframe with rows sorted by chrom and start position
  df$chrom <- factor(df$chrom, levels=c(paste("chr", 1:22, sep=""), "chrX", "chrY"),
                     ordered = TRUE)
  ordered_df <- df[order(df$chrom, df$start),]
  return(ordered_df)
}
