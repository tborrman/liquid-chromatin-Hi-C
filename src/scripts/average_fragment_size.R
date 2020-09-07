# Generate bedfile of average fragment size after 4-hour DpnII-digest for each bin
options(digits=22)


# Sequencing data for DpnII fragments sliced at different size ranges from gel
slice1_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz1-K562-R1/bed/HB-Dp4h-siz1-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                        sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size1_signal"))
slice2_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz2-K562-R1/bed/HB-Dp4h-siz2-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                        sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size2_signal"))
slice3_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz3-K562-R1/bed/HB-Dp4h-siz3-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                        sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size3_signal"))

get_percent <- function(x) {
  # Transform sizex_signal to percent of total signal
  signal_sum = sum(x, na.rm=TRUE)
  p <- (x/signal_sum) * 100
  return(p)
}

get_ave_frag_size <- function(x) {
  # Calculate average fragment size for bin
  
  # Constants#################################################################
  # sx = mean size in bp of fragments from slice x
  s1 = 643
  s2 = 2332
  s3 = 5495
  # qx = quantity of slice x fragments in ng/uL
  q1 = 1.6562
  q2 = 2.544
  q3 = 2.4632
  ############################################################################
  
  if (any(is.na(x))) {
    avfs <- NA
  }
  else {
    # px = percent slice x reads mapped to bin
    p1 = x[1]
    p2 = x[2]
    p3 = x[3]
    
    if ((p1 == 0) & (p2 == 0) & (p3 == 0)) {
      avfs <- NA
    }
    else {
      avfs <- round(((p1*q1*s1) + (p2*q2*s2) + (p3*q3*s3)) / ((p1*q1) + (p2*q2) + (p3*q3)), 3) 
    }
  }
  return(avfs)
}

df <- data.frame(slice1_df, slice2_df["size2_signal"], slice3_df["size3_signal"])
signal_df <- df[4:6]

# Transform signals to percent of total signal
percent_df <- data.frame(apply(signal_df, 2, get_percent))
# Get average fragment size 
ave_frag_size <- apply(percent_df, 1, get_ave_frag_size)

options(digits=7)
out_df <- data.frame(df[1:3], ave_frag_size)

write.table(out_df, 
            "C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/ave_frag_size/average_fragment_size_40kb.bedGraph", 
            sep="\t", quote=FALSE, col.names=FALSE, row.names=FALSE)

