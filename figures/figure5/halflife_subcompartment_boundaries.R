df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt",
  sep="\t", header=TRUE, check.names=FALSE)
sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed", 
  sep="\t", header=TRUE)

tomerge <- df[c("chrom", "start", "end", "half-life_LOS", "H3K27me3_R1")]

m <- cbind(tomerge, sub[4])


get_bin_size <- function(r) {
  bs <- as.numeric(r["end"]) - as.numeric(r["start"])
  return(bs)
}

get_boundary_df <- function(d, feature, windowsize, binsize, C1, C2) {
  # Get dataframe of feature signal for windowsize 
  # upstram and downstream of compartment switch from C1 to C2
  
  # Create dataframe of windowsize*2 window around C1->CB boundary switches
  C1toC2_feature <- data.frame()
  windowbins <- windowsize/binsize
  genomebins <- nrow(d)
  for (i in 1:(genomebins - 1)) {
    # Make sure windowsize from start and end of genome
    if (i > windowbins & i < (genomebins - windowbins)) {
      # Make sure window size is within chromosome
      ichrom <- as.character(d[i, "chrom"])
      if (as.character(d[(i-windowbins),"chrom"]) ==  ichrom &
          as.character(d[(i+windowbins),"chrom"]) == ichrom) { 
        # Check if C1 -> C2 switch
        isub <- as.character(d[i, "sub"])
        nextsub <- as.character(d[i+1, "sub"])
        if ((!is.na(isub)) & (!is.na(nextsub))) {
          if (isub == C1 & nextsub == C2) {
            #if(d[i,"size"] > 1000000 & d[i+1, "size"] > 1000000) {
            # Get data around 2Mb of eigen switch (boundary)
            boundary_feature <- d[(i-windowbins):(i+windowbins), feature]
            C1toC2_feature <- rbind(C1toC2_feature, boundary_feature)
            #}
          }
        }
      }
    } 
  }
  # Name columns
  cols <- ncol(C1toC2_feature)
  colnames(C1toC2_feature) <- c((floor(cols/2):0)*-bs/1000, (1:floor(cols/2))*bs/1000)
  return(C1toC2_feature)
}

boundary_plot <- function(t, h, switch) {
  pdf(paste("subcompartment_boundaries/compartment_boundaries_", switch, "_half-life_LOS_H3K27me3_merge_sizefix.pdf", sep=""), 
      height=6, width=6)
  par(mar = c(5, 4, 4, 4) + 0.1)
  mean_t <- apply(t, 2, mean, na.rm=TRUE)
  mean_h <- apply(h, 2, mean, na.rm=TRUE)
  bp <- as.numeric(colnames(t))
  plot(bp, mean_t, type="l", col="dodgerblue",
       xlab=paste(switch, " compartment boundary (kb)"), 
       ylab="", cex.axis=0.8)
  axis(2, col="dodgerblue", col.axis="dodgerblue", cex.axis=0.8, at=pretty(range(mean_t), 5))
  mtext(bquote("t"[1/2] ~ "(minutes)"), side=2, col="dodgerblue", line=2.5)
  par(new=TRUE)
  plot(bp, mean_h, type="l", col="black",
       axes=FALSE, bty="n", xlab="", ylab="")
  axis(side=4, at=pretty(range(mean_h), 5), cex.axis=0.8)
  mtext("H3K27me3", side=4, line=2.5)
  
  dev.off()
}

bs <- get_bin_size(df[1,])
ws <- 1000000   # window upstream and downstream of boundary switch

A1toB1_thalf <- get_boundary_df(m, "half-life_LOS", ws, bs, "A1", "B1")
A1toB1_H3K27me3 <- get_boundary_df(m, "H3K27me3_R1", ws, bs, "A1", "B1")
boundary_plot(A1toB1_thalf, A1toB1_H3K27me3, "A1-B1")

A2toB1_thalf <- get_boundary_df(m, "half-life_LOS", ws, bs, "A2", "B1")
A2toB1_H3K27me3 <- get_boundary_df(m, "H3K27me3_R1", ws, bs, "A2", "B1")
boundary_plot(A2toB1_thalf, A2toB1_H3K27me3, "A2-B1")

B2toB1_thalf <- get_boundary_df(m, "half-life_LOS", ws, bs, "B2", "B1")
B2toB1_H3K27me3 <- get_boundary_df(m, "H3K27me3_R1", ws, bs, "B2", "B1")
boundary_plot(B2toB1_thalf, B2toB1_H3K27me3, "B2-B1")

B3toB1_thalf <- get_boundary_df(m, "half-life_LOS", ws, bs, "B3", "B1")
B3toB1_H3K27me3 <- get_boundary_df(m, "H3K27me3_R1", ws, bs, "B3", "B1")
boundary_plot(B3toB1_thalf, B3toB1_H3K27me3, "B3-B1")



A1toA2_thalf <- get_boundary_df(m, "half-life_LOS", ws, bs, "A1", "A2")
A1toA2_H3K27me3 <- get_boundary_df(m, "H3K27me3_R1", ws, bs, "A1", "A2")
boundary_plot(A1toA2_thalf, A1toA2_H3K27me3, "A1-A2")


