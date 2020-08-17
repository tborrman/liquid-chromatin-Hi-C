df <- read.table("/cygwin64/home/Tyler/Research/digest/feature_analysis/C-40000/v4/feature_matrix_v4_40kb.txt", sep="\t",
                 header=TRUE, check.names=FALSE)
#data <- df[c("chrom", "start", "end", "half-life_LOS", "PCA_eigen1")]


get_bin_size <- function(r) {
  bs <- as.numeric(r["end"]) - as.numeric(r["start"])
  return(bs)
}

get_AtoB_df <- function(d, feature, windowsize, binsize) {
  # Get dataframe of feature signal for windowsize 
  # upstram and downstream of compartment switch from A to B
  
  # Create dataframe of 2Mb window around A->B boundary switches
  AtoB_feature <- data.frame()
  windowbins <- windowsize/binsize
  genomebins <- nrow(d)
  for (i in 1:(genomebins - 1)) {
    # Make sure 1Mb from start and end of genome
    if (i > windowbins & i < (genomebins - windowbins)) {
      # Make sure window size is within chromosome
      ichrom <- as.character(d[i, "chrom"])
      if (as.character(d[(i-windowbins),"chrom"]) ==  ichrom &
          as.character(d[(i+windowbins),"chrom"]) == ichrom) { 
        # Check if eigen switch
        ieigen <- as.numeric(d[i, "PCA_eigen1"])
        nexteigen <- as.numeric(d[i+1, "PCA_eigen1"])
        if ((!is.na(ieigen)) & (!is.na(nexteigen))) {
          if (ieigen > 0 & nexteigen < 0) {
            # Get data around 2Mb of eigen switch (boundary)
            boundary_feature <- d[(i-windowbins):(i+windowbins), feature]
            AtoB_feature <- rbind(AtoB_feature, boundary_feature)
          }
        }
      }
    } 
  }
  # Name columns
  cols <- ncol(AtoB_feature)
  colnames(AtoB_feature) <- c((floor(cols/2):0)*-bs/1000, (1:floor(cols/2))*bs/1000)
  return(AtoB_feature)
}

get_BtoA_df <- function(d, feature, windowsize, binsize) {
  # Get dataframe of feature signal for windowsize 
  # upstram and downstream of compartment switch from B to A
  
  # Create dataframe of 2Mb window around A->B boundary switches
  BtoA_feature <- data.frame()
  windowbins <- windowsize/binsize
  genomebins <- nrow(d)
  for (i in 1:(genomebins - 1)) {
    # Make sure 1Mb from start and end of genome
    if (i > windowbins & i < (genomebins - windowbins)) {
      # Make sure window size is within chromosome
      ichrom <- as.character(d[i, "chrom"])
      if (as.character(d[(i-windowbins),"chrom"]) ==  ichrom &
          as.character(d[(i+windowbins),"chrom"]) == ichrom) { 
        # Check if eigen switch
        ieigen <- as.numeric(d[i, "PCA_eigen1"])
        nexteigen <- as.numeric(d[i+1, "PCA_eigen1"])
        if ((!is.na(ieigen)) & (!is.na(nexteigen))) {
          if (ieigen < 0 & nexteigen > 0) {
            # Get data around 2Mb of eigen switch (boundary)
            boundary_feature <- d[(i-windowbins):(i+windowbins), feature]
            BtoA_feature <- rbind(BtoA_feature, boundary_feature)
          }
        }
      }
    } 
  }
  # Name columns in kb
  cols <- ncol(BtoA_feature)
  colnames(BtoA_feature) <- c((floor(cols/2):0)*-bs/1000, (1:floor(cols/2))*bs/1000)
  return(BtoA_feature)
}


boundary_plot <- function(t, pc1, switch) {
  pdf(paste(switch,"/compartment_boundaries_", switch, "_half-life_LOS_PC1_merge.pdf", sep=""), 
      height=6, width=6)
    par(mar = c(5, 4, 4, 4) + 0.1)
    mean_t <- apply(t, 2, mean, na.rm=TRUE)
    mean_pc1 <- apply(pc1, 2, mean, na.rm=TRUE)
    bp <- as.numeric(colnames(t))
    plot(bp, mean_t, type="l", col="dodgerblue",
         xlab=paste(switch, " compartment boundary (kb)"), 
         ylab="", ylim=rev(range(mean_t)), cex.axis=0.8)
    axis(2, col="dodgerblue", col.axis="dodgerblue", cex.axis=0.8, at=pretty(range(mean_t), 8))
    mtext(bquote("t"[1/2] ~ "(minutes)"), side=2, col="dodgerblue", line=2.5)
    par(new=TRUE)
    plot(bp, mean_pc1, type="l", col="black",
         axes=FALSE, bty="n", xlab="", ylab="")
    axis(side=4, at=pretty(range(mean_pc1), 5), cex.axis=0.8)
    mtext("EV1", side=4, line=2.5)
    
  dev.off()
}

bs <- get_bin_size(df[1,])
ws <- 1000000   # window upstream and downstream of boundary switch

AtoB_thalf <- get_AtoB_df(df, "half-life_LOS", ws, bs)
AtoB_PC1 <- get_AtoB_df(df, "PCA_eigen1", ws, bs)
BtoA_thalf <- get_BtoA_df(df, "half-life_LOS", ws, bs)
BtoA_PC1 <- get_BtoA_df(df, "PCA_eigen1", ws, bs)
boundary_plot(AtoB_thalf, AtoB_PC1, "A-B")
boundary_plot(BtoA_thalf, BtoA_PC1, "B-A")





