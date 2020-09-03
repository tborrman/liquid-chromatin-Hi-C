
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
  colnames(AtoB_feature) <- c(((25:0)*-40), ((1:25)*40))
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
  colnames(BtoA_feature) <- c(((25:0)*-40), ((1:25)*40))
  return(BtoA_feature)
}


boundary_plot <- function(d, feature, switch, mycolor) {
  png(paste(switch,"/compartment_boundaries_", switch, "_", feature, ".png", sep=""), 
      height=2000, width=4000, res=300)
  par(mfrow=c(1,2))
  boxplot(d,  col=mycolor,
          xlab=paste(switch, " compartment boundary (kb)"),
          ylab=feature)
  mean_feature <- apply(d, 2, mean, na.rm=TRUE)
  bp <- as.numeric(colnames(d))
  plot(bp, mean_feature, type="l", col=mycolor,
       xlab=paste(switch, " compartment boundary (kb)"), 
       ylab=feature)
  dev.off()
}

bs <- get_bin_size(df[1,])
ws <- 1000000   # window upstream and downstream of boundary switch

for (f in colnames(df)[4:ncol(df)]) {
  AtoB <- get_AtoB_df(df, f, ws, bs)
  boundary_plot(AtoB, f, "A-B", "dodgerblue")
  BtoA <- get_BtoA_df(df, f, ws, bs)
  boundary_plot(BtoA, f, "B-A", "dodgerblue")
}




