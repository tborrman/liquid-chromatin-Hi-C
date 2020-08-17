library(ggplot2)
df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/feature_analysis/C-40000/v6/feature_matrix_v6_40kb.txt", sep="\t",
                 header=TRUE, check.names=FALSE)

size_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/compartment_size/filter1000/eigen1_40kb.size.bedGraph", 
                      sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "eigen", "size"))

sub_df <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                     sep="\t", header=TRUE)

df <- cbind(df, size_df["size"], sub_df["sub"])

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
            if(d[i,"size"] > 1000000 & d[i+1, "size"] > 1000000) {
              # Get data around 2Mb of eigen switch (boundary)
              boundary_feature <- d[(i-windowbins):(i+windowbins), feature]
              if(feature == "sub") {
                boundary_feature <- as.character(boundary_feature)
                AtoB_feature <- rbind(AtoB_feature, boundary_feature, stringsAsFactors=FALSE)
              }
              else {
                AtoB_feature <- rbind(AtoB_feature, boundary_feature)
              }
            }
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
            if(d[i,"size"] > 1000000 & d[i+1, "size"] > 1000000) {
              # Get data around 2Mb of eigen switch (boundary)
              boundary_feature <- d[(i-windowbins):(i+windowbins), feature]
              BtoA_feature <- rbind(BtoA_feature, boundary_feature)
            }
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
  pdf(paste("compartment_boundaries_", switch, "_half-life_LOS_PC1_merge_sizefix.pdf", sep=""), 
      height=6, width=6)
    par(mar = c(5, 4, 4, 4) + 0.1)
    mean_t <- apply(t, 2, mean, na.rm=TRUE)
    mean_pc1 <- apply(pc1, 2, mean, na.rm=TRUE)
    bp <- as.numeric(colnames(t))
    plot(bp, mean_t, type="l", col="dodgerblue",
         xlab=paste(switch, " compartment boundary (kb)"), 
         ylab="", cex.axis=0.8)
    axis(2, col="dodgerblue", col.axis="dodgerblue", cex.axis=0.8, at=pretty(range(mean_t), 4))
    mtext(bquote("t"[1/2] ~ "(minutes)"), side=2, col="dodgerblue", line=2.5)
    par(new=TRUE)
    plot(bp, mean_pc1, type="l", col="black",
         axes=FALSE, bty="n", xlab="", ylab="", ylim=rev(range(mean_pc1)))
    axis(side=4, at=pretty(range(mean_pc1), 4), cex.axis=0.8)
    mtext("PC1", side=4, line=2.5)
    
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




# Plot stacked barplot of subcompartments at A->B boundaries
AtoB_sub <- get_AtoB_df(df, "sub", ws, bs)
windowaxis <- factor(rep(colnames(AtoB_sub), 5),
                     levels=colnames(AtoB_sub), ordered=TRUE)
subaxis <- factor(c(rep("A1", 51), rep("A2", 51), rep("B1", 51), 
                    rep("B2", 51), rep("B3", 51)), levels=c("A1", "A2", "B1", "B2", "B3"),
                    ordered=TRUE)

counts <- c()
for (i in 1:length(windowaxis)){
    w = windowaxis[i]
    s = subaxis[i]
    counts <- c(counts,sum(AtoB_sub[w] == as.character(s), na.rm = TRUE))
}
  
stack_df <- data.frame(windowaxis, subaxis, counts)
pdf("subcompartment_PC1_boundaries_stacked_barplot.pdf",
    height=7,width=12)

ggplot(stack_df, aes(fill=subaxis, y=counts, x=windowaxis)) + 
  geom_bar( stat="identity", position="fill") +
  xlab("A-B compartment boundary (kb)") + ylab("Subcompartment %") +
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
  scale_fill_manual(values=c("A1" = "#ED2024",
                              "A2" = "#F06BA8",
                              "B1" = "#8150A0",
                              "B2" = "#6FCCDD",
                              "B3" = "#3953A4"),
                     breaks=c("A1", "A2", "B1", "B2", "B3"))
dev.off()


