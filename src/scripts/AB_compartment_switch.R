library(ggplot2)
NT <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/subcompartments/compartments/DL-20200103-LCHCSeq-K562-PIIBNT-Mock-R1-filter1000__hg19__genome__C-100000-iced__allchr.zScore.eigen1.bedGraph",
                 sep="\t", header=FALSE, col.names= c("chrom", "start", "end", "eigen_NT"))
TD <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/PIIB/subcompartments/compartments/DL-20200103-LCHCSeq-K562-PIIBTD-Mock-R1-filter1000__hg19__genome__C-100000-iced__allchr.zScore.eigen1.bedGraph",
                 sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "eigen_TD"))
subc <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/K562_track_hg19.bed",
                   sep="\t", skip=1, col.names=c("chrom", "start", "end", "sub", "a", "b", "c", "d", "e"))
subc <- subc[1:4]
subc$start <- subc$start + 1
df <- merge(NT, TD, by=c("chrom", "start", "end"))
df <- merge(df, subc, by=c("chrom", "start", "end"))

keep_chr <- function(d, chroms) {
  # Remove chromosomes that couldn't accurately
  # call compartments due to resolution
  # Args
  #   d: dataframe containing eigenvectors for Not Treated and Pol II Block
  #   c: vector of character chroms to keep
  # Return
  #   r: dataframe with chromosomes removed
  r <- data.frame()
  for (chrom in chroms) {
    d_chrom <- d[d$chrom == chrom,]
    r <- rbind(r, d_chrom)
  }
  return(r)
}

get_barplot_df <- function(d, s) {
  # Plot stacked barplot of A, B, NA compartment percentages
  # for Not Treated and Pol II Block samples for given subcompartment
  # called in K562 by SNIPER
  # Args
  #   d: dataframe containing eigenvectors for Not Treated and Pol II Block
  #   s: character for subcompartment
  
  dsub <- d[d$sub == s,]
  
 
  A_NT <- sum(dsub$eigen_NT > 0, na.rm=TRUE)
  B_NT <- sum(dsub$eigen_NT < 0, na.rm =TRUE)
  NA_NT <- sum(is.na(dsub$eigen_NT))
  A_TD <- sum(dsub$eigen_TD > 0, na.rm=TRUE)
  B_TD <- sum(dsub$eigen_TD < 0, na.rm =TRUE)
  NA_TD <- sum(is.na(dsub$eigen_TD))
  counts <- c(A_NT, B_NT, NA_NT, A_TD, B_TD, NA_TD)
  treatment <- c(rep("Not_Treated", 3), rep("PolII_Block", 3))
  compartment <- c("A", "B", "NA", "A", "B", "NA")
  bardf <- data.frame(counts, treatment, compartment)
  return(bardf)
}

plot_barplot <- function(sdf, sc) {
  # Plot stacked barplot giving percentage of A and B compartments
  # in Treated and Pol II Block samples for given SNIPER subcompartment bins
  # Args
  #    sdf: dataframe in stacked barplot format for ggplot
  #    sc: character of SNIPER subcompartment
  
  png(paste(sc, "_AB_stacked_barplot.png", sep=""), height=2000, width=1200, res=300)
    print(ggplot(sdf, aes(fill=compartment, y=counts, x=treatment)) + 
      geom_bar( stat="identity", position="fill") +
      xlab("Treatment") + ylab("Compartment %") +
      theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
            panel.background = element_blank(), axis.line = element_line(colour = "black")) +
      #axis.text.x=element_text(angle=90,hjust=1,vjust=0.5)) +
      scale_fill_manual(values=c("A" = "red", "B"= "blue", "NA"="grey"),
                        breaks=c("A", "B", "NA")) +
      ggtitle(paste(sc, "SNIPER subcompartment bins")))
  dev.off()
}

rdf <- keep_chr(df, c("chr11", "chr14", "chr15", "chr16", "chr17", "chr18", "chr19", "chr20", "chr21"))

subs <- c("A1", "A2", "B1", "B2", "B3")

for (i in subs) {
 stack_df <- get_barplot_df(rdf, i)
 plot_barplot(stack_df, i)
}

