size1 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz1-K562-R1/bed/HB-Dp4h-siz1-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size1_reads"))
size2 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz2-K562-R1/bed/HB-Dp4h-siz2-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size2_reads"))
size3 <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/size_analysis/HB-Dp4h-siz3-K562-R1/bed/HB-Dp4h-siz3-K562-R1_coverage_40kb_readnorm_copycorrect.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "size3_reads"))

sub <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/subcompartments/SNIPER_K562_subcompartments_40kb.bed",
                  sep="\t", header=TRUE)

df <- merge(size1, size2, by = c("chrom", "start", "end"))
df <- merge(df, size3, by = c("chrom", "start", "end"))
df <- merge(df, sub, by = c("chrom", "start", "end"))

df <- na.omit(df)

sub_slices <- function(d, cn) {
  # Get read numbers for each 
  # subcompartment for slice of pie chart
  # d: full dataframe
  # cn: int for column number of DNA size
  
  slices <- c()
  for (s in c("A1", "A2", "B1", "B2", "B3")) {
    slices <- c(slices, sum(d[d$sub == s, cn]))
  }
  return(slices)
}

size1_col <- 4
size2_col <- 5
size3_col <- 6


size1_slices <- sub_slices(df, size1_col)
size1_pct <- round(size1_slices/sum(size1_slices)*100,2)
size2_slices <- sub_slices(df, size2_col)
size2_pct <- round(size2_slices/sum(size2_slices)*100,2)
size3_slices <- sub_slices(df, size3_col)
size3_pct <- round(size3_slices/sum(size3_slices)*100,2)

png("subcompartment_pie_charts.png", width=4000, height=1500, res=500)
par(mfrow=c(1,3))
pie(size1_slices, 
    labels = paste(c("A1 ", "A2 ", "B1 ", "B2 ", "B3 "), size1_pct, rep("%",5),sep=""), 
    col= c("red", "hotpink", "purple", "cyan", "blue"), 
    main = "Size 1")
pie(size2_slices, 
    labels = paste(c("A1 ", "A2 ", "B1 ", "B2 ", "B3 "), size2_pct, rep("%",5),sep=""), 
    col= c("red", "hotpink", "purple", "cyan", "blue"), 
    main = "Size 2")
pie(size3_slices, 
    labels = paste(c("A1 ", "A2 ", "B1 ", "B2 ", "B3 "), size3_pct, rep("%",5),sep=""), 
    col= c("red", "hotpink", "purple", "cyan", "blue"), 
    main = "Size 3")
dev.off()


