source("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/digest-Hi-C/my_functions.R")
thalf <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/filter1000/timecourse1/half-life/half-life_exponential_40kb_removed_outliers_range6Mb_filter1000_timecourse1.bedGraph",
                    sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "hl"))
dseq <- read.table("C:/Users/tyler/Dropbox (UMass Medical School)/digest_092718/DpnII-seq/timecourse/copy_correct_coverage/40kb/HBDpSeqK562-DN1hR1_S5_L002_copy_correct_coverage_40kb.bed",
                   sep="\t", header=FALSE, col.names=c("chrom", "start", "end", "DpnIIseq"))

get_residuals <- function(x, p) {
  Y <- as.numeric(as.character(x))
  X <- as.numeric(as.character(p))
  m <- lm(Y~X)
  r <- m$residuals
  return(r)
}



df <- cbind(thalf, dseq["DpnIIseq"])
df_clean <- na.omit(df)

hl_residuals <- get_residuals(df_clean$hl, df_clean$DpnIIseq)

b <- cbind(df_clean[1:3], hl_residuals)
m <- merge(df[1:3], b, by=c("chrom", "start", "end"), all=TRUE)
om <- OrderChromStart(m)

write.table(om, "half-life_filter1000_timecourse1_residuals_from_1hDpnIIseq_40kb.bedGraph", sep="\t",
            quote=FALSE, row.names=FALSE, col.names=FALSE)





