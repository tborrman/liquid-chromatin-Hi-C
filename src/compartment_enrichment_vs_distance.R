compartment_enrichment <- function(compartment_df, AB_df, distances) {
  ratio_enrich <- c()
  for (distance in distances) {
    comp_subset <- compartment_df[compartment_df$distance == distance,]
    AB_subset <- AB_df[AB_df$distance == distance,]
    ratio <- mean(comp_subset$cScore, na.rm=TRUE) / mean(AB_subset$cScore, na.rm=TRUE)
    ratio_enrich <- c(ratio_enrich, ratio)
  }
  return(ratio_enrich)
}

compartment_mean <- function(compartment_df, distances) {
  mean_list <- c()
  for (distance in distances) {
    comp_subset <- compartment_df[compartment_df$distance == distance,]
    mean_scalar <- mean(comp_subset$cScore, na.rm = TRUE)
    mean_list <- c(mean_list, mean_scalar)
  }
  return(mean_list)
}

wkdir <- getwd()
DN_df<- read.table(paste(wkdir, "/eigenI_eigenJ_DN_interactome.txt", sep=""), sep="\t", header=TRUE)
MN_df<- read.table(paste(wkdir, "/eigenI_eigenJ_MN_interactome.txt", sep=""), sep="\t", header=TRUE)

distance <- abs(DN_df$I - DN_df$J)
DN_df <- cbind(DN_df, distance)

distance <- abs(MN_df$I - MN_df$J)
MN_df <- cbind(MN_df, distance)


#Remove rows with NAs in eigen columns
DN_df <- DN_df[(!is.na(DN_df$eigen_I)) & (!is.na(DN_df$eigen_J)),]
MN_df <- MN_df[(!is.na(MN_df$eigen_I)) & (!is.na(MN_df$eigen_J)),]

DN_AA <- DN_df[DN_df$eigen_I > 0 & DN_df$eigen_J > 0,]
DN_BB <- DN_df[DN_df$eigen_I < 0 & DN_df$eigen_J < 0,]
DN_AB <- DN_df[(DN_df$eigen_I < 0 & DN_df$eigen_J > 0) | (DN_df$eigen_I > 0 & DN_df$eigen_J < 0),]

MN_AA <- MN_df[MN_df$eigen_I > 0 & MN_df$eigen_J > 0,]
MN_BB <- MN_df[MN_df$eigen_I < 0 & MN_df$eigen_J < 0,]
MN_AB <- MN_df[(MN_df$eigen_I < 0 & MN_df$eigen_J > 0) | (MN_df$eigen_I > 0 & MN_df$eigen_J < 0),]

distances <- seq(0,248500000, 500000)

# Mean plots
MN_AA_mean <- compartment_mean(MN_AA, distances)
MN_BB_mean <- compartment_mean(MN_BB, distances)
MN_AB_mean <- compartment_mean(MN_AB, distances)
DN_AA_mean <- compartment_mean(DN_AA, distances)
DN_BB_mean <- compartment_mean(DN_BB, distances)
DN_AB_mean <- compartment_mean(DN_AB, distances)

options(scipen = 1000)
png("compartment_means_vs_distance.png", width=2500, height=2500, res=300)
plot(distances/1000, MN_AA_mean, type='l', col="red", lwd=2, xlab= "Distance (kb)", ylab="Mean Interaction Frequency")
lines(distances/1000, MN_AB_mean, type='l', col="mediumorchid4", lwd=2)
legend("topright", c("Non-fragmented A-A", "Non-fragmented A-B"), lty=c(1,1), lwd=c(2,2), col=c("red", "mediumorchid4"))
dev.off()

png("compartment_log_means_vs_distance.png", width=4500, height=3000, res=300)
par(mfrow=c(1,2))
#MN
plot(distances/1000, log(MN_AA_mean), type='l', col="red", lwd=2, xlab= "Distance (kb)", ylab="Log Mean Interaction Frequency",
     ylim=c(-1,7), main="Non-Fragmented")
lines(distances/1000, log(MN_AB_mean), type='l', col="mediumorchid4", lwd=2)
lines(distances/1000, log(MN_BB_mean), type='l', col="blue", lwd=2)

#DN
plot(distances/1000, log(DN_AA_mean), type='l', col="red", lwd=2, xlab= "Distance (kb)", ylab="Log Mean Interaction Frequency",
     ylim=c(-1,7), main="Fragmented")
lines(distances/1000, log(DN_AB_mean), type='l', col="mediumorchid4", lwd=2)
lines(distances/1000, log(DN_BB_mean), type='l', col="blue", lwd=2)
legend("topright", c("A-A interactions", "B-B interactions", "A-B interactions" ), lty=c(1,1,1), lwd=c(2,2,2), col=c("red", "blue", "mediumorchid4"))
dev.off()





# Ratio of enrichment 
MN_AA_enrich <- compartment_enrichment(MN_AA, MN_AB, distances)
MN_BB_enrich <- compartment_enrichment(MN_BB, MN_AB, distances)
DN_AA_enrich <- compartment_enrichment(DN_AA, DN_AB, distances)
DN_BB_enrich <- compartment_enrichment(DN_BB, DN_AB, distances)


png("compartment_enrich_vs_distance.png", width=3000, height=2500, res=300)
plot(distances/1000, MN_AA_enrich, type='l', lty=1, col="red", lwd=2, xlab= "Distance (kb)", ylab="Ratio of enrichment over AB interactions")
lines(distances/1000, MN_BB_enrich, type='l',lty=1, col="blue", lwd=2)
lines(distances/1000, DN_AA_enrich, type='l', lty=3, col="red", lwd=2 )
lines(distances/1000, DN_BB_enrich, type='l', lty=3, col="blue", lwd=2)

legend("topleft", c("Non-fragmented A-A", "Non-fragmented B-B", "Fragmented A-A", "Fragmented B-B"), lty=c(1,1,3,3), lwd=c(2,2,2,2),
       col=c("red", "blue", "red", "blue"))
dev.off()



