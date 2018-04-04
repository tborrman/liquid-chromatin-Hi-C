hn <- read.table("hindIII_coverage_sorted_500kb_GMM.txt", sep="\t", header=TRUE)
dn <- read.table("dpnII_coverage_sorted_500kb_R1_GMM.txt", sep="\t", header=TRUE)
h40 <- read.table("hindIII_coverage_sorted_40kb_GMM.txt", sep="\t", header=TRUE)
d40 <- read.table("dpnII_coverage_sorted_40kb_R1_GMM.txt", sep="\t", header=TRUE)

png("dpnII_vs_hindIII_500kb.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(hn$reads, dn$reads, pch=21, ylim=c(1000,5000), xlim=c(2000,8000), bg="dodgerblue",
     xlab="HindIII reads (500kb)", ylab= "DpnII reads R1 (500kb)")
text(6500,1500, bquote("r"[s]~ "="~  
                         .(round(cor(hn$reads, dn$reads, method="spearman"),2))), 
     cex=1.5)
dev.off()

png("dpnII_vs_hindIII_40kb.png", height=2000, width=2000, res=250)
par(mar = c(5, 5, 4, 2) + 0.1)
plot(h40$reads, d40$reads, pch=21, ylim=c(0,500), xlim=c(0,800), bg="dodgerblue",
     xlab="HindIII reads (40kb)", ylab= "DpnII reads R1 (40kb)")
text(600, 75, bquote("r"[s]~ "="~  
                         .(round(cor(h40$reads, d40$reads, method="spearman"),2))), 
     cex=1.5)
dev.off()

hn_copy <- hn$copy
dn_copy <- dn$copy

nas <-  sum(is.na(hn_copy) | is.na(dn_copy))
not_eq <- sum(hn_copy!= dn_copy, na.rm=TRUE)

slices <- c(length(hn$copy) - (nas + not_eq), not_eq, nas)
pct <- round(slices/sum(slices)*100)
mylabs <- c("Identical", "Different", "NA")

png("DpnII_HindIII_pie_500kb.png", width=3000, height=3000, res=500)
pie(slices, labels= paste(mylabs, rep(" ",length(mylabs)), pct, rep("%", length(mylabs)), sep=""),
    main="HindIII and DpnII copy number states (500kb)",
    col=c("dodgerblue", "red", "grey50"))
dev.off()

hn_copy <- h40$copy
dn_copy <- d40$copy

nas <-  sum(is.na(hn_copy) | is.na(dn_copy))
not_eq <- sum(hn_copy!= dn_copy, na.rm=TRUE)

slices <- c(length(hn_copy) - (nas + not_eq), not_eq, nas)
pct <- round(slices/sum(slices)*100)
mylabs <- c("Identical", "Different", "NA")

png("DpnII_HindIII_pie_40kb.png", width=3000, height=3000, res=500)
pie(slices, labels= paste(mylabs, rep(" ",length(mylabs)), pct, rep("%", length(mylabs)), sep=""),
    main="HindIII and DpnII copy number states (40kb)",
    col=c("dodgerblue", "red", "grey50"))
dev.off()
