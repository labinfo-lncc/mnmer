library (dplyr)

args = commandArgs(trailingOnly=TRUE)
out <- sub (".rda","", args[1])
print (paste0("metrics_mean_",out,".csv"))

load (args[1])

rtab <- dat %>% group_by (db, mn, algorithm) %>% summarise (homogeneity=mean(homogeneity), commpleteness=mean(commpleteness), vmeasure=mean(vmeasure), adjrand=mean(adjrand), adjmutual=mean(adjmutual), silhouette=mean(silhouette))
rtab <- as.data.frame (rtab)
rtab <- rtab[order(rtab$mn, decreasing = T),] 

write.csv (rtab, file=paste0("metrics_mean_",out,".csv"), row.names=F)

