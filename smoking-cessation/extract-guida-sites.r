library(readxl)
guida <- readxl::read_xlsx("guida_supplement.xlsx")
colnames(guida) <- c("cpg","chr","pos","gene","gene.region","cgi.region","direction","time")
guida <- guida[guida$time %in% c(10,15,20),]
table(guida$time)
#10 15 20 
#36 33 27 
write.csv(guida, file="guida-sites.csv", row.names=F)