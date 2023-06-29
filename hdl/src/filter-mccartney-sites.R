library(data.table)

## Read in McCartney et al. (https://doi.org/10.1186/s13059-018-1514-1) alcohol
## predictor CpGs and weights from Supplementary Table S3
mccartney <- fread("mccartney_supplement_s6.csv")

## Restrict sites to those with an effect of >=0.2
mccartney_restricted <- mccartney[abs(mccartney$Beta)>=0.2,]

## Change column names and write csv
colnames(mccartney_restricted) <- c("cpg", "coefficient")
write.csv(mccartney_restricted, file = "mccartney-sites.csv", quote = F, row.names = F)
