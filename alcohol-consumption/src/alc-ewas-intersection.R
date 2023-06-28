library(data.table)
setwd("alcohol-consumption/")

## Read in supplementary data pertaining to top sites from Dugue et al. and
## Lohoff et al. restrict to top 100 CpGs for each (already ordered by P-value)
dugue <- fread("dugue_supplement_s3.csv", nrows = 100)
lohoff <- fread("lohoff_supplement_s2.csv", nrows = 100)

## Take intersection of top 100 CpG sites from each EWAS
cpgs <- intersect(dugue$CpG, lohoff$ID)

## Effect sizes and pvals in both EWAS for intersected sites
dugue[dugue$CpG %in% cpgs,8:9]
lohoff[lohoff$ID %in% cpgs,2:3]

## Write csv of intersecting CpGs 
write.csv(as.data.frame(cpgs), file = "alc-ewas-sites.csv", row.names = F, quote = F)