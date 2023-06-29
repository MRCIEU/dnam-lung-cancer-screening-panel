library(data.table)

## Read in supplementary data pertaining to top sites from Dugue et al. and
## Lohoff et al. for each (already ordered by P-value)
dugue <- as.data.frame(fread("dugue_supplement_s3.csv"))
lohoff <- as.data.frame(fread("lohoff_supplement_s2.csv"))

dugue <- dugue[dugue[,"Any at 10-7"] == "yes",] ## report in abstract

## Want intersection of about top 100 CpG sites from each EWAS
cpgs <- intersect(dugue$CpG, lohoff$ID)

nrow(dugue) ## 1414
nrow(lohoff) ## 2506
length(cpgs) ## 217

## Need to reduce this number to half
dugue$p <- pmin(dugue[,"P (week)"], dugue[,"P (decade)"], dugue[,"P (life)"])

dugue <- dugue[dugue$p < 1e-10,]
lohoff <- lohoff[lohoff$P.value < 1e-10,]

cpgs <- intersect(dugue$CpG, lohoff$ID)

nrow(dugue) ## 349
nrow(lohoff) ## 1030
length(cpgs) ## 80

## We'll take the top 80

## Effect sizes and pvals in both EWAS for intersected sites
dugue[dugue$CpG %in% cpgs,c("EST (week","P (week)","EST (decade)","P (decade)","EST (life)","P (life)")]
lohoff[lohoff$ID %in% cpgs,c("Zscore","P.value")]

## Write csv of intersecting CpGs 
write.csv(as.data.frame(cpgs), file = "alc-ewas-sites.csv", row.names = F, quote = F)
