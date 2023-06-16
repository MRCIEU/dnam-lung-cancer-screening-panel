sources <- c(
    "age",
    "ancestry",
    "blood-cell-types",
    "bmi",
    "breast-cancer",
    "cadmium",
    "colorectal-cancer",
    "copd",
    "cotinine",
    "dunedin-pace",
    "educational-attainment",
    "episcores",
    "lung-cancer",
    "prostate-cancer",
    "prostate-cancer-t-cells",
    "smoking-cessation",
    "smoking-former",
    "smoking-status")

annotate.450k <- function(sites, manifest=load.450k.manifest()) {
    manifest[match(sites, manifest$IlmnID),]
}

illumina2regions.r 
age/hannum-model.csv  
ancestry/ancestry-sites.csv
bmi/bmi-model.csv 
breast-cancer/xu-sites.csv
cadmium/cadmium-sites.csv
copd/important-sites.csv
cotinine/discovery-sites.csv
dunedin-pace/pace-model.csv
educational-attainment/sites.csv
episcores/episcore-sites.csv
lead/lead-model.csv
lung-cancer/battram-sites.csv
lung-cancer/zhao-sites.csv
prostate-cancer/sites.csv
prostate-cancer-t-cells/sites.csv
smoking-cessation/guida-sites.csv
smoking-cessation/joehanes-sites.csv
smoking-cessation/mccartney-sites.csv
smoking-former/andrayas-sites.csv
smoking-status/maas-sites.csv

liftover.r hg19 hg38 
blood-cell-types/regions-hg19.csv
colorectal-cancer/regions-hg19.csv


files <- file.path(sources, "regions.csv")

cols <- c("chr","start","end","source","details")
panel <- lapply(sources, function(src) {
    filename <- file.path(src, "regions.csv")
    cat("loading ", filename, "\n")
    stopifnot(file.exists(filename))                
    regions <- read.csv(filename)
    regions$source <- src
    stopifnot(cols %in% colnames(regions))
    regions[,cols]
})
panel <- do.call(rbind, panel)

write.csv(panel, file="panel.csv", row.names=F)
