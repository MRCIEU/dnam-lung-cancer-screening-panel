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
