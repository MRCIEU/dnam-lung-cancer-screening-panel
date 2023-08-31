## These packages, especially tidyverse, require several developmental linux libraries: libssl-dev, libcurl4-openssl-dev, libfontconfig-dev, libxml2-dev, libfontconfig1-dev, libharfbuzz-dev, librust-harfbuzz-sys-dev, libfribidi-dev, libtiff-dev

pkgs <- list(
    cran=c("remotes","data.table","tidyverse","vcfR",
           "matrixStats","readxl","BiocManager"),
    git=c("danbelsky/DunedinPACE", "perishky/meffonym"),
    bioc=c("GenomicRanges","MungeSumstats","liftOver"))

is.installed <- function(pkg)
    require(pkg,character.only=T,quietly=T)

for (pkg in pkgs$cran)
    if (!is.installed(pkg))
        install.packages(pkg)

for (pkg in pkgs$git)
    if (!is.installed(basename(pkg)))
        remotes::install_github(pkg)

for (pkg in pkgs$bioc)
    if (!is.installed(pkg))
        BiocManager::install(pkg)

for (pkg in basename(unlist(pkgs)))
    cat("installed =", is.installed(pkg), pkg, "\n")


