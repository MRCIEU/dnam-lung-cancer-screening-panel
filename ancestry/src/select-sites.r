#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
#' For each ancestry, identify top 50 SNPs that
#' 1. are mQTLs
#' 2. among those, top 250 in fst GWAS
#' 3. among those, top 50 GoDMC mQTL effect 

sites.dir <- args[1]
output.filename <- args[2]

cols <- c("chr","pos","id","obs.ct","hudson.fst","coords","z.glm","p.glm","beta.godmc","cpg.godmc","pct.godmc")

filenames <- list.files(sites.dir, ".csv", full.names=T)
names(filenames) <- basename(sub(".csv", "", filenames))

sites <- lapply(names(filenames), function(ancestry) {
    sites <- fread(filenames[[ancestry]])
    sites <- sites[!is.na(sites$p.glm),]
    stopifnot(length(cols) == ncol(sites))
    colnames(sites) <- cols
    if (nrow(sites) > 250)
        sites <- sites[order(sites$hudson.fst,decreasing=T)[1:250],]
    if (nrow(sites) > 50)
        sites <- sites[order(abs(sites$beta.godmc),decreasing=T)[1:50],]
    sites$start <- sites$end <- sites$pos
    sites$cpg <- sites$cpg.godmc
    sites$ancestry <- ancestry
    sites
})
sites <- do.call(rbind, sites)

fwrite(sites, file=output.filename)

