#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
#' For each ancestry, identify top 50 SNPs that
#' 1. are mQTLs
#' 2. have logistic GWAS p < 5e-8
#' 3. among those, top 200 in fst GWAS
#' 4. among those, top 50 GoDMC mQTL effect 

sites.dir <- args[1]
output.filename <- args[2]

cols <- c("chr","pos","id","obs.ct","hudson.fst","coords","pct.fst","z.glm","p.glm","beta.godmc","cpg.godmc","pct.godmc")

filenames <- list.files(sites.dir, ".csv", full.names=T)
names(filenames) <- basename(sub(".csv", "", filenames))

panel.sites <- lapply(names(filenames), function(ancestry) {
    sites <- fread(filenames[[ancestry]])
    if (nrow(sites) > 200)
        sites <- sites[order(sites$hudson.fst,decreasing=T)[1:200],]
    if (nrow(sites) > 50)
        sites <- sites[order(abs(sites$beta.godmc),decreasing=T)[1:50],]
    stopifnot(length(cols) == ncol(sites))
    colnames(sites) <- cols
    sites$start <- sites$end <- sites$pos
    sites$cpg <- sites$cpg.godmc
    sites$ancestry <- ancestry
    sites
})
panel.sites <- do.call(rbind, panel.sites)

fwrite(panel.sites, file=output.filename)

