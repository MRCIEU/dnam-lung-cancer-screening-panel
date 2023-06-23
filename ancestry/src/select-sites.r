#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
#' For each ancestry, identify top SNPs for which
#' 1. logistic GWAS p < 5e-8
#' 2. among those, top 200 in fst GWAS
#' 3. among those, top 50 GoDMC mQTL effect 
#'
#' (about 5 minutes per ancestry, total is 2-2.5 hours runtime)

sites.dir <- args[1]
output.filename <- args[2]

cols <- c("chr","pos","id","obs.ct","hudson.fst","coords","pct.fst","z.glm","p.glm","beta.godmc","cpg.godmc","pct.godmc")

filenames <- list.files(sites.dir, ".csv", full.names=T)
names(filenames) <- basename(sub(".csv", "", filenames))

panel.sites <- lapply(names(filenames), function(ancestry) {
    sites <- fread(filenames[[ancestry]])
    sites <- sites[1:min(50,nrow(sites)),]
    stopifnot(length(cols) == ncol(sites))
    colnames(sites) <- cols
    sites$start <- sites$end <- sites$pos
    sites$cpg <- sites$cpg.godmc
    sites$ancestry <- ancestry
    sites
})
panel.sites <- do.call(rbind, panel.sites)

fwrite(panel.sites, file=output.filename)

