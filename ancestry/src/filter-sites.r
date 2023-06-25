#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
#' For each ancestry, identify top SNPs
#' that are mQTLs and for which
#' GWAS fst > 95th (among mqtls)

ancestry <- args[1]
glm.dir <- args[2]
fst.dir <- args[3]
godmc.file <- args[4]
dir.create(output.dir <- args[5])

cat("filter-sites.r", ancestry, glm.dir, fst.dir, godmc.file, output.dir, "\n")

## mQTL summary statistics (GoDMC)
godmc <- fread(godmc.file)
godmc.pct <- ecdf(abs(godmc$beta))
godmc$beta.pct <- godmc.pct(abs(godmc$beta))

## GWAS output files
fst.files <- list.files(
    fst.dir,
    pattern=paste0(ancestry, "-chr[0-9]+.*.fst.var$"),
    full.names=T)

glm.files <- list.files(
    glm.dir,
    pattern=paste0(ancestry, "-chr[0-9]+.*.glm.logistic$"),
    full.names=T)

fst.files <- sort(fst.files)
glm.files <- sort(glm.files)

stopifnot(length(fst.files) == length(glm.files))

load.gwas <- function(filename) {
    cat(date(), "loading", filename, "...\n")
    dat <- as.data.frame(fread(filename))
    dat$coords <- paste(dat[["#CHROM"]], dat[["POS"]], sep="_")
    dat <- dat[order(dat[["#CHROM"]], as.integer(dat[["POS"]])),]
    dat
}

dat <- lapply(1:length(fst.files), function(i) {
    ## summary stats for GWAS (Fst)
    dat <- load.gwas(fst.files[i])
    
    ## summary stats for GWAS (logistic regression)
    dat.glm <- load.gwas(glm.files[i])

    ## add logistic GWAS summary statistics 
    stopifnot(identical(dat$coords, dat.glm$coords))
    dat$z.glm <- dat.glm$Z_STAT
    dat$p.glm <- dat.glm$P

    ## add mQTL statistics
    idx <- match(dat$coords, godmc$coords)
    dat$beta.godmc <- godmc$beta[idx]
    dat$cpg.godmc <- godmc$cpg[idx]
    dat$pct.godmc <- godmc$beta.pct[idx]

    ## keep only that are also mqtls
    dat <- dat[!is.na(dat$beta.godmc),]

    ## keep those in the top 97.5th percentile of GWAS fst
    dat[which(dat$HUDSON_FST > quantile(dat$HUDSON_FST,probs=0.975,na.rm=T)),]
})

dat <- do.call(rbind, dat)

filename <- file.path(output.dir, paste0(ancestry, ".csv"))

fwrite(dat, file=filename)


