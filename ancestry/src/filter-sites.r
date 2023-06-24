#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
#' For each ancestry, identify top SNPs for which
#' 1. logistic GWAS p < 5e-8
#' 2. among those, top 200 in fst GWAS
#' 3. among those, top 50 GoDMC mQTL effect 

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

    ## keep only with logistic p < 5e-8
    dat <- dat[which(dat$p.glm < 5e-8),]

    ## add mQTL statistics
    idx <- match(dat$coords, godmc$coords)
    dat$beta.godmc <- godmc$beta[idx]
    dat$cpg.godmc <- godmc$cpg[idx]
    dat$pct.godmc <- godmc$beta.pct[idx]

    ## keep only that are also mqtls
    dat[!is.na(dat$beta.godmc),]
})

dat <- do.call(rbind, dat)

## keep the top 200 by FST
if (nrow(dat) > 200)
    dat <- dat[order(dat$HUDSON_FST,decreasing=T)[1:200],]

## sort by mQTL effect
dat <- dat[order(abs(dat$beta.godmc),decreasing=T),]

filename <- file.path(output.dir, paste0(ancestry, ".csv"))

fwrite(dat, file=filename)


