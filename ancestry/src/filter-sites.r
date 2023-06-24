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

load.gwas <- function(path, ancestry, model=c("glm","fst")) {
  files <- list.files(
      path,
      pattern=paste0(ancestry, "-chr[0-9]+[.]+.+", model),
      full.names=T)
  stopifnot(length(files)>0)
  dat <- do.call(rbind, lapply(files, fread))
  dat$coords <- paste(dat[["#CHROM"]], dat[["POS"]], sep="_")
  dat <- dat[order(dat[["#CHROM"]], as.integer(dat[["POS"]])),]
  dat
}

## mQTL summary statistics (GoDMC)
godmc <- fread(godmc.file)
godmc.pct <- ecdf(abs(godmc$beta))
godmc$beta.pct <- godmc.pct(abs(godmc$beta))

filename <- file.path(output.dir, paste0(ancestry, ".csv"))
    
## summary stats for GWAS (Fst)
dat <- load.gwas(fst.dir, ancestry, "fst")
dat$pct.fst <- ecdf(dat$HUDSON_FST)(dat$HUDSON_FST)

## summary stats for GWAS (logistic regression)
dat.glm <- load.gwas(glm.dir, ancestry, "glm")

## add logistic GWAS summary statistics 
stopifnot(identical(dat$coords, dat.glm$coords))
dat$z.glm <- dat.glm$Z_STAT
dat$p.glm <- dat.glm$P

## add mQTL statistics
idx <- match(dat$coords, godmc$coords)
dat$beta.godmc <- godmc$beta[idx]
dat$cpg.godmc <- godmc$cpg[idx]
dat$pct.godmc <- godmc$beta.pct[idx]

## top 200 hudson fst with logistic p < 5e-8
iseligible <- dat$p.glm < 5e-8 & !is.na(dat$beta.godmc)
if (sum(iseligible) >= 200)
    threshold <- sort(dat$HUDSON_FST[iseligible],decreasing=T)[200]
else
    threshold <- min(dat$HUDSON_FST[iseligible],na.rm=T)
dat <- dat[iseligible & dat$HUDSON_FST >= threshold,]

dat <- dat[order(abs(dat$beta.godmc),decreasing=T),]

fwrite(dat, file=filename)


