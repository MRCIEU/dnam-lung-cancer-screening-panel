#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output.filename <- args[1]
reduced.filename <- args[2]

library(data.table)

## load functions for manipulating GenomicRanges (specifically the 'reduce()' function)

source("src/granges.r")

## load illumina 450k/epic annotation

source("src/load-illumina-manifest-function.r")
manifest <- load.illumina.manifest(".")

## sources obtained from Illumina bead chips

illumina.sources <- read.csv(text="id,filename
age,hannum-model.csv
alcohol-consumption,dnamalci-sites.csv
alcohol-consumption,chamberlain-sites.csv
alcohol-consumption,alc-ewas-sites.csv
alcohol-consumption,mccartney-sites.csv
ancestry,output/sites.csv
bmi,bmi-model.csv
breast-cancer,xu-sites.csv
cadmium,cadmium-sites.csv
copd,important-sites.csv
cotinine,discovery-sites.csv
crp,crp-sites.csv
crp,reduced-hillary-model.csv
dunedin-pace,pace-model.csv
dunedin-poam38,dunedinpoam38-sites.csv
educational-attainment,sites.csv
egfr,breeze.csv
egfr,ammous.csv
episcores,episcore-sites.csv
hdl,mccartney-sites.csv
il6,stevenson-sites.csv
lead,lead-model.csv
lung-cancer,battram-sites.csv
lung-cancer,zhao-sites.csv
prostate-cancer,sites.csv
prostate-cancer-t-cells,sites.csv
smoking-cessation,guida-sites.csv
smoking-cessation,joehanes-sites.csv
smoking-cessation,mccartney-sites.csv
smoking-former,andrayas-sites.csv
smoking-status,maas-sites.csv", stringsAsFactors=F)

## sources that supply genomic regions

region.sources <- read.csv(text="id,filename
blood-cell-types,regions-hg19.csv
colorectal-cancer,regions-hg19.csv", stringsAsFactors=F)

## load illumina sources

illumina.sites <- lapply(
    1:nrow(illumina.sources),
    function(i) {
        id <- illumina.sources$id[i]
        filename <- file.path(id, illumina.sources$filename[i])
        cat("loading", filename, "\n")
        dat <- read.csv(filename,stringsAsFactors=F)
        data.frame(source=id, details=dat$cpg)
    })
illumina.sites <- do.call(rbind, illumina.sites)

## add genomic coordinates for illumina sites

idx <- match(illumina.sites$details, manifest$IlmnID)
stopifnot(!any(is.na(idx)))

illumina.sites <- data.frame(
    illumina.sites,
    chr=manifest$CHR[idx],
    start=manifest$MAPINFO[idx],
    end=manifest$MAPINFO[idx],
    stringsAsFactors=F)

## load regions 

regions <- lapply(
    1:nrow(region.sources),
    function(i) {
        id <- region.sources$id[i]
        filename <- file.path(id, region.sources$filename[i])
        cat("loading", filename, "\n")
        dat <- read.csv(filename, stringsAsFactors=F)
        data.frame(
            source=id,
            dat[,c("chr","start","end","details")],
            stringsAsFactors=F)
    })
regions <- do.call(rbind, regions)

## merge all regions into a single panel and save it

panel <- rbind(
    illumina.sites,
    regions)

fwrite(panel, file=output.filename)

## merge overlapping regions
panel.gr <- to.granges(panel$chr,panel$start,panel$end)
panel.reduced <- from.granges(reduce(panel.gr))

mapping <- lapply(1:nrow(panel.reduced), function(i) {
    with(panel.reduced, {
        which(
            panel$chr==chr[i]
            & (panel$start >= start[i] & panel$start <= end[i]
               | panel$end >= start[i] & panel$end <= end[i]))
    })
})

panel.reduced$source <- sapply(mapping, function(idx) paste(unique(panel$source[idx]),collapse="/"))

panel.reduced$details <- sapply(mapping, function(idx) paste(unique(panel$details[idx]),collapse="/"))

fwrite(panel.reduced, file=reduced.filename)

