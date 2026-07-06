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

illumina.sources <- as.data.frame(fread("sources.csv"))

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
        data.frame(
			source=id, 
			pmid=illumina.sources$pmid[i],
			type=illumina.sources$type[i],
			details=dat$cpg)
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
            pmid=NA,
            type=NA,
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

