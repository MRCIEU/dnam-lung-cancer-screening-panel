#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output.filename <- args[1]
reduced.filename <- args[2]
dir.create(annot.dir <- args[3], showWarnings=F)

library(data.table)

## load illumina 450k/epic annotation

source("src/load-illumina-manifest-function.r")
manifest <- load.illumina.manifest(annot.dir)

## load chain file for converting hg19 to hg38 genomic coordinates

## BiocManager::install("liftOver")
library(liftOver)
source("src/granges.r")
source("src/liftover-function.r")
chain.file <- retrieve.chain.file(from="hg19",to="hg38",path=annot.dir)

## sources obtained from Illumina bead chips

illumina.sources <- read.csv(text="id,filename
age,hannum-model.csv
ancestry,output/panel-sites/panel-sites.csv
bmi,bmi-model.csv
breast-cancer,xu-sites.csv
cadmium,cadmium-sites.csv
copd,important-sites.csv
cotinine,discovery-sites.csv
dunedin-pace,pace-model.csv
educational-attainment,sites.csv
episcores,episcore-sites.csv
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

## sources with hg19 coordinates

hg19.sources <- read.csv(text="id,filename
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

## load regions with hg19 coordinates

hg19.regions <- lapply(
    1:nrow(hg19.sources),
    function(i) {
        id <- hg19.sources$id[i]
        filename <- file.path(id, hg19.sources$filename[i])
        cat("loading", filename, "\n")
        dat <- read.csv(filename, stringsAsFactors=F)
        data.frame(
            source=id,
            dat[,c("chr","start","end","details")],
            stringsAsFactors=F)
    })
hg19.regions <- do.call(rbind, hg19.regions)

## convert to hg38 coordinates

illumina.sites$chr <- paste0("chr", illumina.sites$chr)
illumina.sites.hg38 <- cbind(
    liftover(illumina.sites,chain.file),
    illumina.sites[,c("source","details")])

hg19.regions.hg38 <- cbind(
    liftover(hg19.regions, chain.file),
    hg19.regions[,c("source","details")])

## fix one region

idx <- which(is.na(illumina.sites.hg38$chr))
stopifnot(
    length(idx)==1 &
    illumina.sites.hg38$details[idx]=="cg23997508")
illumina.sites.hg38$chr[idx] <- "chr22"
illumina.sites.hg38$start[idx] <- 12097170
illumina.sites.hg38$end[idx] <- 12097170

## merge all regions into a single panel and save it

panel <- rbind(
    illumina.sites.hg38,
    hg19.regions.hg38)

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



## check the panels

nrow(panel)
## 3880

nrow(panel.reduced)
## 2754

quantile(panel$end-panel$start)
##  0%  25%  50%  75% 100% 
##   0    0    0    0 1457

quantile(setdiff(panel$end-panel$start,0))
##     0%     25%     50%     75%    100% 
##  23.00  174.25  281.50  417.50 1457.00 

quantile(setdiff(panel.reduced$end-panel.reduced$start,0))
##     0%     25%     50%     75%    100% 
##  23.00  174.25  281.50  417.50 1457.00 

table(sapply(mapping,length))
##  1    2    3    4    5    6    7    8    9   10   11   12   13   14   16 
##2258  253  106   54   33   16    8    8    4    5    2    2    3    1    1

panel.reduced$source[which(sapply(mapping,length) > 10)]
## [1] "ancestry" 
## [2] "ancestry" 
## [3] "ancestry" 
## [4] "bmi/cadmium/cotinine/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"
## [5] "ancestry" 
## [6] "ancestry" 
## [7] "ancestry" 
## [8] "ancestry" 
## [9] "ancestry"
## (indicates that some CpG sites are being used for multiple ancestries)

table(sapply(strsplit(panel.reduced$source, "/"), length))
##   1    2    3    4    5    6    7    8 
##2657   77   11    3    2    1    1    2 

panel.reduced[which(sapply(strsplit(panel.reduced$source, "/"),length) > 5),"source"]
## [1] "cadmium/cotinine/dunedin-pace/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"
## [2] "cadmium/cotinine/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"             
## [3] "bmi/cadmium/cotinine/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"         
## [4] "cadmium/cotinine/episcores/lung-cancer/smoking-cessation/smoking-status"

table(sapply(strsplit(panel.reduced$details, "/"), length))
##   1 
##2754 

