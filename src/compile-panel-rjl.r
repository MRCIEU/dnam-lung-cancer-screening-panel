library(data.table)

## load illumina 450k/epic annotation (function doesn't work for me - commandArgs() pulls up NAs)
## have used slightly different code to do the same thing - un-comment lines 7,8 and 11 if running first time

## Download EPIC manifest
#online_file <- "https://webdata.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip"
#download.file(url = online_file, destfile = "epic_manifest.zip")

## Unzip EPIC manifest
#unzip("epic_manifest.zip")

## Read in EPIC manifest
manifest_epic <- fread("infinium-methylationepic-v-1-0-b5-manifest-file.csv", skip = 7)
manifest_epic$chip <- "epic"
manifest_epic <- as.data.frame(manifest_epic)

## Read in 450K manifest
manifest_450K <- fread("https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv", skip = 7)
manifest_450K$chip <- "450k"
manifest_450K <- as.data.frame(manifest_450K)

## Keep consistent column names and rbind
cols <- intersect(colnames(manifest_450K), colnames(manifest_epic))

manifest <- rbind(manifest_450K[,cols], manifest_epic[,cols])

## load chain file for converting hg19 to hg38 genomic coordinates

## BiocManager::install("liftOver")
library(liftOver)
#source("granges.r")

to.granges <- function(chr,start,end) {
  txt <- paste0(chr,":",start,"-",end)
  as(txt, "GRanges")
}

from.granges <- function(gr) {
  if ("CompressedGRangesList" %in% class(gr)) {
    len <- sapply(gr, length)
    if (any(len > 1)) {
      for (i in which(len > 1)) {
        chr <- as.character(seqnames(gr[[i]]))[1]
        start <- min(as.integer(start(gr[[i]])))
        end <- max(as.integer(end(gr[[i]])))
        gr[[i]] <- to.granges(chr,start,end)
      }
      gr <- unlist(gr)
    }
  }
  data.frame(
    chr=as.character(seqnames(gr)),
    start=as.integer(start(gr)),
    end=as.integer(end(gr)),
    stringsAsFactors=F)
  
}


# source("src/liftover-function.r")

retrieve.chain.file <- function(from,to,path=".") {
  capitalize <- function(str)
    paste0(toupper(substring(str,1,1)), tolower(substring(str,2)))
  from <- tolower(from)
  to <- capitalize(to)
  url <- paste0(
    "https://hgdownload.soe.ucsc.edu/goldenPath/",
    from,
    "/liftOver/",
    from,"To",to,".over.chain.gz")
  chain.file.gz <- file.path(path, basename(url))
  chain.file <- sub(".gz$", "", chain.file.gz)
  if (!file.exists(chain.file)) {
    download.file(
      url,
      destfile=chain.file.gz)
    system(paste("gunzip", chain.file.gz))
  }
  chain.file
}

liftover <- function(coords, chain.file) {
  stopifnot(all(c("chr","start","end") %in% colnames(coords)))    
  chain <- import.chain(chain.file)
  coords.from <- with(coords, to.granges(chr,start,end))
  coords.to <- liftOver(coords.from, chain)
  from.granges(coords.to)
}

chain.file <- retrieve.chain.file(from="hg19",to="hg38",path="../")

## sources obtained from Illumina bead chips

illumina.sources <- read.csv(text="id,filename
age,hannum-model.csv
ancestry,output/sites.csv
bmi,bmi-model.csv
breast-cancer,xu-sites.csv
cadmium,cadmium-sites.csv
copd,important-sites.csv
cotinine,discovery-sites.csv
crp,crp-sites.csv
dunedin-pace,pace-model.csv
dunedin-pace-poam38,dunedinpoam38-sites.csv
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
setwd("C:/Users/Ryan/Documents/GitHub/dnam-lung-cancer-screening-panel")
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
!any(is.na(idx))

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
    liftover(illumina.sites,chain.file =  "C://Users/Ryan/Documents/hg19ToHg38.chain"),
    illumina.sites[,c("source","details")])

hg19.regions.hg38 <- cbind(
    liftover(hg19.regions, chain.file = "C://Users/Ryan/Documents/hg19ToHg38.chain"),
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

fwrite(panel, file="panel-rjl.csv")

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

fwrite(panel.reduced, file="panel-reduced-rjl.csv")

## check the panels

nrow(panel)
## 3934

nrow(panel.reduced)
## 3409

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
##   1    2    3    4    5    6    7    8    9   10   15 
##3074  237   59   22    9    2    1    1    1    1    2

panel.reduced$source[which(sapply(mapping,length) > 10)]
#[1] "bmi/cadmium/cotinine/dunedin-pace-poam38/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"
#[2] "ancestry"                                                                                             
## (indicates that some CpG sites are being used for multiple ancestries)

table(sapply(strsplit(panel.reduced$source, "/"), length))
##    1    2    3    4    5    6    7    9 
## 3311   72   16    4    1    2    1    2 

panel.reduced[which(sapply(strsplit(panel.reduced$source, "/"),length) > 5),"source"]
## [1] "cadmium/cotinine/dunedin-pace/dunedin-pace-poam38/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"
## [2] "cadmium/cotinine/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"                                 
## [3] "bmi/cadmium/cotinine/dunedin-pace-poam38/educational-attainment/episcores/lung-cancer/smoking-cessation/smoking-status"         
## [4] "ancestry/cotinine/educational-attainment/lung-cancer/smoking-cessation/smoking-status"                                          
## [5] "cadmium/cotinine/episcores/lung-cancer/smoking-cessation/smoking-status"                                    

table(sapply(strsplit(panel.reduced$details, "/"), length))
##    1 
## 3409