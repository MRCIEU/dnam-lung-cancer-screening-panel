#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

panel.filename <- args[1]
ancestry.filename <- args[2]
episcores.filename <- args[3]
twist.bedfile <- args[4]
rmd.filename <- args[5]

if (F) {
    panel.filename <- "../panel.csv"
    ancestry.filename <- "../ancestry/output/sites.csv"
    episcores.filename <- "../episcores/episcore-sites.csv"
    twist.bedfile <- "all_target_segments_not_covered_by_probes_withshifting.bed"
    rmd.filename <- "output/stats.rmd"
}


library(knitr)
library(data.table)

## load panel
panel <- fread(panel.filename)
## convert to ucsc format including 0-based coords
panel$chr <- paste0("chr",panel$chr)
panel$start <- panel$start - 1

## load twist bed file
uncovered <- fread(twist.bedfile)
colnames(uncovered) <- c("chr","start","end")

## map twist bed file regions to panel
uncovered$idx <- sapply(1:nrow(uncovered), function(i) 
    with(uncovered, 
         which(
             panel$chr==uncovered$chr[i]
             & panel$start <= uncovered$start[i]
             & panel$end >= uncovered$end[i])))

## make sure all bed file regions map to panel regions
stopifnot(all(sapply(uncovered$idx,length) > 0))

## identify panel regions that cannot be 100% targeted
panel$uncovered <- 1:nrow(panel) %in% unlist(uncovered$idx)

ancestry <- fread(ancestry.filename)
episcores <- fread(episcores.filename)
ancestry$uncovered <- ancestry$cpg %in% panel$details[panel$uncovered]
episcores$uncovered <- episcores$cpg %in% panel$details[panel$uncovered]

stats.table <- function(status, group) {
    stats <- t(sapply(unique(group[status]), function(name) {
        affected <- status[group == name]
        c(affected=sum(affected),
          total=length(affected),
          pct=100*round(mean(affected),3))
    }))
    stats[order(stats[,"pct"]),]
}

knit(rmd.filename,output=sub("rmd$","md",rmd.filename))





