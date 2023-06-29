#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

sites.filename <- args[1]
datasets.filename <- args[2]
output.dir <- args[3]

library(data.table)

sites <- unique(fread(sites.filename)$details)
sites <- sites[grepl("^cg[0-9]+$",sites)]

datasets <- read.csv(datasets.filename)
geo.url <- "https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc="
datasets$dataset.url <- paste0(geo.url, datasets$gse)

source("src/render-ancestry-dnam-report-function.r")

ilogit2 <- function(x) 2^x / (1 + 2^x) ## needed by geograbi to load GEO datasets

for (i in 1:nrow(datasets)) {
    report.filename <- file.path(output.dir, paste0("ancestry-dnam-", datasets$gse[i], ".html"))
    if (file.exists(report.filename))
        next    
    cat(date(), "processing dataset", datasets$gse[i], "...\n")
    load.filename <- file.path("src", paste0("load-", datasets$gse[i], ".r"))
    stopifnot(file.exists(load.filename))
    script.env <- new.env()    
    source(load.filename,local=script.env)
    script.env <- as.list(script.env)
    render.ancestry.dnam.report(
        vars=script.env$vars,
        meth=script.env$meth,
        sites=sites,
        title=datasets$title[i],
        dataset.citation=datasets$dataset.citation[i],
        dataset.url=datasets$dataset.url[i],
        report.filename=report.filename)
}


