#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggfortify)
library(geograbi) ## remotes::install_github("yousefi138/geograbi")

library(rmarkdown)
library(bookdown)
library(knitr)
library(kableExtra)
library(tableone)

sites.filename <- args[1] ## sites.filename="output/sites.csv"
input.filename <-  args[2] ## input.filename="src/check-sites-GSE77716.rmd"
report.filename <- args[3]  ## report.filename="output/check-sites-GSE77716.html"

sites <- fread(sites.filename)

dir.create(geo.dir <- "geo", showWarnings=F)

gse <- "GSE77716"

filename <- geograbi.download.series.files(path=geo.dir, gses=gse)
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
vars <- vars[,c("gender","sample type","ethnicity")]

filename <- geograbi.download.supplementary.file(
    path=geo.dir,
    gse=gse,
    filename=paste0(gse, "_Matrix_processed.tsv.gz"))

meth <- as.data.frame(fread(filename))
rownames(meth) <- meth[,"ID_REF"]
meth <- meth[,!grepl("^(Detection|ID_REF)",colnames(meth))]

stopifnot(all(samples$title %in% colnames(meth)))

meth <- meth[,samples$title]

knitr:::opts_chunk$set(
    fig.align="center",
    fig.dpi=320,
    fig.height=5,
    fig.width=5,
    message=FALSE,
    warning=FALSE,
    collapse=FALSE)

source("src/lmable.r")

render(
    input=input.filename,
    output_file=basename(report.filename),
    output_dir=dirname(report.filename),
    output_format="html_document")

