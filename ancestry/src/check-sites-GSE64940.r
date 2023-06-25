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
input.filename <-  args[2] ## input.filename="src/check-sites-GSE64940.rmd"
report.filename <- args[3]  ## report.filename="output/check-sites-GSE64940.html"

sites <- fread(sites.filename)

dir.create(geo.dir <- "geo", showWarnings=F)

filename <- geograbi.download.series.files(path=geo.dir, gses="GSE64940")
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
ilogit2 <- function(x) 2^x / (1 + 2^x)
meth <- geograbi.read.gse.matrix(filename)$data

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

