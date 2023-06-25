#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggfortify)
library(geograbi) ## remotes::install_github("yousefi138/geograbi")

library(rmarkdown)
library(bookdown)
library(knitr)

sites.filename <- args[1] ## output/sites.csv
input.filename <-  args[2] ## src/check-sites-GSE40279.rmd
report.filename <- args[3]  ## output/check-sites-GSE40279.html

sites <- fread(sites.filename)

filename <- geograbi.download.series.files(gses="GSE40279")
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
meth <- geograbi.read.gse.matrix(filename)$data

knitr:::opts_chunk$set(
    fig.align="center",
    fig.dpi=320,
    fig.height=5,
    fig.width=5,
    message=FALSE,
    warning=FALSE,
    collapse=FALSE)

render(
    input=input.filename,
    output_file=basename(report.filename),
    output_dir=dirname(report.filename),
    output_format="html_document")

