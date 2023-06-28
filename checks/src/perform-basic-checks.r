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

panel.filename <- args[1] 
panel.reduced.filename <- args[2] 
report.filename <- args[3]

panel <- fread(panel.filename)
panel.reduced <- fread(panel.reduced.filename)

knitr:::opts_chunk$set(
    fig.align="center",
    fig.dpi=320,
    fig.height=5,
    fig.width=5,
    message=FALSE,
    warning=FALSE,
    collapse=FALSE)

render(
    input="src/perform-basic-checks.rmd",
    output_file=basename(report.filename),
    output_dir=dirname(report.filename),
    output_format="html_document")



