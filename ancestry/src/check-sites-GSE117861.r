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
input.filename <-  args[2]## input.filename="src/check-sites-GSE117861.rmd"
report.filename <- args[3]## report.filename="output/check-sites-GSE117861.html"

sites <- fread(sites.filename)

dir.create(geo.dir <- "geo", showWarnings=F)

gses <- c("GSE117859","GSE117860")

filenames <- geograbi.download.series.files(
    path=geo.dir,
    gses=gses)
samples <- lapply(filenames, geograbi.get.samples)
ids <- lapply(samples, function(samples) sub("_.*", "", samples$title))
vars <- lapply(samples, geograbi.extract.characteristics)
for (i in 1:length(vars))
    rownames(vars[[i]]) <- ids[[i]]
cols <- c("age","Sex","race","hiv","smoking","artadherence")
stopifnot(all(sapply(vars, function(vars) all(cols %in% colnames(vars)))))
vars <- lapply(vars, function(vars) vars[,cols])
for (i in 1:length(vars))
    vars[[i]]$age <- as.numeric(vars[[i]]$age)

meth <- lapply(gses, function(gse) {
    filename <- geograbi.download.supplementary.file(
        path=geo.dir,
        gse=gse,
        filename=paste0(gse, "_MatrixProcessed.txt.gz"))
    meth <- as.data.frame(fread(filename))
    rownames(meth) <- meth[,"ID_REF"]
    meth[,grepl("^Sample",colnames(meth))]
})

common.sites <- intersect(rownames(meth[[1]]), rownames(meth[[2]]))
meth <- lapply(meth, function(meth) meth[common.sites,])

impute.mean <- function(x) {
    idx <- which(is.na(x),arr.ind=T)
    if (length(idx) == 0) return(x)
    m <- rowMeans(x,na.rm=T)
    x[idx] <- m[idx[,2]]
    x
}

for (i in 1:length(meth))
    meth[[i]] <- impute.mean(meth[[i]])

stopifnot(all(sapply(1:length(vars), function(i) {
    identical(rownames(vars[[i]]), colnames(meth[[i]]))
})))

vars[[1]]$chip <- "450k"
vars[[2]]$chip <- "epic"

vars <- do.call(rbind, vars)
meth <- do.call(cbind, meth)

rownames(vars) <- paste(rownames(vars), vars$chip, sep="_")
colnames(meth) <- rownames(vars)

## adjust data for 450k/epic difference
diff <- rowMeans(meth[,vars$chip=="450k"],na.rm=T)-rowMeans(meth[,vars$chip=="epic"],na.rm=T)
meth[,vars$chip=="epic"] <- meth[,vars$chip=="epic"] + diff

## remove all missing values
meth <- impute.mean(meth)


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

