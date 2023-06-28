library(ggfortify)
library(rmarkdown)
library(bookdown)
library(knitr)
library(kableExtra)
library(tableone)
library(matrixStats)

source("src/lmable.r")

render.ancestry.dnam.report <- function(
    vars,
    meth,
    sites,
    title,
    dataset.citation,
    dataset.url,
    report.filename
    ) {

    stopifnot(is.matrix(meth))
    stopifnot(is.data.frame(vars))
    stopifnot(nrow(vars) == ncol(meth))
    stopifnot("ancestry" %in% colnames(vars))
    stopifnot(any(sites %in% rownames(meth)))
    
    knitr:::opts_chunk$set(
        fig.align="center",
        fig.dpi=320,
        fig.height=5,
        fig.width=5,
        echo=FALSE,
        message=FALSE,
        warning=FALSE,
        collapse=FALSE)   

    params <- list(
        report.title=title,
        dataset.citation=dataset.citation,
        dataset.url=dataset.url)
    
    render(
        input="src/ancestry-dnam-report.rmd",
        output_file=basename(report.filename),
        output_dir=dirname(report.filename),
        output_format="html_document",
        params=params)
}
