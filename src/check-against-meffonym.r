#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

selected.filename <- args[1]
output.filename <- args[2]

library(meffonym) ## remotes::install_github("perishky/meffonym")

selection <- read.csv(selected.filename)

models <- meffonym.models()

results <- t(sapply(models, function(model) {
    sites <- names(meffonym.get.model(model)$coefs)
    c(model=model,n=length(sites), overlap=length(intersect(sites,selection$details)))
}))

included.models <- c("CXCL9","EN.RAGE","hannum","HGF","IGFBP-1","MMP-12","OSM","TGF.alpha","VEGFA")
stopifnot(
    results[results[,"model"] %in% included.models,"n"]
    ==       
    results[results[,"model"] %in% included.models,"overlap"])


write.csv(results, file=output.filename, row.names=F)

