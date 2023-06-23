#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input.filename <- args[1]
output.filename <- args[2]

library(data.table)

panel.sites <- fread(input.filename)
fwrite(unique(panel.sites[,c("chr","pos")]),
       file=output.filename, sep="\t",col.names=F)
