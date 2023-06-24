#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

input.filename <- args[1]
output.filename <- args[2]

library(data.table)

sites <- fread(input.filename)
sites$chr <- paste0("chr",sites$chr)
fwrite(unique(sites[,c("chr","pos")]),
       file=output.filename, sep="\t",col.names=F)
