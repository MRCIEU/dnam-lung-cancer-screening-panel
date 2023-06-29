#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

stopifnot(length(args)==2)
stopifnot(file.exists(args[1]))

x <- read.csv(args[1])
x <- x[x$Predictor == "Elnet",]

stopifnot(nrow(x) == 1468)

x <- x[abs(x$Beta) > quantile(abs(x$Beta),probs=0.9),]

stopifnot(nrow(x) == 147)

colnames(x) <- tolower(colnames(x))

stopifnot("cpg" %in% colnames(x))

write.csv(x, file=args[2],quote=F,row.names=F)
