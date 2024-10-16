library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

filename <- geograbi.download.series.files(path=geo.dir, gses="GSE64940")
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
meth <- geograbi.read.gse.matrix(filename)$data
meth <- as.matrix(meth)

vars$ancestry <- vars$race

