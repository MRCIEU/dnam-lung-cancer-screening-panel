library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

filename <- geograbi.download.series.files(path=geo.dir, gses="GSE40279")
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
vars[,"age (y)"] <- as.numeric(vars[,"age (y)"])
meth <- geograbi.read.gse.matrix(filename)$data
meth <- as.matrix(meth)

vars$ancestry <- vars$ethnicity
