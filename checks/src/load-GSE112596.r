library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

gse <- "GSE112596"

filename <- geograbi.download.series.files(path=geo.dir, gses=gse)
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
vars <- vars[,c("cell type","diagnosis","age (years)", "gender","race","therapy","disease duration (months)")]

filename <- geograbi.download.supplementary.file(
    path=geo.dir,
    gse=gse,
    filename=paste0(gse, "_GEOmatrixprocessed.csv.gz"))

meth <- as.data.frame(fread(filename))
rownames(meth) <- meth[,1]
meth <- meth[,-1]
meth <- meth[,!grepl("p<",meth[2,])]
meth <- as.matrix(meth)

stopifnot(all(samples$title %in% colnames(meth)))

meth <- meth[,samples$title]

vars$ancestry <- vars$race
