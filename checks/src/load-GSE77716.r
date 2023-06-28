library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

gse <- "GSE77716"

filename <- geograbi.download.series.files(path=geo.dir, gses=gse)
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
vars <- vars[,c("gender","sample type","ethnicity")]

filename <- geograbi.download.supplementary.file(
    path=geo.dir,
    gse=gse,
    filename=paste0(gse, "_Matrix_processed.tsv.gz"))

meth <- as.data.frame(fread(filename))
rownames(meth) <- meth[,"ID_REF"]
meth <- meth[,!grepl("^(Detection|ID_REF)",colnames(meth))]
meth <- as.matrix(meth)

stopifnot(all(samples$title %in% colnames(meth)))

meth <- meth[,samples$title]

vars$ancestry <- vars$ethnicity
