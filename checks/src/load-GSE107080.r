library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

gse <- "GSE107080"

filename <- geograbi.download.series.files(path=geo.dir, gses=gse)
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
vars <- vars[,c("age","Sex","race","hiv","idu","hcv_dx","smoking","artadherence")]

vars[,"age"] <- as.numeric(vars[,"age"])

filename <- geograbi.download.supplementary.file(
    path=geo.dir,
    gse=gse,
    filename=paste0(gse, "_MatrixProcessed.txt.gz"))

meth <- as.data.frame(fread(filename))
rownames(meth) <- meth[,"ID_REF"]
meth <- meth[,!grepl("^(Detection|ID_REF)",colnames(meth))]
meth <- as.matrix(meth)

stopifnot(all(samples$description.2 %in% colnames(meth)))

meth <- meth[,samples$description.2]

vars$ancestry <- vars$race
