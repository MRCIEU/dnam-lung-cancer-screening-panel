library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

gse <- "GSE111629"

filename <- geograbi.download.series.files(path=geo.dir, gses=gse)
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)

vars[,"age"] <- as.numeric(vars[,"age"])


filename <- geograbi.download.supplementary.file(
    path=geo.dir,
    gse=gse,
    filename=paste0(gse, "_PEGblood_450kMethylationDataBackgroundNormalized.txt.gz"))

meth <- as.data.frame(fread(filename))
rownames(meth) <- meth[,1]
meth <- meth[,-1]
meth <- as.matrix(meth)

stopifnot(all(samples$source_name_ch1 %in% colnames(meth)))

meth <- meth[,samples$source_name_ch1]

vars$ancestry <- vars$ethnicity

