library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

gses <- c("GSE117859","GSE117860")

geo.dir <- "."

filenames <- geograbi.download.series.files(
    path=geo.dir,
    gses=gses)
samples <- lapply(filenames, geograbi.get.samples)
ids <- lapply(samples, function(samples) sub("_.*", "", samples$title))
vars <- lapply(samples, geograbi.extract.characteristics)
for (i in 1:length(vars))
    rownames(vars[[i]]) <- ids[[i]]
cols <- c("age","Sex","race","hiv","smoking","artadherence")
stopifnot(all(sapply(vars, function(vars) all(cols %in% colnames(vars)))))
vars <- lapply(vars, function(vars) vars[,cols])
for (i in 1:length(vars))
    vars[[i]]$age <- as.numeric(vars[[i]]$age)

meth <- lapply(gses, function(gse) {
    filename <- geograbi.download.supplementary.file(
        path=geo.dir,
        gse=gse,
        filename=paste0(gse, "_MatrixProcessed.txt.gz"))
    meth <- as.data.frame(fread(filename))
    rownames(meth) <- meth[,"ID_REF"]
    meth <- meth[,grepl("^Sample",colnames(meth))]
    meth <- as.matrix(meth)
    meth - rowMeans(meth,na.rm=T) ## harmonize datasets by mean-cent 
})

common.sites <- intersect(rownames(meth[[1]]), rownames(meth[[2]]))
meth <- lapply(meth, function(meth) meth[common.sites,])

stopifnot(all(sapply(1:length(vars), function(i) {
    identical(rownames(vars[[i]]), colnames(meth[[i]]))
})))

vars[[1]]$chip <- "450k"
vars[[2]]$chip <- "epic"

vars <- do.call(rbind, vars)
meth <- do.call(cbind, meth)

rownames(vars) <- paste(rownames(vars), vars$chip, sep="_")
colnames(meth) <- rownames(vars)

vars$ancestry <- vars$race






