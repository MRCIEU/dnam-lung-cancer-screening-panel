library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

geo.dir <- "."

gses <- c("GSE72773","GSE72775")

filenames <- geograbi.download.series.files(path=geo.dir, gses=gses)
samples <- lapply(filenames, geograbi.get.samples)
vars <- lapply(samples, geograbi.extract.characteristics)

names(samples) <- names(vars) <- gses

for (i in 1:length(vars)) {
    vars[[i]] <- vars[[i]][,c("Sex","age","ethnicity")]
    vars[[i]]$batch <- gses[i]
}

meth <- lapply(gses, function(gse) {
    filename <- geograbi.download.supplementary.file(
        path=geo.dir,
        gse=gse,
        filename=paste0(gse, "_datBetaNormalized.csv.gz"))
    meth <- as.data.frame(fread(filename))
    rownames(meth) <- meth[,1]
    meth <- meth[,-1]
    for (i in 1:ncol(meth))
        if (!is.numeric(meth[[i]]))
            meth[[i]] <- as.numeric(meth[[i]])
    meth <- as.matrix(meth)
    id <- sub(".* (Sample[0-9]+)$", "\\1", samples[[gse]]$title)
    stopifnot(all(id %in% colnames(meth)))
    meth <- meth[,id]
    ## zero to normalize the datasets together
    meth - rowMeans(meth,na.rm=T)
})

stopifnot(identical(rownames(meth[[1]]), rownames(meth[[2]])))

meth <- do.call(cbind, meth)
vars <- do.call(rbind, vars)

vars$ancestry <- vars$ethnicity

