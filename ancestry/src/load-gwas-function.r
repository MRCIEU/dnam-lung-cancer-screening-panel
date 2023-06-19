library(data.table)

load.gwas <- function(path, ancestry, model=c("glm","fst")) {
  files <- list.files(
      path,
      pattern=paste0(ancestry, "-.*", model),
      full.names=T)
  stopifnot(length(files)>0)
  dat <- do.call(rbind, lapply(files, fread))
  dat$coords <- paste(dat[["#CHROM"]], dat[["POS"]], sep="_")
  dat <- dat[order(dat[["#CHROM"]], as.integer(dat[["POS"]])),]
  dat
}

