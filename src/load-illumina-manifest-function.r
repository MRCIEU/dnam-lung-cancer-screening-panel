load.illumina.manifest <- function(path) {
    url <- "https://webdata.illumina.com/downloads/productfiles/humanmethylation450/humanmethylation450_15017482_v1-2.csv"
    filename <- file.path(path, basename(url))
    filename.gz <- paste0(filename, ".gz")
    if (!file.exists(filename.gz)) {
        download.file(url, destfile=filename)
        system(paste("gzip", filename))
    }
    manifest.450k <- fread(filename.gz, skip=7,fill=TRUE)
    manifest.450k <- as.data.frame(manifest.450k)
    manifest.450k$chip <- "450k"

    url <- "ftp://webdata2:webdata2@ussd-ftp.illumina.com/downloads/productfiles/methylationEPIC/infinium-methylationepic-v-1-0-b5-manifest-file-csv.zip"
    filename.zip <- file.path(path, basename(url))
    if (!file.exists(filename.zip)) 
        download.file(url, destfile=filename.zip)
    filename <- unzip(filename.zip,list=T)$Name[1]
    manifest.epic <- fread(
        cmd=paste("unzip -p",filename.zip,filename),
        skip=7,fill=T)
    manifest.epic <- as.data.frame(manifest.epic)
    manifest.epic$chip <- "epic"

    cols <- intersect(
        colnames(manifest.450k),
        colnames(manifest.epic))

    rbind(manifest.450k[,cols], manifest.epic[,cols])
}
