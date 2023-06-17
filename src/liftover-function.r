
retrieve.chain.file <- function(from,to,path=".") {
    capitalize <- function(str)
        paste0(toupper(substring(str,1,1)), tolower(substring(str,2)))
    from <- tolower(from)
    to <- capitalize(to)
    url <- paste0(
        "https://hgdownload.soe.ucsc.edu/goldenPath/",
        from,
        "/liftOver/",
        from,"To",to,".over.chain.gz")
    chain.file.gz <- file.path(path, basename(url))
    chain.file <- sub(".gz$", "", chain.file.gz)
    if (!file.exists(chain.file)) {
        download.file(
            url,
            destfile=chain.file.gz)
        system(paste("gunzip", chain.file.gz))
    }
    chain.file
}

liftover <- function(coords, chain.file) {
    stopifnot(all(c("chr","start","end") %in% colnames(coords)))    
    chain <- import.chain(chain.file)
    coords.from <- with(coords, to.granges(chr,start,end))
    coords.to <- liftOver(coords.from, chain)
    from.granges(coords.to)
}
