to.granges <- function(chr,start,end) {
    txt <- paste0(chr,":",start,"-",end)
    as(txt, "GRanges")
}

from.granges <- function(gr) {
    if ("CompressedGRangesList" %in% class(gr)) {
        len <- sapply(gr, length)
        if (any(len > 1)) {
            for (i in which(len > 1)) {
                chr <- as.character(seqnames(gr[[i]]))[1]
                start <- min(as.integer(start(gr[[i]])))
                end <- max(as.integer(end(gr[[i]])))
                gr[[i]] <- to.granges(chr,start,end)
            }
            gr <- unlist(gr)
        }
    }
    data.frame(
        chr=as.character(seqnames(gr)),
        start=as.integer(start(gr)),
        end=as.integer(end(gr)),
        stringsAsFactors=F)

}
