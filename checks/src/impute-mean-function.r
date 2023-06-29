## for each missing value,
## imputes the mean from that row
impute.mean <- function(x) {
    idx <- which(is.na(x),arr.ind=T)
    if (length(idx) == 0) return(x)
    m <- rowMeans(x,na.rm=T)
    x[idx] <- m[idx[1,]]
    x
}
