lmable <- function(x, ...) {
    x <- as.data.frame(x,check.names=F)
    x[,ncol(x)] <- format(x[,ncol(x)], scientific=T,digits=3)
    kable(x,...)
}
