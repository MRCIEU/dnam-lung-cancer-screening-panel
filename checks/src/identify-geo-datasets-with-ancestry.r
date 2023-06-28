library(geograbi) ## remotes::install_github("yousefi138/geograbi")
library(data.table)

## "GSE117859","GSE117860",

gses <- c("GSE40279","GSE64940","GSE77716", "GSE111629", "GSE107080", "GSE53740", "GSE72775", "GSE72773", "GSE54882", "GSE67705", "GSE122408", "GSE112596")

retrieve.samples <- function(gse) {
    cat("retrieving", gse, "...\n")
    geograbi.read.gse.matrix(gse, data=F)$samples
}

samples <- sapply(gses, retrieve.samples, simplify=F)

vars <- sapply(samples, geograbi.extract.characteristics, simplify=F)

sapply(vars, function(vars) {
    varname <- "race"
    if ("ethnicity" %in% colnames(vars))
        varname <- "ethnicity"
    freqs <- sort(table(vars[[varname]]),decreasing=T)
    n <- min(length(freqs),3)
    paste(names(freqs)[1:n],freqs[1:n],sep="=",collapse=";")
})

##                                          GSE40279 
## "Caucasian - European=426;Hispanic - Mexican=230" 
##                                          GSE64940 
##                           "AA=112;EA=91;Mixed=12" 
##                                          GSE77716 
##    "Mexican=276;Puerto Rican=220;Other Latino=61" 
##                                         GSE111629 
##                       "Caucasian=508;Hispanic=64" 
##                                         GSE107080 
##                     "White=329;Other=45;Black=31" 
##                                          GSE53740 
##               "White=311;Chinese=34;Other Race=7" 
##                                          GSE72775 
##                       "Caucasian=289;Hispanic=46" 
##                                          GSE72773 
##            "Caucasian=235;Hispanic=38;Tsimane=37" 
##                                          GSE54882 
##                                         "CEU=305" 
##                                          GSE67705 
##                                   "white=188;=96" 
##                                         GSE122408 
##                            "African American=180" 
##                                         GSE112596 
##                "Caucasian=78;Black=18;Hispanic=9" 
