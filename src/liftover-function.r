## BiocManager::install("liftOver")
library(liftOver)

library(rtracklayer)
path = system.file(package="liftOver", "extdata", "hg38ToHg19.over.chain")
ch = import.chain(path)
ch