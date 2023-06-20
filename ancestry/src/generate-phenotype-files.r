#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

dir.create(output.dir <- args[1])

pheno.url <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

library(data.table)
pheno <- fread(pheno.url, select = 1:3)
pheno <- as.data.frame(pheno)

for (varname in c("pop","super_pop")) {
    for (ancestry in unique(pheno[[varname]])) {
        design <- pheno[,c("sample",varname)]
        colnames(design) <- c("IID","phenotype")
        design$phenotype <- sign(design$phenotype == ancestry)
        write.table(
            x=design, 
            file=file.path(output.dir, paste0(ancestry, ".txt")),
            quote = F, row.names = F, col.names = T)
    }
}

