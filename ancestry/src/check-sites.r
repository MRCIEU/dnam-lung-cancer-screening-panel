#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(vcfR)
library(ggfortify)

library(rmarkdown)
library(bookdown)
library(knitr)

sites.filename <- args[1]
genotypes.dir <- args[2]
input.filename <- args[3]
report.filename <- args[4]

pheno.url <- "http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel"

filenames <- list.files(genotypes.dir, "genotypes-.*.vcf.gz$", full.names=T)
vcf <- lapply(filenames,read.vcfR)

genotypes <- lapply(vcf, function(vcf) vcf@gt[,-1])
genotypes <- do.call(rbind, genotypes)

coords <- lapply(vcf, function(vcf) vcf@fix)
coords <- do.call(rbind, coords)
coords <- paste(coords[,"CHROM"], coords[,"POS"], sep="_")

genotypes[genotypes=="0|0"] <- 0
genotypes[genotypes=="0|1" | genotypes=="1|0"] <- 1
genotypes[genotypes=="1|1"] <- 2
genotypes <- apply(genotypes,2,as.integer)
rownames(genotypes) <- coords

pheno <- fread(pheno.url, select = 1:3)
pheno <- as.data.frame(pheno)
pheno <- pheno[pheno$sample %in% colnames(genotypes),]

sites <- fread(file=sites.filename)
sites$coords <- paste0("chr",sites$coords)
stopifnot(all(sites$coords %in% rownames(genotypes)))

common.samples <- intersect(colnames(genotypes), pheno$sample)

genotypes <- genotypes[match(sites$coords, rownames(genotypes)),
                       match(common.samples, colnames(genotypes))]

pheno <- pheno[match(common.samples, pheno$sample),]

knitr:::opts_chunk$set(
    fig.align="center",
    fig.dpi=320,
    fig.height=5,
    fig.width=5,
    message=FALSE,
    warning=FALSE,
    collapse=FALSE)

render(
    input=input.filename,
    output_file=basename(report.filename),
    output_dir=dirname(report.filename),
    output_format="html_document")

