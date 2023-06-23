#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(vcfR)
library(ggfortify)

panel.filename <- args[1]
genotypes.dir <- args[2]
plot.filename <- args[3]

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

panel.sites <- fread(file=panel.filename)
stopifnot(all(panel.sites$coords %in% rownames(genotypes)))

common.samples <- intersect(colnames(genotypes), pheno$sample)

genotypes <- genotypes[match(panel.sites$coords, rownames(genotypes)),
                       match(common.samples, colnames(genotypes))]

pheno <- pheno[match(common.samples, pheno$sample),]


pdf(plot.filename)
pca.ret <- prcomp(t(genotypes), scale=T)
autoplot(pca.ret, x=1, y=2, data=pheno, colour = 'super_pop')
autoplot(pca.ret, x=2, y=3, data=pheno, colour = 'super_pop')
autoplot(pca.ret, x=1, y=3, data=pheno, colour = 'super_pop')

for (super in unique(pheno$super_pop)) {
    insuper <- pheno$super_pop==super
    isvar <- apply(genotypes[,insuper],1,var,na.rm=T)>0
    pca.ret <- prcomp(t(genotypes[isvar,insuper]), scale=T)
    print(autoplot(pca.ret, x=1, y=2, data=pheno[insuper,], colour = 'pop'))
    print(autoplot(pca.ret, x=2, y=3, data=pheno[insuper,], colour = 'pop'))
    print(autoplot(pca.ret, x=1, y=3, data=pheno[insuper,], colour = 'pop'))
}
dev.off()


