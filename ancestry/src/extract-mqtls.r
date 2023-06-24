#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

output.filename <- args[1]

library(MungeSumstats) ## BiocManager::install("MungeSumstats")
library(data.table)

godmc.url <- "http://fileserve.mrcieu.ac.uk/mqtl/assoc_meta_all.csv.gz"
lookup.url <- "ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz"

## Download mqtl information from: http://mqtldb.godmc.org.uk/downloads
godmc.filename <- basename(godmc.url)
if (!file.exists(godmc.filename)) {
    tmp.filename <- file.path("/tmp", godmc.filename)
    download.file(godmc.url, destfile=tmp.filename)
    mqtls <- fread(tmp.filename)
    unlink(tmp.filename)
    mqtls <- mqtls[mqtls$clumped == T,]
    mqtls$CHR <- gsub(pattern = ":.*", "", mqtls$snp)
    mqtls$CHR <- sub("chr", "", mqtls$CHR)
    mqtls$BP <- gsub(pattern = "*:SNP", "", mqtls$snp)
    mqtls$BP <- gsub(".*:", "", mqtls$BP)
    mqtls$BP <- as.numeric(mqtls$BP)
    mqtls <- mqtls[!is.na(mqtls$BP),]
    nrow(mqtls) ## 264456
    mqtls$coords <- paste(mqtls$CHR,mqtls$BP,sep="_")
    mqtls$beta <- mqtls$beta_a1
    mqtls <- mqtls[,c("CHR","BP","beta","cpg","coords")]
    fwrite(mqtls, file=godmc.filename)
} else
    mqtls <- fread(godmc.filename)

## Download rsid lookup from the UCSC Genome Browser
lookup.filename <- basename(lookup.url)
if (!file.exists(lookup.filename)) {
    tmp.filename <- file.path("/tmp", lookup.filename)
    download.file(lookup.url, destfile=tmp.filename)
    lookup <- fread(tmp.filename,header=F)
    unlink(tmp.filename)
    cols <- c(
        "id","chr","start","end","rsid","score","strand",
        "ncbi","ucsc","observed","moltype","class","valid")
    lookup <- lookup[,1:length(cols)]
    colnames(lookup) <- cols
    nrow(lookup) ## 14831956
    lookup <- lookup[grepl("by-1000genomes",lookup$valid),]
    nrow(lookup) ## 14676828
    lookup$chr <- sub("chr", "", lookup$chr)    
    lookup$coords <- paste(lookup$chr,lookup$end,sep="_")
    mean(mqtls$coords %in% lookup$coords) ## 0.9613433
    lookup <- lookup[lookup$coords %in% mqtls$coords,]
    nrow(lookup) ## 209312
    fwrite(lookup,file=lookup.filename)
} else {
    lookup <- fread(lookup.filename)
}

## add rsid for each mqtl
mqtls$SNP <- lookup$rsid[match(mqtls$coords, lookup$coords)]
mqtls <- mqtls[!is.na(mqtls$SNP),]
nrow(mqtls) ## 254233

## convert to hg38 coordinates
mqtls.hg38 <- MungeSumstats::liftover(
    sumstats_dt=mqtls,
    convert_ref_genome="GRCh38",
    ref_genome="GRCh37",
    chrom_col="CHR",
    start_col="BP",
    end_col="BP",
    as_granges=F,
    verbose=T)
mqtls.hg38$coords <- with(mqtls.hg38, paste(CHR,BP,sep="_"))

nrow(mqtls.hg38) ## 253959

fwrite(mqtls.hg38, file=output.filename)

