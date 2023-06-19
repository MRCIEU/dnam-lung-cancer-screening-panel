rm(list = ls())

library(MungeSumstats)
library(data.table)

## Download mqtl information from: http://mqtldb.godmc.org.uk/downloads
mqtls <- fread("assoc_meta_all.csv.gz")
mqtls <- mqtls[mqtls$clumped == T,]
mqtls$CHR <- gsub(pattern = ":.*", "", mqtls$snp)
mqtls$CHR <- substr(x = mqtls$CHR, start = 4, stop = 4)
mqtls$BP <- gsub(pattern = "*:SNP", "", mqtls$snp)
mqtls$BP <- gsub(".*:", "", mqtls$BP)
mqtls$BP <- as.numeric(mqtls$BP)
mqtls <- mqtls[!is.na(mqtls$BP),]

save(mqtls, file = "godmc.Robj")

## Lookup table is chr_bp and rsids for GRCh37 from CIMBA (TNBC GWAS) - author emailed this to help with an MR
lookup <- fread("lookuptable.txt", header = F)
lookup$coords <- gsub("_[A-Z].*", "", lookup$V1)

## Keep relevant mQTL stats
mqtl_coords <- mqtls[,c("CHR", "BP", "beta_a1", "cpg")]
mqtl_coords$coords <- paste0(mqtl_coords$CHR, "_", mqtl_coords$BP)

## Merge lookup table and mQTL stats
mqtl_coords <- merge(mqtl_coords, lookup[,c("V2", "coords")], all.x = T, by = "coords")
mqtl_coords$coords <- NULL
colnames(mqtl_coords) <- c("CHR", "BP", "beta", "cpg" ,"SNP")

## Use MungeSumstats package function to convert GRCh37 to 38 so we have both assembly builds
liftover <- MungeSumstats::liftover(sumstats_dt = mqtl_coords, convert_ref_genome = "GRCh38", ref_genome = "GRCh37", chrom_col = "CHR", start_col = "BP", end_col = "BP", as_granges = F, verbose = T)
liftover$coords <- paste0(liftover$CHR, "_", liftover$BP)
rm(lookup)

save(liftover, file = "godmc_hg38.Robj")