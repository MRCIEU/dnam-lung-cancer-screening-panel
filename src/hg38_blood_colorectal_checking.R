#rm(list = ls())
library(data.table)
#BiocManager::install("MungeSumstats")
library(MungeSumstats)

## Read in the panel
panel <- fread("panel-20230607.csv")

## Take blood and colorectal cancer (SEPT9) positions
blood <- panel[grep("Blood", panel$details),]
colo <- panel[grep("SEPT9", panel$details),]

## MungeSumstats expects certain column names...
colnames(blood) <- c("CHR", "BP", "end", "source", "SNP")

## Create id from "SNP" 
blood$SNP <- paste0(blood$SNP, 1:nrow(blood))

## Loop through and compare to Loyfer supplementary table S15
loc_blood <- MungeSumstats::liftover(sumstats_dt = blood,
                               convert_ref_genome = "GRCh37",
                               ref_genome = "GRCh38",
                               chrom_col = "CHR", start_col = "BP", 
                               end_col = "end")
  
loyfer <- fread("loyfer_supp_table_15.csv", skip = 2)
loyfer <- loyfer[grep("Blood", loyfer$Type),]

loc_blood$position <- paste0("chr", loc_blood$CHR, ":", loc_blood$BP, "-", loc_blood$end)

table(loc_blood$position %in% loyfer$position)
# FALSE  TRUE 
# 56    81 

checking <- loc_blood[!loc_blood$position %in% loyfer$position,]
## Appears to be where MungeSumstats liftover function rounded the start/end

## MungeSumstats expects certain column names...
colnames(colo) <- c("CHR", "BP", "end", "source", "SNP")
loc_colo <- MungeSumstats::liftover(sumstats_dt = colo,
                               convert_ref_genome = "GRCh37",
                               ref_genome = "GRCh38",
                               chrom_col = "CHR", start_col = "BP", 
                               end_col = "end")

## Manually check against wasserkort supplement (colorectal-cancer/wasserkort_supplement.pdf)

## SEPT9_AMP4 17 75368814 75369020
## SEPT9_AMP5 17 75369420 75369648
## SEPT9_AMP6 17 75370258 75370549

# These are correct