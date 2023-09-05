#!/usr/bin/env Rscript
args <- commandArgs(trailingOnly=TRUE)

panel.filename <- args[1]
targets.filename <- args[2]
covered.filename <- args[3]
uncovered.filename <- args[4]
probes.filename <- args[5]
rmd.filename <- args[6]

if (F) {
    
    panel.filename <- "../panel.csv"
    
    targets.filename <- "data/all_target_segments_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed"
    
    covered.filename <- "data/all_target_segments_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed"
    
    uncovered.filename <- "data/all_target_segments_not_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed"
    
    probes.filename <- "data/merged_probe_file_shareable_Methyl_UniversityofBristol_lung-cancer-risk-panel_shift100_1X_MTE-93452736_hg19_LOWStringency.bed"

    rmd.filename <- "output/stats.rmd"
}


library(knitr)
library(data.table)

## load panel
panel <- fread(panel.filename)
## convert to ucsc format including 0-based coords
panel$chr <- paste0("chr",panel$chr)
#panel$start <- panel$start - 1

cols <- c("chr","start","end")
targets <- fread(targets.filename,col.names=c(cols,"id"))
covered <- fread(covered.filename,col.names=cols)
uncovered <- fread(uncovered.filename,col.names=cols)
probes <- fread(probes.filename,col.names=cols)

## map uncovered regions to panel
panel$uncovered <- sapply(1:nrow(panel), function(i) { 
    any(panel$chr[i]==uncovered$chr
        & panel$start[i] <= uncovered$start+1
        & panel$end[i] >= uncovered$end-1)
})

## map covered regions to panel
panel$covered <- sapply(1:nrow(panel), function(i) {
    any(covered$chr==panel$chr[i]
        & covered$start<=panel$start[i]+1
        & covered$end>=panel$end[i]-1)
})

panel$has.coverage <- sapply(1:nrow(panel), function(i) {
    any(probes$chr==panel$chr[i]
        & (probes$start >= panel$start[i]-1 & probes$start <= panel$end[i]+1
           | probes$end >= panel$start[i]-1 & probes$end <= panel$end[i]+1
           | panel$start[i] >= probes$start-1 & panel$start[i] <= probes$end+1
           | panel$end[i] >= probes$start-1 & panel$end[i] <= probes$end+1))
})


with(panel,table(uncovered,covered,has.coverage))
## , , has.coverage = FALSE

##          covered
## uncovered FALSE TRUE
##     FALSE     0    0
##     TRUE     20    0

## , , has.coverage = TRUE

##          covered
## uncovered FALSE TRUE
##     FALSE     0 4506
##     TRUE      9    0

## 20 regions have no coverage
## 9 regions have partial coverage

stats.table <- function(status, group) {
    stats <- t(sapply(unique(group[status]), function(name) {
        affected <- status[group == name]
        c(affected=sum(affected),
          total=length(affected),
          pct=100*round(mean(affected),3))
    }))
    stats[order(stats[,"pct"]),]
}

knit(rmd.filename,output=sub("rmd$","md",rmd.filename))
