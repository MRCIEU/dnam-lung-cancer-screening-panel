# Blood cell-type specific DNA methylation

> Loyfer, N., … Kaplan, T. (2023). A DNA methylation atlas of normal human cell types. Nature, 613(7943), 355–364
> (PMID [https://pubmed.ncbi.nlm.nih.gov/36599988/](36599988)).

**Supplementary Dataset 1**

"Genome-wide set of unmethylated regions per cell type, annotated. Zip
file contains 39 bed files, each with all genomic regions (blocks of
at least four CpGs in which at least 85% of sequenced fragments are
unmethylated in at least 85% of covered CpGs)."

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC9811898/bin/41586_2022_5580_MOESM5_ESM.zip

- Blood-B.4CpGs.U85.f85.bed
- Blood-Granul.4CpGs.U85.f85.bed
- Blood-Mono+Macro.4CpGs.U85.f85.bed
- Blood-NK.4CpGs.U85.f85.bed
- Blood-T.4CpGs.U85.f85.bed

*Note that coordinates refer to the hg19 genome assembly.*

`regions-hg19.csv` was generated from the Loyfer et al. file as follows:

```r
loyfer.table <- read.csv("loyfer_supp_table_15.csv",skip=2)
loyfer.table <- loyfer.table[grep("Blood",loyfer.table$Type),]
loyfer.table <- loyfer.table[,c("Type","chr","start","end","position","Number.of.CpGs","Length","Gene")]
colnames(loyfer.table)[colnames(loyfer.table) == "Type"] <- "Details"
colnames(loyfer.table) <- tolower(colnames(loyfer.table))
loyfer.table$chr <- sub("chr","",loyfer.table$chr)
loyfer.table$position <- sub("chr","",loyfer.table$position)
write.csv(loyfer.table, file="regions-hg19.csv",row.names=F)
```
