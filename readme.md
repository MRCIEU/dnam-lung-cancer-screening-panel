# Lung Cancer Screening Panel Pipeline

## Ancestry

Before attempting to compile the final panel, verify that the
genomic regions for handling ancestry have been selected, 
i.e. [ancestry/output/sites.csv](ancestry/output/sites.csv) has been created.

The instructions to create this file can be found
in [ancestry/readme.md](ancestry/readme.md).

## Compiling the final panel

```
Rscript src/compile-panel.r panel.csv panel-reduced.csv annotation/
```

* `panel.csv` will provide a list of regions with some information about why they were selected (see 'source' and 'details').

* `panel-reduced.csv` will be the same as `panel.csv` but will have merged any overlapping regions.

* `annotation/` will contain any genome or microarray annotation files that were used in the computation.

Genomic coordinates will refer to the hg38 human genome assembly.

## Checking the final panel

```
Rscript src/check-against-meffonym.r panel-reduced.csv 
```

Note: We perhaps conspicuously don't include PhenoAge, CRP and DunedinPoam38.

## Adding additional sites to the panel

1. Add a folder containing documentation and files relevant to the addition.  In particular, include one or more csv files providing either genomic regions (chr,start,end columns) or Illumina Beadchip CpG site identifiers (cpg column).

2. Add the csv file(s) to the appropriate list in `src/compile-panel.csv`.
