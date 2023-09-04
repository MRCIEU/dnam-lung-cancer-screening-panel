# Lung Cancer Screening Panel Pipeline

## Dependencies

All required R packages can be installed as follows:

```
Rscript src/install-r-packages.r
```

## Ancestry

Before attempting to compile the final panel, verify that the
genomic regions for handling ancestry have been selected, 
i.e. [ancestry/output/sites.csv](ancestry/output/sites.csv) has been created.

The instructions to create this file can be found
in [ancestry/readme.md](ancestry/readme.md).

## Adding additional sites to the panel

1. Add a folder containing documentation and files relevant to the addition.  In particular, include one or more csv files providing either genomic regions (chr,start,end columns) or Illumina Beadchip CpG site identifiers (cpg column).

2. Add the csv file(s) to the appropriate list in `src/compile-panel.csv`.

## Compiling the final panel

```
Rscript src/compile-panel.r panel.csv panel-reduced.csv 
```

* `panel.csv` will provide a list of regions with some information about why they were selected (see 'source' and 'details').

* `panel-reduced.csv` will be the same as `panel.csv` but will have merged any overlapping regions.

Genomic coordinates will refer to the **hg19** human genome assembly.

## Checking the final panel

After creating the panel, run a few checks to reduce risk of errors.
See [checks/readme.md](checks/readm.md) for more details.

## Twist Clarification

Twist bioinformatics responded to our request 
noting that about 20% of our requested regions were difficult to 
target due to repetitive DNA.

Evaluation of these regions in terms of our planned use of them
suggested that we should try to target them but that, 
if targeting failed, then the impact would be small.
For example, most DNAm models would at most lose 2-3% of CpG sites, 
each ancestry would retain more than 90% of ancestry-specific sites, 
and each cell type would retain more than 50% of cell-type specific regions.

For more details, see 
[twist-clarification/readme.md](twist-clarification/readme.md).

The targeting statistics can be found here: 
[twist-clarification/output/stats.md](twist-clarification/output/stats.md).
