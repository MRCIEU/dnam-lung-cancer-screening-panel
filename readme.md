# DNA methylation Panel for Lung Cancer Screening 

The panel includes CpG sites known to be associated with factors relevant to lung cancer risk:

- demographics ([age](age), [ancestry](ancestry), [educational attainment](educational-attainment))
- exposures (smoking history--[status](smoking-status), [former](smoking-former), [cessation](smoking-cessation), [cotinine](cotinine)--[alcohol consumption](alcohol-consumption), [cadmium](cadmium), [lead](lead))
- aging ([Hannum](age) and Dunedin [PACE](dunedin-pace)/[PoAm](dunedin-poam38) clocks)
- [blood cell types](blood-cell-types)
- cancer risk ([lung](lung-cancer), [breast](breast-cancer), [colorectal](colorectal-cancer), [prostate](prostate))
- cardiovascular disease risk ([COPD](copd), [HDL cholesterol](hdl), [BMI](bmi))
- protein abundance ([OSM, EN-RAGE, CXCL9, VEGFA, TGFa, IGFBP1, MMP12, HGF](episcores), [IL6](il6), [EGFR](egfr), [CRP](crp))

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

Details can be found here:
[twist-clarification/readme.md](twist-clarification/readme.md).

Additional coverage statistics we generated can be found here: 
[twist-clarification/output/stats.md](twist-clarification/output/stats.md).

## Twist Panel

Quality control analyses and outputs
for the probes used in the final panel are described here: [twist-panel/readme.md](twist-panel/readme.md)

Outputs can be found here: 
[twist-panel/output/stats.md](twist-panel/output/stats.md).

