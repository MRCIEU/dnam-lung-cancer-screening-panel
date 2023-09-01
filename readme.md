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

From Twist:

> Using our default filters it looks like we would be able to cover 80.39% 
> of your targets, with the following regions filtered out: 
> [bed file](twist-clarification/Target_bases_not_covered_by_probes_Methyl_UniversityofBristol_lung-cancer-risk-panel_1X_MTE-93452736_hg19_230901090027.bed)
>
> We can recover a few more targets if we allow shifting of probes but the following would still remain uncovered: 
> [bed file](twist-clarification/all_target_segments_not_covered_by_probes_withshifting.bed)

We assessed the sites the could not be reliably targeted (without 'shifting'). 
Below we refer to these as 'uncovered'.

```
Rscript src/clarify.r
```

That script shows that

* all uncovered regions that are more than a single CpG site were included to estimate blood cell counts
* for each blood cell type, more than 50% of the regions specific to that cell type can be targeted entirely
* for each ancestry, more than 90% of the CpG sites specific to that ancestry can be targeted
* for each episcore, more than 97% of the CpG sites in the model can be targeted
* for each 'source' represented by the remaining uncovered CpG sites (these are mainly predictive models), 
  more than 93% of contributing CpG sites can be targeted

Here are the proportions of CpG sites that cannot be targeted for the last item (remaining sources):

|source              | proportion|
|:-------------------|-------:|
|bmi                 | 0.0030|
|alcohol-consumption | 0.0038|
|crp                 | 0.0088|
|smoking-cessation   | 0.0124|
|hdl                 | 0.0128|
|dunedin-pace        | 0.0195|
|breast-cancer       | 0.0400|
|dunedin-poam38      | 0.0541|
|cotinine            | 0.0625|


