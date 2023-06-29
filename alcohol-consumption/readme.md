# Alcohol consumption

There are multiple EWAS/prediction models which exist as biomarkers of alcohol consumption. We have opted for:

* `dnamalci-sites.csv` 144 CpGs from the dnamalci model developed from Liu et al. (10.1038/mp.2016.192) by Yousefi et al. (10.1101/591404). CpG names and coefficients were recovered via the `dnamalci.get.model("dnamalc.144cpg")` command in the dnamalci R package (https://github.com/yousefi138/dnamalci)

* `chamberlain-sites.csv` 3 CpGs identified as predictive of alcohol consumption by Chamberlain et al. 2022 (10.1186/s13148-022-01376-7), mentioned in the Methods section of their text (cg06690548, cg03497652, and cg00716257)

* `alc-ewas-sites.csv` 80 CpGs with p < 1e-10 in two large EWAS of alcohol exposure: Lohoff et al. (10.1038/s41380-021-01378-6) and Dugue et al. (10.1111/adb.12855).

```
Rscript src/alc-ewas-intersections.R
```

* `mccartney-sites.csv` 87 CpGs from McCartney et al.'s paper on prediction of complex traits and death (https://doi.org/10.1186/s13059-018-1514-1). The original 451 CpG predictor [mccartney-alc-predictor](mccartney_supplement_s3.csv) has been filtered by model coefficient > 0.4.

```
Rscript src/filter-mccartney-sites.R
```
