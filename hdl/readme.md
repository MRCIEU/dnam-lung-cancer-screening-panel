# HDL

* 89 CpGs from McCartney et al.'s paper on prediction of complex traits and death (https://doi.org/10.1186/s13059-018-1514-1). The original 737 CpG predictor [mccartney-hdl-predictor](mccartney_supplement_s6.csv) has been filtered by effect size to retain the probes with coefficients larger than 0.2

```
Rscript src/filter-mccartney-sites.R
```
