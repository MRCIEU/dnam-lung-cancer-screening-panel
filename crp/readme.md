# CRP 

* `crp-sites.csv` A DNAm score of 8 CpGs was shown to explain ~5% variance in CRP levels in CHARGE, in a 2016 paper by Ligthart et al. (https://doi.org/10.1186/s13059-016-1119-5). The authors report these findings in the text of their results section

* `reduced-hillary-model.csv` A model of CRP was developed in Generation Scotland by Robert Hillary using Elastic Net. The model included 1468 CpG sites, however, evaluation of the model in ALSPAC showed identical performance using only the CpG sites with the largest coefficients (top 10%, i.e. 147 CpG sits).

```
Rscript src/reduce-hillary-model.r hillary-crp-model-20230629.csv reduced-hillary-model.csv
```
