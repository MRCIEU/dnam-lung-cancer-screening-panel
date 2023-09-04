## Impact of repetitive DNA on proposed panel

### All affected regions of length at least 2 are for estimating blood cell counts


```r
stopifnot(all(!(panel$uncovered & panel$end-panel$start > 1) | panel$source == "blood-cell-types"))
```

### For each blood cell type, less than half of regions affected


```r
kable(with(panel[panel$source == "blood-cell-types",], stats.table(uncovered, details)))
```



|                 | affected| total| pct|
|:----------------|--------:|-----:|---:|
|Blood-NK         |        3|    25|  12|
|Blood-T          |        6|    25|  24|
|Blood-Mono+Macro |        9|    25|  36|
|Blood-B          |       10|    25|  40|
|Blood-Granul     |       10|    25|  40|

### For each ancestry, less than one-tenth of CpG sites affected


```r
kable(with(ancestry,stats.table(uncovered, ancestry)))
```



|    | affected| total| pct|
|:---|--------:|-----:|---:|
|ACB |        1|    50|   2|
|AMR |        1|    50|   2|
|CDX |        1|    50|   2|
|CHB |        1|    50|   2|
|CHS |        1|    50|   2|
|CLM |        1|    50|   2|
|EAS |        1|    50|   2|
|ESN |        1|    50|   2|
|GWD |        1|    50|   2|
|ITU |        1|    50|   2|
|JPT |        1|    50|   2|
|MSL |        1|    50|   2|
|SAS |        1|    50|   2|
|YRI |        1|    50|   2|
|ASW |        2|    50|   4|
|FIN |        2|    50|   4|
|IBS |        2|    50|   4|
|KHV |        2|    50|   4|
|LWK |        2|    50|   4|
|TSI |        2|    50|   4|
|PUR |        3|    50|   6|
|PEL |        4|    50|   8|

### For each episcore, less than 3% of CpG sites affected


```r
kable(with(episcores,stats.table(uncovered,gene)))
```



|        | affected| total| pct|
|:-------|--------:|-----:|---:|
|IGFBP1  |        1|   120| 0.8|
|TGFA    |        1|   100| 1.0|
|VEGFA   |        2|   159| 1.3|
|HGF     |        1|    72| 1.4|
|OSM     |        2|   125| 1.6|
|CXCL9   |        1|    56| 1.8|
|MMP12   |        5|   257| 1.9|
|S100A12 |        2|    82| 2.4|

### For most remaining sources, less than 2% of CpG sites affected


```r
kable(with(panel[! panel$source %in% c("blood-cell-types","ancestry","episcores"),],
           stats.table(uncovered,source)))
```



|                    | affected| total| pct|
|:-------------------|--------:|-----:|---:|
|alcohol-consumption |        1|   314| 0.3|
|bmi                 |        1|   397| 0.3|
|crp                 |        1|   155| 0.6|
|smoking-cessation   |        2|   222| 0.9|
|hdl                 |        1|    89| 1.1|
|dunedin-pace        |        3|   173| 1.7|
|cotinine            |        2|    55| 3.6|
|breast-cancer       |        4|   100| 4.0|
|dunedin-poam38      |        2|    46| 4.3|
