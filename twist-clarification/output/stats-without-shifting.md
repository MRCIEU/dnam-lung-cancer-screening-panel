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
|Blood-NK         |        4|    25|  16|
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
|AFR |        1|    50|   2|
|AMR |        1|    50|   2|
|BEB |        1|    50|   2|
|EUR |        1|    50|   2|
|GBR |        1|    50|   2|
|GIH |        1|    50|   2|
|MXL |        1|    50|   2|
|PJL |        1|    50|   2|
|CDX |        2|    50|   4|
|EAS |        2|    50|   4|
|IBS |        2|    50|   4|
|ITU |        2|    50|   4|
|JPT |        2|    50|   4|
|SAS |        2|    50|   4|
|STU |        2|    50|   4|
|YRI |        2|    50|   4|
|CLM |        3|    50|   6|
|ESN |        3|    50|   6|
|FIN |        3|    50|   6|
|KHV |        3|    50|   6|
|LWK |        3|    50|   6|
|MSL |        3|    50|   6|
|PUR |        3|    50|   6|
|CHS |        4|    50|   8|
|GWD |        4|    50|   8|
|TSI |        4|    50|   8|
|ACB |        5|    50|  10|
|ASW |        5|    50|  10|
|PEL |        5|    50|  10|
|CHB |        6|    50|  12|
|CEU |        7|    50|  14|

### For each episcore, less than 3% of CpG sites affected


```r
kable(with(episcores,stats.table(uncovered,gene)))
```



|        | affected| total| pct|
|:-------|--------:|-----:|---:|
|CXCL9   |        1|    56| 1.8|
|TGFA    |        2|   100| 2.0|
|IGFBP1  |        3|   120| 2.5|
|HGF     |        2|    72| 2.8|
|VEGFA   |        5|   159| 3.1|
|MMP12   |        9|   257| 3.5|
|OSM     |        7|   125| 5.6|
|S100A12 |        5|    82| 6.1|

### For most remaining sources, less than 2% of CpG sites affected


```r
kable(with(panel[! panel$source %in% c("blood-cell-types","ancestry","episcores"),],
           stats.table(uncovered,source)))
```



|                    | affected| total|  pct|
|:-------------------|--------:|-----:|----:|
|alcohol-consumption |        4|   314|  1.3|
|smoking-cessation   |        4|   222|  1.8|
|crp                 |        3|   155|  1.9|
|hdl                 |        2|    89|  2.2|
|bmi                 |        9|   397|  2.3|
|age                 |        3|    71|  4.2|
|dunedin-pace        |        8|   173|  4.6|
|dunedin-poam38      |        3|    46|  6.5|
|cotinine            |        4|    55|  7.3|
|prostate-cancer     |        2|    25|  8.0|
|egfr                |        2|    16| 12.5|
|breast-cancer       |       14|   100| 14.0|
