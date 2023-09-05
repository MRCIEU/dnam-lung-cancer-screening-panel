## Proposed panel vs desired targets

### Logical checks


```r
stopifnot(!panel$uncovered | !panel$covered)
stopifnot(!panel$covered | panel$has.coverage)
```

### All affected regions of length at least 2 are for estimating blood cell counts


```r
stopifnot(all(!(panel$uncovered & panel$end-panel$start > 1) | panel$source == "blood-cell-types"))
```

### For each blood cell type, at most 20% of regions are affected


```r
kable(with(panel[panel$source == "blood-cell-types",], stats.table(uncovered, details)))
```



|                 | affected| total| pct|
|:----------------|--------:|-----:|---:|
|Blood-Granul     |        2|    25|   8|
|Blood-T          |        2|    25|   8|
|Blood-B          |        4|    25|  16|
|Blood-Mono+Macro |        5|    25|  20|

### For each blood cell type, at most 8% of the regions have no coverage at all


```r
kable(with(panel[panel$source == "blood-cell-types",], stats.table(!has.coverage, details)))
```



|                 | affected| total| pct|
|:----------------|--------:|-----:|---:|
|Blood-Granul     |        2|    25|   8|
|Blood-Mono+Macro |        2|    25|   8|

### For most remaining sources, less than 2% of CpG sites affected


```r
kable(with(panel[panel$source != "blood-cell-types",],
           stats.table(uncovered,source)))
```



|                    | affected| total| pct|
|:-------------------|--------:|-----:|---:|
|alcohol-consumption |        1|   314| 0.3|
|ancestry            |        6|  1550| 0.4|
|episcores           |        5|   971| 0.5|
|dunedin-pace        |        2|   173| 1.2|
|dunedin-poam38      |        2|    46| 4.3|


```r
kable(with(panel[panel$source != "blood-cell-types",],
           stats.table(!has.coverage,source)))
```



|                    | affected| total| pct|
|:-------------------|--------:|-----:|---:|
|alcohol-consumption |        1|   314| 0.3|
|ancestry            |        6|  1550| 0.4|
|episcores           |        5|   971| 0.5|
|dunedin-pace        |        2|   173| 1.2|
|dunedin-poam38      |        2|    46| 4.3|

