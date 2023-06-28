# Checking the panel

Perform basic checks.
The script will generate an error if the reduced panel does not cover all panel regions.

* Input: Panel `../panel.csv` and non-overlapping panel `../panel-reduced.csv`
* Output: Report describing results of the check `outputs/perform-basic-checks.html`.
```
mkdir outputs
Rscript src/perform-basic-checks.r ../panel.csv ../panel-reduced.csv outputs/perform-basic-checks.html
```

Compare sites in the panel to DNAm models in [https://github.com/perishky/meffonym](meffonym).
The script will generate an error if selected models are not fully included in the panel.

* Input: Panel `../panel.csv`
* Output: Spreadsheet showing number of sites in each meffonym model and the number in the panel (`outputs/meffonym.csv`)
```
Rscript src/check-against-meffonym.r ../panel.csv outputs/meffonym.csv
```

Generate a randomly ordered spreadsheet for manually spot-checking panel regions.
```
Rscript src/make-check-panel.r ../panel-reduced.csv outputs/panel-reduced-py.csv
```

Compare DNA methylation clusters to ancestry

* Input: `../panel.csv`, `ancestry-dnam-datasets.csv`, `src/load-GSE*.r`, `src/ancestry-dnam-report.rmd`
* Output: `outputs/ancestry-dnam-GSE*.html`

```{r}
Rscript src/generate-ancestry-dnam-reports.r ../panel.csv ancestry-dnam-datasets.csv outputs
```

