# Lung Cancer Screening Panel Pipeline


## Ancestry

Before attempting to compile the final panel, verify that the
genomic regions for handling ancestry have been selected, 
i.e. [ancestry/output/panel-sites/panel-sites.csv](ancestry/output/panel-sites/panel-sites.csv) has been created.

The instructions to create this file can be found
in [ancestry/readme.md](ancestry/readme.md).

## Compiling the final panel

```
Rscript src/compile-panel.r ./panel.csv annotation/
```


```



```

