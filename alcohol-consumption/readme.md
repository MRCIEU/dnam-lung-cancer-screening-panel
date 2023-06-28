# Alcohol consumption

There are multiple EWAS/prediction models which exist as biomarkers of alcohol consumption. We have opted for:

* 144 CpGs from the dnamalci model developed from Liu et al. (10.1038/mp.2016.192) by Yousefi et al. (10.1101/591404). CpG names and coefficients were recovered via the `dnamalci.get.model("dnamalc.144cpg")` command in the dnamalci R package (https://github.com/yousefi138/dnamalci)
* 3 CpGs identified as predictive of alcohol consumption by Chamberlain et al. 2022 (10.1186/s13148-022-01376-7), mentioned in the Methods section of their text (cg06690548, cg03497652, and cg00716257)
* 14 CpGs which intersect between the respective top 100 sites from Lohoff et al. (10.1038/s41380-021-01378-6) and Dugue et al. (10.1111/adb.12855), as seen in `alc-ewas-intersection.R` [here](src/alc-ewas-intersection.R)
