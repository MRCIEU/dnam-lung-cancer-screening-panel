#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
library(ggfortify)
library(geograbi) ## remotes::install_github("yousefi138/geograbi")

sites.filename <- args[1] ## sites.csv
plot.filename <- args[2]  ## pca-of-dnam-GSE40279.pdf

sites <- fread(sites.filename)

filename <- geograbi.download.series.files(gses="GSE40279")
samples <- geograbi.get.samples(filename)
vars <- geograbi.extract.characteristics(samples)
meth <- geograbi.read.gse.matrix(filename)$data

all.sites <- sites$cpg
all.ret <- prcomp(t(meth[rownames(meth)%in% all.sites,]), scale=T)
t(apply(all.ret$x[,1:10], 2, function(p) coef(summary(lm(p~vars$ethnicity)))[2,]))
##         Estimate Std. Error     t value     Pr(>|t|)
## PC1   8.77697255  0.9404041   9.3331929 1.572673e-19
## PC2  -1.34525245  0.8168217  -1.6469353 1.000517e-01
## PC3  -9.00974440  0.4580811 -19.6684481 5.201032e-68
## PC4   5.78106279  0.4307882  13.4197323 1.913859e-36
## PC5   0.37315982  0.4086015   0.9132610 3.614418e-01
## PC6  -2.16861200  0.3737944  -5.8016168 1.023418e-08
## PC7  -0.05008514  0.3547369  -0.1411895 8.877637e-01
## PC8   1.28047108  0.3120230   4.1037712 4.578561e-05
## PC9   0.24839456  0.2901048   0.8562237 3.921877e-01
## PC10 -1.58446623  0.2753130  -5.7551442 1.330921e-08

sel.sites <- sites$cpg[sites$ancestry %in% c("AMR","CLM","MXL","PEL","PUR","EUR","CEU","FIN","GBR","IBS","TSI")]
sel.ret <- prcomp(t(meth[rownames(meth)%in% sel.sites,]), scale=T)
t(apply(sel.ret$x[,1:10], 2, function(p) coef(summary(lm(p~vars$ethnicity)))[2,]))
##         Estimate Std. Error      t value     Pr(>|t|)
## PC1  -5.30066809  0.5689955  -9.31583406 1.817113e-19
## PC2   1.53984369  0.5056855   3.04506196 2.419812e-03
## PC3  -6.33067792  0.2881565 -21.96958028 1.518814e-80
## PC4  -3.09193048  0.2880232 -10.73500545 7.147714e-25
## PC5   0.21089195  0.2533413   0.83244208 4.054633e-01
## PC6  -1.02022917  0.2404116  -4.24367677 2.516849e-05
## PC7   0.79634871  0.2294745   3.47031390 5.540253e-04
## PC8   0.60275353  0.1984072   3.03796203 2.476733e-03
## PC9  -1.02982128  0.1795564  -5.73536253 1.487589e-08
## PC10  0.01505591  0.1820531   0.08270063 9.341149e-01

pdf(plot.filename)
autoplot(all.ret, x=3, y=4, data=vars, colour = 'ethnicity', title="all sites")
autoplot(sel.ret, x=3, y=4, data=vars, colour = 'ethnicity', title="selected sites")
dev.off()
