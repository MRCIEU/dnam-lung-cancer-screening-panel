#!/usr/bin/env Rscript
args = commandArgs(trailingOnly=TRUE)

library(data.table)
args <-   c("output/pheno",
  "output/gwas-glm",
  "output/gwas-fst",
  "godmc-hg38.csv.gz",
  "output/panel-sites")
#' For each ancestry, identify the top 50 SNPs
#' in terms of meQTL effect such that
#' the GWAS Fst > 95th percentile of meQTLs
#' with logistic GWAS p < 5e-8.
#'
#' (about 5 minutes per ancestry)

pheno.dir <- args[1]
glm.dir <- args[2]
fst.dir <- args[3]
godmc.file <- args[4]
dir.create(output.dir <- args[5])

source("src/load-gwas-function.r")
## ------------------------------
## load data

ancestries <- list.files(pheno.dir, pattern="^[A-Z]+.txt$")
ancestries <- sub(".txt", "", ancestries)

## mQTL summary statistics (GoDMC)
godmc <- fread(godmc.file)
godmc.pct <- ecdf(abs(godmc$beta))
godmc$beta.pct <- godmc.pct(abs(godmc$beta))

for (ancestry in ancestries) {
    cat("----------------------------------------\n")
    cat(date(), ancestry, "\n")

    ## summary stats for GWAS (Fst)
    dat <- load.gwas(fst.dir, ancestry, "fst")
    
    ## summary stats for GWAS (logistic regression)
    dat.glm <- load.gwas(glm.dir, ancestry, "glm")
    
    ## add logistic GWAS summary statistics 
    stopifnot(identical(dat$coords, dat.glm$coords))
    dat$z.glm <- dat.glm$Z_STAT
    dat$p.glm <- dat.glm$P
        
    ## add mQTL statistics
    idx <- match(dat$coords, godmc$coords)
    dat$beta.godmc <- godmc$beta[idx]
    dat$cpg.godmc <- godmc$cpg[idx]
    dat$pct.godmc <- godmc$beta.pct[idx]

    ## Fst threshold at 95% percentile of Fst for mQTLs
    ## with logistic p < 5e-8
    threshold <- with(dat, quantile(HUDSON_FST[p.glm < 5e-8 & !is.na(beta.godmc)], probs=0.95, na.rm=T))
    
    top.pairs <- dat[which(dat$HUDSON_FST > threshold & !is.na(beta.godmc)),]
    top.pairs$is.top50 <- 1:nrow(top.pairs) %in% order(abs(top.pairs$beta.godmc),decreasing=T)[1:50]

    fwrite(top.pairs, file=file.path(output.dir, paste0(ancestry, ".csv")))
    
    ## output results

    ## correlation between logistic z-score and Fst
    cat("cor(z.glm,fst) =", with(dat, cor(abs(z.glm), HUDSON_FST, use="pairwise.complete")), "\n")
    
    ## Fst relative to logistic p-value = 5e-8
    cat("quantiles(fst[p.glm<5e-8]) =", with(dat, quantile(HUDSON_FST[p.glm < 5e-8], na.rm=T)), "\n")
    cat("quantiles(fst[p.glm>5e-8]) =", with(dat, quantile(HUDSON_FST[p.glm > 5e-8], na.rm=T)), "\n")
    
    ## number logistic associations at p < 5e-8
    cat("number p.glm < 5e-8 =", sum(dat$p.glm < 5e-8,na.rm=T), "\n")

    ## fst threshold
    cat("fst threshold =", threshold, "\n")

    ## GWAS associations with Fst > threshold and meQTL
    cat("number fst > threshold and meQTL =", nrow(top.pairs), "\n")
    
    ## distribution of mQTL effects among top 50 sites
    cat("quantiles(top 50 mQTL effects) =", quantile(with(top.pairs, abs(beta.godmc[is.top50]))), "\n")
    
    ## percentiles of mQTL effects among top 50 sites
    cat("... corresponding percentiles =", quantile(ecdf(abs(godmc$beta))(abs(top.pairs$beta.godmc))), "\n")
}


cols <- c("chr","pos","id","obs.ct","hudson.fst","coords","z.glm","p.glm","beta.godmc","cpg.godmc","pct.godmc","is.top50")

panel.sites <- lapply(ancestries, function(ancestry) {
    filename <- file.path(output.dir, paste0(ancestry, ".csv"))
    sites <- fread(filename)
    stopifnot(sum(sites$is.top50,na.rm=T)>0)
    sites <- sites[sites$is.top50,]
    stopifnot(length(cols) == ncol(sites))
    colnames(sites) <- cols
    sites$start <- sites$end <- sites$pos
    sites$cpg <- sites$cpg.godmc
    sites$ancestry <- ancestry
    sites
})
panel.sites <- do.call(rbind, panel.sites)

fwrite(panel.sites, file=file.path(output.dir, "panel-sites.csv"))



## ----------------------------------------
##     ACB
## cor(z.glm,fst) = 0.8057804
## quantiles(fst[p.glm<5e-8]) = -0.000911166 0.0610505 0.0923413 0.1414502 0.602116
## quantiles(fst[p.glm>5e-8]) = -0.00289138 -0.00106996 0.00392139 0.0172468 0.361634
## number p.glm < 5e-8 = 2085812
## fst threshold = 0.366137
## number fst > threshold and meQTL = 285
## quantiles(top 50 mQTL effects) = 0.3589414 0.4335772 0.5041716 0.6758411 1.532532
## ... corresponding percentiles = 0.0004279701 0.2563998 0.4813916 0.7020123 0.9984294
## Loading objects:
##       dat
## ----------------------------------------
##     AFR
## cor(z.glm,fst) = 0.9090676
## quantiles(fst[p.glm<5e-8]) = 0.000186468 0.0394018 0.0710868 0.132892 0.897636
## quantiles(fst[p.glm>5e-8]) = -0.000519405 0.000501096 0.00463758 0.012219 0.23092
## number p.glm < 5e-8 = 8923152
## fst threshold = 0.5135426
## number fst > threshold and meQTL = 1293
## quantiles(top 50 mQTL effects) = 0.6012073 0.633715 0.7083299 0.9071415 1.532532
## ... corresponding percentiles = 0.0002285277 0.2396716 0.4689098 0.6859447 0.9984294
## Loading objects:
##       dat
## ----------------------------------------
##     AMR
## cor(z.glm,fst) = 0.9106141
## quantiles(fst[p.glm<5e-8]) = 0.00342562 0.03572135 0.0498023 0.0707248 0.378841
## quantiles(fst[p.glm>5e-8]) = -0.000872405 -0.00039968 0.00157086 0.00758021 0.194131
## number p.glm < 5e-8 = 1886375
## fst threshold = 0.1732978
## number fst > threshold and meQTL = 452
## quantiles(top 50 mQTL effects) = 0.4083244 0.4356362 0.467019 0.5460069 1.127537
## ... corresponding percentiles = 5.401565e-05 0.1916517 0.4062974 0.6742815 0.9922591
## Loading objects:
##       dat
## ----------------------------------------
##     ASW
## cor(z.glm,fst) = 0.7737309
## quantiles(fst[p.glm<5e-8]) = 0.00725419 0.0805743 0.112407 0.157941 0.495862
## quantiles(fst[p.glm>5e-8]) = -0.00451006 -0.00226252 0.00322114 0.0183232 0.343529
## number p.glm < 5e-8 = 512837
## fst threshold = 0.3268863
## number fst > threshold and meQTL = 55
## quantiles(top 50 mQTL effects) = 0.0705155 0.1282891 0.1707521 0.258788 1.294406
## ... corresponding percentiles = 0.004362802 0.2579102 0.410415 0.6281147 0.9958699
## Loading objects:
##       dat
## ----------------------------------------
##     BEB
## cor(z.glm,fst) = 0.7509986
## quantiles(fst[p.glm<5e-8]) = 0.0238947 0.0554142 0.083516 0.113937 0.705645
## quantiles(fst[p.glm>5e-8]) = -0.00321903 -0.00166169 0.00185287 0.013371 0.7256
## number p.glm < 5e-8 = 182463
## fst threshold = 0.2679628
## number fst > threshold and meQTL = 42
## quantiles(top 50 mQTL effects) = 0.0442842 0.09442485 0.134654 0.2274855 0.48377
## ... corresponding percentiles = 0.005434805 0.1547995 0.3313465 0.5845418 0.8847763
## Loading objects:
##       dat
## ----------------------------------------
##     CDX
## cor(z.glm,fst) = 0.7987636
## quantiles(fst[p.glm<5e-8]) = 0.00424579 0.08025 0.12114 0.178094 0.850444
## quantiles(fst[p.glm>5e-8]) = -0.00298217 -0.00039968 0.00737566 0.0293675 0.545763
## number p.glm < 5e-8 = 1179326
## fst threshold = 0.4012553
## number fst > threshold and meQTL = 282
## quantiles(top 50 mQTL effects) = 0.325042 0.3754365 0.4235452 0.5252375 1.133741
## ... corresponding percentiles = 0.00420491 0.2065808 0.4100702 0.6313806 0.9924337
## Loading objects:
##       dat
## ----------------------------------------
##     CEU
## cor(z.glm,fst) = 0.8022591
## quantiles(fst[p.glm<5e-8]) = 0.00738666 0.0738138 0.103041 0.142041 0.723754
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.00102915 0.00424947 0.0203776 0.539383
## number p.glm < 5e-8 = 595692
## fst threshold = 0.295337
## number fst > threshold and meQTL = 191
## quantiles(top 50 mQTL effects) = 0.2808921 0.3368212 0.4206102 0.489851 0.8987742
## ... corresponding percentiles = 0.0003947297 0.2577793 0.4541303 0.6810521 0.9814186
## Loading objects:
##       dat
## ----------------------------------------
##     CHB
## cor(z.glm,fst) = 0.8041446
## quantiles(fst[p.glm<5e-8]) = 0.00500394 0.0781358 0.11626 0.170888 0.716683
## quantiles(fst[p.glm>5e-8]) = -0.00270008 -0.000460825 0.00616481 0.0264521 0.61094
## number p.glm < 5e-8 = 1232033
## fst threshold = 0.4051191
## number fst > threshold and meQTL = 291
## quantiles(top 50 mQTL effects) = 0.325042 0.3483962 0.3878103 0.4845272 1.005911
## ... corresponding percentiles = 0.0004612105 0.1659693 0.3864529 0.6305371 0.9876803
## Loading objects:
##       dat
## ----------------------------------------
##     CHS
## cor(z.glm,fst) = 0.8135435
## quantiles(fst[p.glm<5e-8]) = 0.00597759 0.0794297 0.117574 0.171891 0.652948
## quantiles(fst[p.glm>5e-8]) = -0.00265094 -0.00039968 0.00601137 0.0257593 0.611625
## number p.glm < 5e-8 = 1265349
## fst threshold = 0.393728
## number fst > threshold and meQTL = 340
## quantiles(top 50 mQTL effects) = 0.3410868 0.3720385 0.4206638 0.5248395 1.053211
## ... corresponding percentiles = 0.0004612105 0.1544203 0.3762065 0.6148539 0.9897495
## Loading objects:
##       dat
## ----------------------------------------
##     CLM
## cor(z.glm,fst) = 0.8135482
## quantiles(fst[p.glm<5e-8]) = 0.0224016 0.0556305 0.0796269 0.105781 0.328758
## quantiles(fst[p.glm>5e-8]) = -0.00295166 -0.00182493 0.00149705 0.01147 0.235044
## number p.glm < 5e-8 = 167505
## fst threshold = 0.20147
## number fst > threshold and meQTL = 39
## quantiles(top 50 mQTL effects) = 0.0523305 0.1183302 0.2255698 0.3139539 0.8748545
## ... corresponding percentiles = 0.01728916 0.2619655 0.5898218 0.7216491 0.9794283
## Loading objects:
##       dat
## ----------------------------------------
##     EAS
## cor(z.glm,fst) = 0.8680956
## quantiles(fst[p.glm<5e-8]) = 0.000439302 0.0397592 0.0720857 0.129926 0.775876
## quantiles(fst[p.glm>5e-8]) = -0.000639519 4.44618e-05 0.00345937 0.0116368 0.421182
## number p.glm < 5e-8 = 4376034
## fst threshold = 0.403382
## number fst > threshold and meQTL = 1046
## quantiles(top 50 mQTL effects) = 0.4849106 0.5246431 0.597621 0.7490058 1.744615
## ... corresponding percentiles = 0.0004612105 0.1845413 0.3754399 0.6227547 0.9995014
## Loading objects:
##       dat

## ----------------------------------------
##     ESN
## cor(z.glm,fst) = 0.7773469
## quantiles(fst[p.glm<5e-8]) = 0.00370301 0.0566427 0.0913793 0.153014 0.814392
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000753388 0.00513833 0.0187187 0.712695
## number p.glm < 5e-8 = 3574918
## fst threshold = 0.4717238
## number fst > threshold and meQTL = 476
## quantiles(top 50 mQTL effects) = 0.4148146 0.451976 0.5083369 0.620679 1.532532
## ... corresponding percentiles = 0.001932098 0.2440988 0.4487682 0.6648319 0.9984294
## Loading objects:
##       dat
## ----------------------------------------
##     EUR
## cor(z.glm,fst) = 0.8603913
## quantiles(fst[p.glm<5e-8]) = 0.000310886 0.0300617 0.0503502 0.0874053 0.816891
## quantiles(fst[p.glm>5e-8]) = -0.00064042 -6.40085e-05 0.00254531 0.00937606 0.288423
## number p.glm < 5e-8 = 4088735
## fst threshold = 0.2642562
## number fst > threshold and meQTL = 977
## quantiles(top 50 mQTL effects) = 0.5711828 0.6390664 0.7078051 0.8045485 1.357909
## ... corresponding percentiles = 0.0003947297 0.2030947 0.4419186 0.7026397 0.9967383
## Loading objects:
##       dat
## ----------------------------------------
##     FIN
## cor(z.glm,fst) = 0.8010559
## quantiles(fst[p.glm<5e-8]) = 0.00427405 0.0618447 0.0946288 0.134698 0.735407
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000877228 0.0055154 0.0226193 0.566318
## number p.glm < 5e-8 = 713294
## fst threshold = 0.2932548
## number fst > threshold and meQTL = 202
## quantiles(top 50 mQTL effects) = 0.2809742 0.3121592 0.3814366 0.4832262 1.038299
## ... corresponding percentiles = 0.01796228 0.2363569 0.4372172 0.6806844 0.9891512
## Loading objects:
##       dat
## ----------------------------------------
##     GBR
## cor(z.glm,fst) = 0.7985387
## quantiles(fst[p.glm<5e-8]) = 0.0227825 0.0776656 0.108308 0.148631 0.718496
## quantiles(fst[p.glm>5e-8]) = -0.0030461 -0.00107636 0.00453612 0.021334 0.447205
## number p.glm < 5e-8 = 511826
## fst threshold = 0.2976773
## number fst > threshold and meQTL = 171
## quantiles(top 50 mQTL effects) = 0.2809742 0.3433469 0.398712 0.4826658 0.8987742
## ... corresponding percentiles = 0.02026833 0.2412048 0.4687893 0.7106132 0.9814186
## Loading objects:
##       dat
## ----------------------------------------
##     GIH
## cor(z.glm,fst) = 0.7798501
## quantiles(fst[p.glm<5e-8]) = 0.0144332 0.0470435 0.0761546 0.108442 0.572979
## quantiles(fst[p.glm>5e-8]) = -0.00270008 -0.00108846 0.00334556 0.0160284 0.417406
## number p.glm < 5e-8 = 465842
## fst threshold = 0.2425602
## number fst > threshold and meQTL = 98
## quantiles(top 50 mQTL effects) = 0.1753247 0.2131099 0.3068365 0.4291675 0.9358379
## ... corresponding percentiles = 0.007799029 0.2632442 0.4723107 0.7164532 0.9840986
## Loading objects:
##       dat
## ----------------------------------------
##     GWD
## cor(z.glm,fst) = 0.7832207
## quantiles(fst[p.glm<5e-8]) = 0.00416893 0.0522506 0.0854989 0.143436 0.774292
## quantiles(fst[p.glm>5e-8]) = -0.00246896 -0.00074702 0.00433632 0.017041 0.724503
## number p.glm < 5e-8 = 3996377
## fst threshold = 0.4578288
## number fst > threshold and meQTL = 489
## quantiles(top 50 mQTL effects) = 0.4498002 0.504342 0.6025052 0.6913516 1.532532
## ... corresponding percentiles = 0.0002285277 0.2336966 0.4601593 0.6723785 0.9984294
## Loading objects:
##       dat
## ----------------------------------------
##     IBS
## cor(z.glm,fst) = 0.8119167
## quantiles(fst[p.glm<5e-8]) = 0.00334584 0.0704774 0.0974031 0.133982 0.502471
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.00102551 0.00396781 0.019292 0.551637
## number p.glm < 5e-8 = 654417
## fst threshold = 0.2748774
## number fst > threshold and meQTL = 198
## quantiles(top 50 mQTL effects) = 0.2987255 0.3401357 0.394695 0.5125635 0.8987742
## ... corresponding percentiles = 0.005310154 0.1998704 0.4637867 0.7036182 0.9814186
## Loading objects:
##       dat
## ----------------------------------------
##     ITU
## cor(z.glm,fst) = 0.7694463
## quantiles(fst[p.glm<5e-8]) = 0.00380342 0.0513284 0.079879 0.111762 0.733193
## quantiles(fst[p.glm>5e-8]) = -0.00272549 -0.00106232 0.00327329 0.0154813 0.71281
## number p.glm < 5e-8 = 442720
## fst threshold = 0.2625657
## number fst > threshold and meQTL = 95
## quantiles(top 50 mQTL effects) = 0.1431411 0.1678257 0.2084033 0.2992637 0.5270747
## ... corresponding percentiles = 0.001204964 0.1666113 0.3854266 0.5766565 0.9071346
## Loading objects:
##       dat
## ----------------------------------------
##     JPT
## cor(z.glm,fst) = 0.8129155
## quantiles(fst[p.glm<5e-8]) = 0.00421465 0.0770171 0.115577 0.171133 0.712755
## quantiles(fst[p.glm>5e-8]) = -0.00267542 -0.00039968 0.00740554 0.0280856 0.557409
## number p.glm < 5e-8 = 1323398
## fst threshold = 0.4044467
## number fst > threshold and meQTL = 332
## quantiles(top 50 mQTL effects) = 0.3383479 0.3733089 0.4137815 0.487811 0.7492857
## ... corresponding percentiles = 0.0004612105 0.1605252 0.3747003 0.6300219 0.9659535
## Loading objects:
##       dat
## ----------------------------------------
##     KHV
## cor(z.glm,fst) = 0.8034166
## quantiles(fst[p.glm<5e-8]) = 0.00342659 0.0784908 0.116719 0.17009 0.624182
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000548613 0.00631223 0.0273524 0.487064
## number p.glm < 5e-8 = 1174041
## fst threshold = 0.3839192
## number fst > threshold and meQTL = 285
## quantiles(top 50 mQTL effects) = 0.325042 0.3574118 0.4219692 0.4845272 1.031467
## ... corresponding percentiles = 0.007666067 0.1622755 0.3999734 0.6385813 0.9888603
## Loading objects:
##       dat
## ----------------------------------------
##     LWK
## cor(z.glm,fst) = 0.7715679
## quantiles(fst[p.glm<5e-8]) = 0.00177698 0.0539557 0.0874899 0.145897 0.709108
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000753388 0.00477485 0.0186299 0.650314
## number p.glm < 5e-8 = 3306296
## fst threshold = 0.4429894
## number fst > threshold and meQTL = 395
## quantiles(top 50 mQTL effects) = 0.4025743 0.4419365 0.4792388 0.6122277 1.30071
## ... corresponding percentiles = 0.004362802 0.2253533 0.4862946 0.6847917 0.9959779
## Loading objects:
##       dat
## ----------------------------------------
##     MSL
## cor(z.glm,fst) = 0.7597348
## quantiles(fst[p.glm<5e-8]) = 0.00337633 0.0594384 0.0950847 0.158405 0.804431
## quantiles(fst[p.glm>5e-8]) = -0.00325594 -0.000759799 0.00582043 0.0209458 0.744438
## number p.glm < 5e-8 = 3155828
## fst threshold = 0.4868222
## number fst > threshold and meQTL = 383
## quantiles(top 50 mQTL effects) = 0.4148146 0.4504042 0.5104922 0.6632349 1.532532
## ... corresponding percentiles = 0.004362802 0.2564372 0.4622451 0.6726028 0.9984294
## Loading objects:
##       dat
## ----------------------------------------
##     MXL
## cor(z.glm,fst) = 0.8326785
## quantiles(fst[p.glm<5e-8]) = 0.00799677 0.0753799 0.108913 0.144643 0.439354
## quantiles(fst[p.glm>5e-8]) = -0.00430149 -0.00263193 0.00241212 0.0189771 0.381654
## number p.glm < 5e-8 = 251339
## fst threshold = 0.2787501
## number fst > threshold and meQTL = 43
## quantiles(top 50 mQTL effects) = 0.0448705 0.1033102 0.158474 0.2962548 1.079595
## ... corresponding percentiles = 0.006016512 0.1942901 0.4183512 0.6992139 0.9908173
## Loading objects:
##       dat
## ----------------------------------------
##     PEL
## cor(z.glm,fst) = 0.8858191
## quantiles(fst[p.glm<5e-8]) = 0.00857764 0.0990033 0.141126 0.195231 0.619049
## quantiles(fst[p.glm>5e-8]) = -0.00325594 -0.00108409 0.00810047 0.0328757 0.379503
## number p.glm < 5e-8 = 1176600
## fst threshold = 0.3909085
## number fst > threshold and meQTL = 275
## quantiles(top 50 mQTL effects) = 0.3841934 0.4113355 0.4432731 0.5203216 1.066026
## ... corresponding percentiles = 0.005684108 0.2257002 0.4662257 0.7602682 0.9902855
## Loading objects:
##       dat
## ----------------------------------------
##     PJL
## cor(z.glm,fst) = 0.754471
## quantiles(fst[p.glm<5e-8]) = 0.0223527 0.0519653 0.080065 0.111088 0.716103
## quantiles(fst[p.glm>5e-8]) = -0.00289138 -0.0013762 0.00232414 0.0137701 0.717641
## number p.glm < 5e-8 = 244473
## fst threshold = 0.2514837
## number fst > threshold and meQTL = 61
## quantiles(top 50 mQTL effects) = 0.0974742 0.1299357 0.167993 0.2604534 0.5041718
## ... corresponding percentiles = 0.01893872 0.2045905 0.3879279 0.6465216 0.895675
## Loading objects:
##       dat
## ----------------------------------------
##     PUR
## cor(z.glm,fst) = 0.8116685
## quantiles(fst[p.glm<5e-8]) = 0.0186987 0.0496655 0.0735995 0.0966547 0.314517
## quantiles(fst[p.glm>5e-8]) = -0.00267542 -0.0018301 0.000565806 0.00814426 0.213104
## number p.glm < 5e-8 = 149601
## fst threshold = 0.1790687
## number fst > threshold and meQTL = 43
## quantiles(top 50 mQTL effects) = 0.046676 0.1027021 0.148784 0.2204236 0.8795
## ... corresponding percentiles = 0.00834334 0.191552 0.3846787 0.579216 0.9798522
## Loading objects:
##       dat
## ----------------------------------------
##     SAS
## cor(z.glm,fst) = 0.8723926
## quantiles(fst[p.glm<5e-8]) = 0.000346706 0.0291871 0.044212 0.0708242 0.569777
## quantiles(fst[p.glm>5e-8]) = -0.000653082 -0.000266788 0.00177419 0.00768108 0.310464
## number p.glm < 5e-8 = 2994159
## fst threshold = 0.2061791
## number fst > threshold and meQTL = 705
## quantiles(top 50 mQTL effects) = 0.4549423 0.4993654 0.5622248 0.7097431 1.30071
## ... corresponding percentiles = 3.32404e-05 0.2062442 0.4047351 0.6434302 0.9959779
## Loading objects:
##       dat
## ----------------------------------------
##     STU
## cor(z.glm,fst) = 0.7676886
## quantiles(fst[p.glm<5e-8]) = 0.0213233 0.0519845 0.0802204 0.110622 0.722511
## quantiles(fst[p.glm>5e-8]) = -0.00272549 -0.00106232 0.00327652 0.0154053 0.720592
## number p.glm < 5e-8 = 429097
## fst threshold = 0.2701353
## number fst > threshold and meQTL = 93
## quantiles(top 50 mQTL effects) = 0.148971 0.1732978 0.249014 0.3400254 0.5041718
## ... corresponding percentiles = 0.01051228 0.1840687 0.4177155 0.6360343 0.895675
## Loading objects:
##       dat
## ----------------------------------------
##     TSI
## cor(z.glm,fst) = 0.8116795
## quantiles(fst[p.glm<5e-8]) = 0.00967857 0.0701062 0.097481 0.134857 0.707252
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.00105311 0.00375804 0.0191462 0.578389
## number p.glm < 5e-8 = 662990
## fst threshold = 0.2851198
## number fst > threshold and meQTL = 201
## quantiles(top 50 mQTL effects) = 0.3034816 0.3429669 0.4152193 0.5308845 0.8987742
## ... corresponding percentiles = 0.006053908 0.2284072 0.4930465 0.7056355 0.9814186
## Loading objects:
##       dat
## ----------------------------------------
##     YRI
## cor(z.glm,fst) = 0.7856097
## quantiles(fst[p.glm<5e-8]) = 0.00348068 0.0548192 0.0891361 0.149388 0.751767
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.000749744 0.00480104 0.0175365 0.720375
## number p.glm < 5e-8 = 3732957
## fst threshold = 0.4711672
## number fst > threshold and meQTL = 493
## quantiles(top 50 mQTL effects) = 0.4320349 0.4653298 0.5083369 0.6272007 1.532532
## ... corresponding percentiles = 0.004362802 0.2352008 0.4414616 0.6748383 0.9984294
## Warning message:
##     In dir.create(output.dir <- "workflow/top-pairs") :
##       'workflow/top-pairs' already exists


