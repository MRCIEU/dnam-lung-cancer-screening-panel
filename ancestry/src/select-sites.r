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
#' the GWAS Fst > 97.5th percentile of meQTLs
#' with logistic GWAS p < 5e-8.
#'
#' (about 5 minutes per ancestry, total is 2-2.5 hours runtime)

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

    ## Fst threshold at 97.5% percentile of Fst for mQTLs
    ## with logistic p < 5e-8
    threshold <- with(dat, quantile(HUDSON_FST[p.glm < 5e-8 & !is.na(beta.godmc)], probs=0.975, na.rm=T))
    
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
## Wed Jun 21 15:57:28 2023 ACB 
## cor(z.glm,fst) = 0.8057804 
## quantiles(fst[p.glm<5e-8]) = -0.000911166 0.0610505 0.0923413 0.1414502 0.602116 
## quantiles(fst[p.glm>5e-8]) = -0.00289138 -0.00106996 0.00392139 0.0172468 0.361634 
## number p.glm < 5e-8 = 2085812 
## fst threshold = 0.324955 
## number fst > threshold and meQTL = 756 
## quantiles(top 50 mQTL effects) = 0.4881662 0.5160562 0.6154027 0.8273261 1.294406 
## ... corresponding percentiles = 0.0004213279 0.2333103 0.4566465 0.7125629 0.9974523 
## ----------------------------------------
## Wed Jun 21 16:04:24 2023 AFR 
## cor(z.glm,fst) = 0.9090676 
## quantiles(fst[p.glm<5e-8]) = 0.000186468 0.0394018 0.0710868 0.132892 0.897636 
## quantiles(fst[p.glm>5e-8]) = -0.000519405 0.000501096 0.00463758 0.012219 0.23092 
## number p.glm < 5e-8 = 8923152 
## fst threshold = 0.432723 
## number fst > threshold and meQTL = 3435 
## quantiles(top 50 mQTL effects) = 0.7621843 0.8448877 0.9680321 1.090138 1.424667 
## ... corresponding percentiles = 0.0002244457 0.2180608 0.4564871 0.6954587 0.9985824 
## ----------------------------------------
## Wed Jun 21 16:09:00 2023 AMR 
## cor(z.glm,fst) = 0.9106141 
## quantiles(fst[p.glm<5e-8]) = 0.00342562 0.03572135 0.0498023 0.0707248 0.378841 
## quantiles(fst[p.glm>5e-8]) = -0.000872405 -0.00039968 0.00157086 0.00758021 0.194131 
## number p.glm < 5e-8 = 1886375 
## fst threshold = 0.1472127 
## number fst > threshold and meQTL = 1196 
## quantiles(top 50 mQTL effects) = 0.521535 0.5635022 0.6319406 0.7726194 1.419211 
## ... corresponding percentiles = 5.906465e-05 0.195822 0.4157443 0.6910239 0.9985313 
## ----------------------------------------
## Wed Jun 21 16:13:50 2023 ASW 
## cor(z.glm,fst) = 0.7737309 
## quantiles(fst[p.glm<5e-8]) = 0.00725419 0.0805743 0.112407 0.157941 0.495862 
## quantiles(fst[p.glm>5e-8]) = -0.00451006 -0.00226252 0.00322114 0.0183232 0.343529 
## number p.glm < 5e-8 = 512837 
## fst threshold = 0.2987961 
## number fst > threshold and meQTL = 146 
## quantiles(top 50 mQTL effects) = 0.1941827 0.2315371 0.3411171 0.4430264 1.294406 
## ... corresponding percentiles = 0.004252655 0.2491416 0.4214145 0.6368882 0.9974523 
## ----------------------------------------
## Wed Jun 21 16:18:19 2023 BEB 
## cor(z.glm,fst) = 0.7509986 
## quantiles(fst[p.glm<5e-8]) = 0.0238947 0.0554142 0.083516 0.113937 0.705645 
## quantiles(fst[p.glm>5e-8]) = -0.00321903 -0.00166169 0.00185287 0.013371 0.7256 
## number p.glm < 5e-8 = 182463 
## fst threshold = 0.226581 
## number fst > threshold and meQTL = 114 
## quantiles(top 50 mQTL effects) = 0.1514779 0.1830882 0.2508346 0.3541425 0.6771178 
## ... corresponding percentiles = 0.001622309 0.166648 0.3678704 0.6018914 0.9618639 
## ----------------------------------------
## Wed Jun 21 16:22:33 2023 CDX 
## cor(z.glm,fst) = 0.7987636 
## quantiles(fst[p.glm<5e-8]) = 0.00424579 0.08025 0.12114 0.178094 0.850444 
## quantiles(fst[p.glm>5e-8]) = -0.00298217 -0.00039968 0.00737566 0.0293675 0.545763 
## number p.glm < 5e-8 = 1179326 
## fst threshold = 0.3551125 
## number fst > threshold and meQTL = 782 
## quantiles(top 50 mQTL effects) = 0.4422355 0.4774985 0.5147881 0.6164063 1.269122 
## ... corresponding percentiles = 0.003902205 0.188555 0.3885509 0.6327242 0.9971176 
## ----------------------------------------
## Wed Jun 21 16:26:46 2023 CEU 
## cor(z.glm,fst) = 0.8022591 
## quantiles(fst[p.glm<5e-8]) = 0.00738666 0.0738138 0.103041 0.142041 0.723754 
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.00102915 0.00424947 0.0203776 0.539383 
## number p.glm < 5e-8 = 595692 
## fst threshold = 0.252837 
## number fst > threshold and meQTL = 513 
## quantiles(top 50 mQTL effects) = 0.4743543 0.4975288 0.6019193 0.7323469 1.247426 
## ... corresponding percentiles = 0.0003898267 0.2594198 0.4760729 0.736237 0.9968735 
## ----------------------------------------
## Wed Jun 21 16:31:05 2023 CHB 
## cor(z.glm,fst) = 0.8041446 
## quantiles(fst[p.glm<5e-8]) = 0.00500394 0.0781358 0.11626 0.170888 0.716683 
## quantiles(fst[p.glm>5e-8]) = -0.00270008 -0.000460825 0.00616481 0.0264521 0.61094 
## number p.glm < 5e-8 = 1232033 
## fst threshold = 0.3495569 
## number fst > threshold and meQTL = 799 
## quantiles(top 50 mQTL effects) = 0.4712405 0.4917546 0.5417302 0.6639378 1.031467 
## ... corresponding percentiles = 0.0002086951 0.1691454 0.3671616 0.6289598 0.9920145 
## ----------------------------------------
## Wed Jun 21 16:35:14 2023 CHS 
## cor(z.glm,fst) = 0.8135435 
## quantiles(fst[p.glm<5e-8]) = 0.00597759 0.0794297 0.117574 0.171891 0.652948 
## quantiles(fst[p.glm>5e-8]) = -0.00265094 -0.00039968 0.00601137 0.0257593 0.611625 
## number p.glm < 5e-8 = 1265349 
## fst threshold = 0.3438824 
## number fst > threshold and meQTL = 909 
## quantiles(top 50 mQTL effects) = 0.457947 0.4889516 0.54803 0.6499327 1.053211 
## ... corresponding percentiles = 0.000452829 0.1744061 0.3671616 0.6256758 0.992735 
## ----------------------------------------
## Wed Jun 21 16:39:30 2023 CLM 
## cor(z.glm,fst) = 0.8135482 
## quantiles(fst[p.glm<5e-8]) = 0.0224016 0.0556305 0.0796269 0.105781 0.328758 
## quantiles(fst[p.glm>5e-8]) = -0.00295166 -0.00182493 0.00149705 0.01147 0.235044 
## number p.glm < 5e-8 = 167505 
## fst threshold = 0.1833278 
## number fst > threshold and meQTL = 106 
## quantiles(top 50 mQTL effects) = 0.2040587 0.2316076 0.2771103 0.4217372 0.8757084 
## ... corresponding percentiles = 0.0003898267 0.1973232 0.5429203 0.6979148 0.9842179 
## ----------------------------------------
## Wed Jun 21 16:43:39 2023 EAS 
## cor(z.glm,fst) = 0.8680956 
## quantiles(fst[p.glm<5e-8]) = 0.000439302 0.0397592 0.0720857 0.129926 0.775876 
## quantiles(fst[p.glm>5e-8]) = -0.000639519 4.44618e-05 0.00345937 0.0116368 0.421182 
## number p.glm < 5e-8 = 4376034 
## fst threshold = 0.3378494 
## number fst > threshold and meQTL = 2789 
## quantiles(top 50 mQTL effects) = 0.6412484 0.6993055 0.8086507 0.9451949 1.744615 
## ... corresponding percentiles = 0.0002086951 0.180222 0.3892557 0.6569564 0.9997913 
## ----------------------------------------
## Wed Jun 21 16:47:57 2023 ESN 
## cor(z.glm,fst) = 0.7773469 
## quantiles(fst[p.glm<5e-8]) = 0.00370301 0.0566427 0.0913793 0.153014 0.814392 
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000753388 0.00513833 0.0187187 0.712695 
## number p.glm < 5e-8 = 3574918 
## fst threshold = 0.4192236 
## number fst > threshold and meQTL = 1224 
## quantiles(top 50 mQTL effects) = 0.5135676 0.5810549 0.6254706 0.845682 1.294406 
## ... corresponding percentiles = 0.0002244457 0.2124014 0.4359759 0.6824556 0.9974523 
## ----------------------------------------
## Wed Jun 21 16:52:16 2023 EUR 
## cor(z.glm,fst) = 0.8603913 
## quantiles(fst[p.glm<5e-8]) = 0.000310886 0.0300617 0.0503502 0.0874053 0.816891 
## quantiles(fst[p.glm>5e-8]) = -0.00064042 -6.40085e-05 0.00254531 0.00937606 0.288423 
## number p.glm < 5e-8 = 4088735 
## fst threshold = 0.2196835 
## number fst > threshold and meQTL = 2534 
## quantiles(top 50 mQTL effects) = 0.7419954 0.8153078 0.9170538 1.033946 1.433932 
## ... corresponding percentiles = 0.0001181293 0.2012658 0.4519903 0.7175902 0.9986612 
## ----------------------------------------
## Wed Jun 21 16:56:31 2023 FIN 
## cor(z.glm,fst) = 0.8010559 
## quantiles(fst[p.glm<5e-8]) = 0.00427405 0.0618447 0.0946288 0.134698 0.735407 
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000877228 0.0055154 0.0226193 0.566318 
## number p.glm < 5e-8 = 713294 
## fst threshold = 0.2507255 
## number fst > threshold and meQTL = 525 
## quantiles(top 50 mQTL effects) = 0.4743543 0.4997518 0.6134008 0.7637683 1.357909 
## ... corresponding percentiles = 0.0003898267 0.2516036 0.4907249 0.7328348 0.9980509 
## ----------------------------------------
## Wed Jun 21 17:00:56 2023 GBR 
## cor(z.glm,fst) = 0.7985387 
## quantiles(fst[p.glm<5e-8]) = 0.0227825 0.0776656 0.108308 0.148631 0.718496 
## quantiles(fst[p.glm>5e-8]) = -0.0030461 -0.00107636 0.00453612 0.021334 0.447205 
## number p.glm < 5e-8 = 511826 
## fst threshold = 0.2623679 
## number fst > threshold and meQTL = 455 
## quantiles(top 50 mQTL effects) = 0.4443228 0.4887278 0.5618418 0.7328116 1.357909 
## ... corresponding percentiles = 0.005930091 0.23137 0.4796562 0.7373277 0.9980509 
## ----------------------------------------
## Wed Jun 21 17:05:19 2023 GIH 
## cor(z.glm,fst) = 0.7798501 
## quantiles(fst[p.glm<5e-8]) = 0.0144332 0.0470435 0.0761546 0.108442 0.572979 
## quantiles(fst[p.glm>5e-8]) = -0.00270008 -0.00108846 0.00334556 0.0160284 0.417406 
## number p.glm < 5e-8 = 465842 
## fst threshold = 0.2112379 
## number fst > threshold and meQTL = 286 
## quantiles(top 50 mQTL effects) = 0.3648004 0.410123 0.478124 0.6018003 1.30071 
## ... corresponding percentiles = 3.937643e-05 0.2397759 0.4873641 0.7207561 0.9975311 
## ----------------------------------------
## Wed Jun 21 17:09:29 2023 GWD 
## cor(z.glm,fst) = 0.7832207 
## quantiles(fst[p.glm<5e-8]) = 0.00416893 0.0522506 0.0854989 0.143436 0.774292 
## quantiles(fst[p.glm>5e-8]) = -0.00246896 -0.00074702 0.00433632 0.017041 0.724503 
## number p.glm < 5e-8 = 3996377 
## fst threshold = 0.4085866 
## number fst > threshold and meQTL = 1298 
## quantiles(top 50 mQTL effects) = 0.5264864 0.6143798 0.6792073 0.8941232 1.294406 
## ... corresponding percentiles = 0.0002244457 0.2303265 0.4610921 0.6852248 0.9974523 
## ----------------------------------------
## Wed Jun 21 17:13:50 2023 IBS 
## cor(z.glm,fst) = 0.8119167 
## quantiles(fst[p.glm<5e-8]) = 0.00334584 0.0704774 0.0974031 0.133982 0.502471 
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.00102551 0.00396781 0.019292 0.551637 
## number p.glm < 5e-8 = 654417 
## fst threshold = 0.2401014 
## number fst > threshold and meQTL = 520 
## quantiles(top 50 mQTL effects) = 0.4302148 0.4751803 0.5317872 0.647932 1.357909 
## ... corresponding percentiles = 0.0003898267 0.1868215 0.4490961 0.7066849 0.9980509 
## ----------------------------------------
## Wed Jun 21 17:18:09 2023 ITU 
## cor(z.glm,fst) = 0.7694463 
## quantiles(fst[p.glm<5e-8]) = 0.00380342 0.0513284 0.079879 0.111762 0.733193 
## quantiles(fst[p.glm>5e-8]) = -0.00272549 -0.00106232 0.00327329 0.0154813 0.71281 
## number p.glm < 5e-8 = 442720 
## fst threshold = 0.2201977 
## number fst > threshold and meQTL = 263 
## quantiles(top 50 mQTL effects) = 0.2701145 0.3054972 0.3613085 0.4410852 0.6771178 
## ... corresponding percentiles = 0.001173418 0.1863313 0.3991629 0.6285109 0.9618639 
## ----------------------------------------
## Wed Jun 21 17:22:23 2023 JPT 
## cor(z.glm,fst) = 0.8129155 
## quantiles(fst[p.glm<5e-8]) = 0.00421465 0.0770171 0.115577 0.171133 0.712755 
## quantiles(fst[p.glm>5e-8]) = -0.00267542 -0.00039968 0.00740554 0.0280856 0.557409 
## number p.glm < 5e-8 = 1323398 
## fst threshold = 0.348451 
## number fst > threshold and meQTL = 902 
## quantiles(top 50 mQTL effects) = 0.4385042 0.4861561 0.5521897 0.6476541 1.247972 
## ... corresponding percentiles = 0.000452829 0.1679139 0.3700853 0.626107 0.9968853 
## ----------------------------------------
## Wed Jun 21 17:26:41 2023 KHV 
## cor(z.glm,fst) = 0.8034166 
## quantiles(fst[p.glm<5e-8]) = 0.00342659 0.0784908 0.116719 0.17009 0.624182 
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000548613 0.00631223 0.0273524 0.487064 
## number p.glm < 5e-8 = 1174041 
## fst threshold = 0.335856 
## number fst > threshold and meQTL = 764 
## quantiles(top 50 mQTL effects) = 0.4305975 0.4721697 0.5097326 0.5784254 1.031467 
## ... corresponding percentiles = 0.0002086951 0.167089 0.3747002 0.6389388 0.9920145 
## ----------------------------------------
## Wed Jun 21 17:30:55 2023 LWK 
## cor(z.glm,fst) = 0.7715679 
## quantiles(fst[p.glm<5e-8]) = 0.00177698 0.0539557 0.0874899 0.145897 0.709108 
## quantiles(fst[p.glm>5e-8]) = -0.00280661 -0.000753388 0.00477485 0.0186299 0.650314 
## number p.glm < 5e-8 = 3306296 
## fst threshold = 0.3954229 
## number fst > threshold and meQTL = 1062 
## quantiles(top 50 mQTL effects) = 0.5083987 0.6041707 0.641942 0.814583 1.30071 
## ... corresponding percentiles = 0.0002244457 0.2055588 0.4607889 0.6940549 0.9975311 
## ----------------------------------------
## Wed Jun 21 17:35:15 2023 MSL 
## cor(z.glm,fst) = 0.7597348 
## quantiles(fst[p.glm<5e-8]) = 0.00337633 0.0594384 0.0950847 0.158405 0.804431 
## quantiles(fst[p.glm>5e-8]) = -0.00325594 -0.000759799 0.00582043 0.0209458 0.744438 
## number p.glm < 5e-8 = 3155828 
## fst threshold = 0.4330515 
## number fst > threshold and meQTL = 1010 
## quantiles(top 50 mQTL effects) = 0.4979061 0.533725 0.6353121 0.8445124 1.424667 
## ... corresponding percentiles = 0.0002244457 0.2412752 0.4748503 0.6919542 0.9985824 
## ----------------------------------------
## Wed Jun 21 17:39:33 2023 MXL 
## cor(z.glm,fst) = 0.8326785 
## quantiles(fst[p.glm<5e-8]) = 0.00799677 0.0753799 0.108913 0.144643 0.439354 
## quantiles(fst[p.glm>5e-8]) = -0.00430149 -0.00263193 0.00241212 0.0189771 0.381654 
## number p.glm < 5e-8 = 251339 
## fst threshold = 0.251976 
## number fst > threshold and meQTL = 115 
## quantiles(top 50 mQTL effects) = 0.1680102 0.2109927 0.3332549 0.4119499 1.127537 
## ... corresponding percentiles = 0.005894652 0.1961084 0.4265255 0.687613 0.9946921 
## ----------------------------------------
## Wed Jun 21 17:43:54 2023 PEL 
## cor(z.glm,fst) = 0.8858191 
## quantiles(fst[p.glm<5e-8]) = 0.00857764 0.0990033 0.141126 0.195231 0.619049 
## quantiles(fst[p.glm>5e-8]) = -0.00325594 -0.00108409 0.00810047 0.0328757 0.379503 
## number p.glm < 5e-8 = 1176600 
## fst threshold = 0.3481224 
## number fst > threshold and meQTL = 730 
## quantiles(top 50 mQTL effects) = 0.5092973 0.5492072 0.6246609 0.7573736 1.160249 
## ... corresponding percentiles = 0.001712875 0.2165871 0.4603676 0.7524295 0.9954087 
## ----------------------------------------
## Wed Jun 21 17:48:03 2023 PJL 
## cor(z.glm,fst) = 0.754471 
## quantiles(fst[p.glm<5e-8]) = 0.0223527 0.0519653 0.080065 0.111088 0.716103 
## quantiles(fst[p.glm>5e-8]) = -0.00289138 -0.0013762 0.00232414 0.0137701 0.717641 
## number p.glm < 5e-8 = 244473 
## fst threshold = 0.2108025 
## number fst > threshold and meQTL = 173 
## quantiles(top 50 mQTL effects) = 0.2403384 0.2784439 0.3915556 0.4541523 1.30071 
## ... corresponding percentiles = 0.001173418 0.2166413 0.422773 0.673715 0.9975311 
## ----------------------------------------
## Wed Jun 21 17:52:18 2023 PUR 
## cor(z.glm,fst) = 0.8116685 
## quantiles(fst[p.glm<5e-8]) = 0.0186987 0.0496655 0.0735995 0.0966547 0.314517 
## quantiles(fst[p.glm>5e-8]) = -0.00267542 -0.0018301 0.000565806 0.00814426 0.213104 
## number p.glm < 5e-8 = 149601 
## fst threshold = 0.1617316 
## number fst > threshold and meQTL = 116 
## quantiles(top 50 mQTL effects) = 0.1778718 0.2120857 0.2810823 0.3518412 1.295521 
## ... corresponding percentiles = 0.008166673 0.2033645 0.4438413 0.644229 0.9974642 
## ----------------------------------------
## Wed Jun 21 17:56:28 2023 SAS 
## cor(z.glm,fst) = 0.8723926 
## quantiles(fst[p.glm<5e-8]) = 0.000346706 0.0291871 0.044212 0.0708242 0.569777 
## quantiles(fst[p.glm>5e-8]) = -0.000653082 -0.000266788 0.00177419 0.00768108 0.310464 
## number p.glm < 5e-8 = 2994159 
## fst threshold = 0.168248 
## number fst > threshold and meQTL = 1874 
## quantiles(top 50 mQTL effects) = 0.6256475 0.6671933 0.6996167 0.8140509 1.309731 
## ... corresponding percentiles = 3.937643e-05 0.2038981 0.4352868 0.6852415 0.9976453 
## ----------------------------------------
## Wed Jun 21 18:00:45 2023 STU 
## cor(z.glm,fst) = 0.7676886 
## quantiles(fst[p.glm<5e-8]) = 0.0213233 0.0519845 0.0802204 0.110622 0.722511 
## quantiles(fst[p.glm>5e-8]) = -0.00272549 -0.00106232 0.00327652 0.0154053 0.720592 
## number p.glm < 5e-8 = 429097 
## fst threshold = 0.221008 
## number fst > threshold and meQTL = 251 
## quantiles(top 50 mQTL effects) = 0.2751712 0.3286895 0.3972233 0.4535075 1.196595 
## ... corresponding percentiles = 0.00143724 0.1791529 0.4206585 0.6564229 0.9961135 
## ----------------------------------------
## Wed Jun 21 18:05:08 2023 TSI 
## cor(z.glm,fst) = 0.8116795 
## quantiles(fst[p.glm<5e-8]) = 0.00967857 0.0701062 0.097481 0.134857 0.707252 
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.00105311 0.00375804 0.0191462 0.578389 
## number p.glm < 5e-8 = 662990 
## fst threshold = 0.2467944 
## number fst > threshold and meQTL = 535 
## quantiles(top 50 mQTL effects) = 0.4534575 0.4978016 0.5856063 0.6745594 1.038299 
## ... corresponding percentiles = 0.0003898267 0.2071653 0.4680795 0.7144263 0.9922547 
## ----------------------------------------
## Wed Jun 21 18:09:27 2023 YRI 
## cor(z.glm,fst) = 0.7856097 
## quantiles(fst[p.glm<5e-8]) = 0.00348068 0.0548192 0.0891361 0.149388 0.751767 
## quantiles(fst[p.glm>5e-8]) = -0.0026025 -0.000749744 0.00480104 0.0175365 0.720375 
## number p.glm < 5e-8 = 3732957 
## fst threshold = 0.4167254 
## number fst > threshold and meQTL = 1301 
## quantiles(top 50 mQTL effects) = 0.5135676 0.5810549 0.6556206 0.8089115 1.30071 
## ... corresponding percentiles = 0.0002244457 0.2265169 0.4501357 0.6835749 0.9975311 



