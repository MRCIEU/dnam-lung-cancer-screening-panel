# DNAm sites associated with ancestry

Sites will be selected by first running
a GWAS of each ancestry in 1000 Genomes (https://www.internationalgenome.org/).
For each ancestry, associated genetic variants that are also mQTLs (http://mqtldb.godmc.org.uk/)
will be ordered by mQTL effect size
and the top 50 selected for the panel.

## Create phenotype files for GWAS

* Input: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
* Output: GWAS phenotype files for each ancestry in `output/pheno/`

'Super-populations' (i.e. AFR AMR EAS EUR SAS) are each compared to all others.
Other populations are compared to all other populations within the *same* super-population.

```
Rscript src/generate-phenotype-files.r output/pheno
```

## Perform GWAS for each ancestry

* Input: `output/pheno/*.txt`, 1000 Genomes data
* Output: GWAS outputs for each ancestry in `gwas-fst/` and `gwas-glm/` 

On the compute cluster:
1. create a folder for the scripts, inputs and outputs (in scratch space), call it `BASE`
2. copy `output/pheno` folder to `BASE`
3. copy `src/gwas-*.sh` scripts to `BASE`
4. set `BASE_DIR` in the scripts to the value of `BASE`
5. download 1000 Genomes data (http://hgdownload.cse.ucsc.edu/gbdb/hg38/1000Genomes/) to `BASE/1000G`
6. submit the GWAS jobs to the system as follows:

```
sbatch src/gwas-glm.sh
sbatch src/gwas-fst.sh
```

When finished, copy the `gwas-fst` and `gwas-glm` output folders
to the `output` folder for panel selection in the next step.

## Download GoDMC

* Input:
  - mQTL summary statistics from GoDMC http://fileserve.mrcieu.ac.uk/mqtl/assoc_meta_all.csv.gz
  - mapping between rsid and genomic coordinates ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz

* Output: GoDMC summary statistics with hg38 coordinates (`output/godmc-hg38.csv.gz`)

```
Rscript src/extract-mqtls.r output/godmc-hg38.csv.gz
```

## Select top ancestry mQTLs

* Input: `output/pheno`, `output/gwas-fst/`, `output/gwas-glm`, `output/godmc-hg38.csv.gz`
* Output: Up to 50 DNAm sites associated with each ancestry in `output/panel-sites`

```
Rscript src/select-sites.r \
  output/pheno \
  output/gwas-glm \
  output/gwas-fst \
  output/godmc-hg38.csv.gz \
  output/panel-sites
```

## Check ancestry mQTLs

Determine to what extent the SNPs identified actually
capture genetic variation of ancestry.

### Extract mQTL genotypes from 1000 Genomes

* Input: `output/panel-sites/panel-sites.csv`, 1000 Genomes data
* Output: Genotypes for each mQTL in `panel-sites.csv`

1. create a file `snps.txt` containing the mQTL genomic coordinates
```r
panel.sites <- fread("output/panel-sites/panel-sites.csv")
fwrite(unique(panel.sites[,c("chr","pos")]),
       file="snps.txt", sep="\t",col.names=F)
```

On the compute cluster:
2. create a folder for the scripts, inputs and outputs (in scratch space), call it `BASE`
3. copy `snps.txt` to `BASE`
4. copy `src/extract-genotypes.sh` scripts to `BASE`
5. set `BASE_DIR` in the scripts to the value of `BASE`
6. download 1000 Genomes data (http://hgdownload.cse.ucsc.edu/gbdb/hg38/1000Genomes/) to `BASE/1000G`
7. submit the jobs to the system as follows:

```
sbatch src/extract-genotypes.sh
```

Once started, each job will take ~5 minutes (~90 minutes if run sequentially).

When finished, copy the `genotypes-chr*.vcf.gz` output files
to the `output/genotypes` folder.

### Compare genetic clusters to ancestry

Plot principal components of selected mQTLs and compare to ancestry.

* Input: `output/panel-sites/panel-sites.csv` and `output/genotypes/*.vcf.gz` 
* Output: `output/pca-of-genotype.pdf`

```
Rscript src/check-sites.r \
  output/panel-sites/panel-sites.csv \
  output/genotypes \
  output/pca-of-genotype.pdf
```
