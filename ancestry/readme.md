# DNAm sites associated with ancestry

Sites will be selected by first running
a GWAS of each ancestry in 1000 Genomes (https://www.internationalgenome.org/).
For each ancestry, associated genetic variants that are also mQTLs (http://mqtldb.godmc.org.uk/)
will be ordered by mQTL effect size
and the top 50 selected for the panel.

## Download and prepare GoDMC

* Input:
  - mQTL summary statistics from GoDMC http://fileserve.mrcieu.ac.uk/mqtl/assoc_meta_all.csv.gz
  - mapping between rsid and genomic coordinates ftp://hgdownload.soe.ucsc.edu/goldenPath/hg19/database/snp151Common.txt.gz

* Output: GoDMC summary statistics with hg38 coordinates in `godmc-hg38.csv.gz`

```
Rscript src/extract-mqtls.r godmc-hg38.csv.gz
```

## Prepare to run ancestry GWAS

* Input: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
* Output: GWAS phenotype files for each ancestry in `pheno/`

'Super-populations' (i.e. AFR AMR EAS EUR SAS) will each be compared
to all other super-populations.
Other populations will be compared to all other populations
within the *same* super-population.

On the compute cluster:
1. create a folder for the scripts, inputs and outputs (in scratch space),
    call it `BASE`
2. copy `ancestry/src/` to `BASE`
3. generate GWAS phenotype files as follows:

```
Rscript src/generate-phenotype-files.r pheno
```

## Perform GWAS for each ancestry

* Input: `pheno/`, 1000 Genomes data
* Output: GWAS outputs for each ancestry in `gwas-fst/` and `gwas-glm/` 

On the compute cluster:
1. copy `ancestry/ancestries.txt` to `BASE`
2. download 1000 Genomes data (http://hgdownload.cse.ucsc.edu/gbdb/hg38/1000Genomes/) to `BASE/1000G`
3. submit the GWAS jobs to the system as follows:

```
sbatch src/gwas-glm.sh
sbatch src/gwas-fst.sh
```

## Select top ancestry mQTLs

* Input: `gwas-fst/`, `gwas-glm/`, `godmc-hg38.csv.gz`, `ancestries.txt`
* Output: Up to 50 DNAm sites associated with each ancestry in `sites/sites.csv`.

On the compute cluster:
1. copy `godmc-hg38.csv.gz` to `BASE`
2. submit the jobs to cluster as follows:

```
sbatch src/filter-sites.sh gwas-glm gwas-fst godmc-hg38.csv.gz sites
```

Afterward, collate the top 50 for each ancestry into a single file:

```
Rscript src/select-sites.r sites sites/sites.csv
```

## Check ancestry mQTLs

Determine to what extent the SNPs identified actually
capture genetic variation of ancestry.

### Extract mQTL genotypes from 1000 Genomes

* Input: `sites/sites.csv`, 1000 Genomes data
* Output: Genotype vcf files in `genotypes/` for each mQTL in `sites.csv`

Submit the jobs to extract genotypes to the system as follows:

```
Rscript src/extract-sites.r sites/sites.csv sites.txt
sbatch extract-genotypes.sh
```

### Compare genetic clusters to ancestry

Plot principal components of selected mQTLs and compare to ancestry.

* Input: `sites/sites.csv` and `genotypes/*.vcf.gz` 
* Output: `pca-of-genotype.pdf`

```
Rscript src/check-sites.r \
  sites/sites.csv \
  genotypes \
  pca-of-genotype.pdf
```
