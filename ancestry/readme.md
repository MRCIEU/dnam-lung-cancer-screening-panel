# DNAm sites associated with ancestry

Sites will be selected by first running
a GWAS of each ancestry in 1000 Genomes (https://www.internationalgenome.org/).
For each ancestry, associated genetic variants that are also mQTLs (http://mqtldb.godmc.org.uk/)
will be ordered by mQTL effect size
and the top 50 selected for the panel.

## Download GoDMC

http://mqtldb.godmc.org.uk/downloads

Save to the `godmc/` folder.

[Convert to hg38 coordinates]

## Create phenotype files for GWAS

Input: http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/integrated_call_samples_v3.20130502.ALL.panel
Output: `output/pheno/*.txt`
```
Rscript src/generate-phenotype-files.r output/pheno
```

## Perform GWAS for each ancestry

On the compute cluster:
1. create a folder for the scripts, inputs and outputs (in scratch space), call it `BASE`
2. copy `output/pheno` folder to `BASE`
3. copy `src/gwas-*.sh` scripts to `BASE`
4. set `BASE_DIR` in the scripts to the value of `BASE`
5. download 1000 Genomes data (http://ftp.1000genomes.ebi.ac.uk/vol1/ftp/release/20130502/) to `BASE/1000G`
6. submit the GWAS jobs to the system as follows:
```
sbatch src/gwas-glm.sh
sbatch src/gwas-fst.sh
```

When finished, copy the `gwas-fst` and `gwas-glm` output folders
to the `output` folder for panel selection in the next step.

## Select top ancestry mQTLs

```
Rscript src/select-sites.r \
  output/pheno \
  output/gwas-glm \
  output/gwas-fst \
  godmc/godmc-hg38.txt \
  output/panel-sites
```
