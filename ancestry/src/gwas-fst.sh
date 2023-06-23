#!/bin/bash

#SBATCH --job-name=gwas-fst
#SBATCH --account=sscm009461
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:60:00
#SBATCH --mem=16G
#SBATCH --array=0-681
  ## 31 ancestries
  ## 22 chromosomes
  ## 31 x 22 = 682

ANC_NUM=$((SLURM_ARRAY_TASK_ID / 22))
CHR_NUM=$((SLURM_ARRAY_TASK_ID % 22))

ANCESTRIES=(CLM ESN FIN GIH IBS JPT LWK MXL PJL SAS TSI EAS EUR GBR GWD ITU KHV MSL PEL PUR STU YRI)
ANCESTRY=${ANCESTRIES[$ANC_NUM]}

CHR=chr$((CHR_NUM+1))

BASE_DIR=$SLURM_SUBMIT_DIR

GENO_DIR=$BASE_DIR/1000G
CHR_FILE=$GENO_DIR/ALL.$CHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

PHENO_DIR=$BASE_DIR/pheno
PHENO_FILE=$PHENO_DIR/$ANCESTRY.txt

OUTPUT=gwas-fst/$ANCESTRY-$CHR

## bc4 plink2 installation did not support the 'fst' option 
## used local version obtained from the plink2 website
## https://www.cog-genomics.org/plink/2.0/
## https://s3.amazonaws.com/plink2-assets/plink2_linux_x86_64_20230607.zip
srun ./plink2 --vcf $CHR_FILE --pheno $PHENO_FILE --fst phenotype report-variants --1 --out $OUTPUT
