#!/bin/bash

#SBATCH --job-name=gwas-glm
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

ANCESTRIES=(CLM ESN FIN GIH IBS JPT LWK MXL PJL SAS TSI EAS EUR GBR GWD ITU KHV MSL PEL PUR STU YRI ACB AFR AMR ASW BEB CDX CEU CHB CHS)

ANCESTRY=${ANCESTRIES[$ANC_NUM]}

CHR=chr$((CHR_NUM+1))

cd $SLURM_SUBMIT_DIR

GENO_DIR=1000G
CHR_FILE=$GENO_DIR/ALL.$CHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

PHENO_DIR=pheno
PHENO_FILE=$PHENO_DIR/$ANCESTRY.txt

OUTPUT=gwas-glm/$ANCESTRY-$CHR

module load apps/plink/2

srun plink2 --vcf $CHR_FILE --pheno $PHENO_FILE --glm --1 --out $OUTPUT
