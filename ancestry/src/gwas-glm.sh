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

ANCESTRIES_FILE=$1
GENO_DIR=$2
PHENO_DIR=$3
OUTPUT_DIR=$4

readarray -t ANCESTRIES < $ANCESTRIES_FILE

ANCESTRY=${ANCESTRIES[$ANC_NUM]}

CHR=chr$((CHR_NUM+1))

cd $SLURM_SUBMIT_DIR

CHR_FILE=$GENO_DIR/ALL.$CHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

PHENO_FILE=$PHENO_DIR/$ANCESTRY.txt

OUTPUT=$OUTPUT_DIR/$ANCESTRY-$CHR

module load apps/plink/2

srun plink2 --vcf $CHR_FILE --pheno $PHENO_FILE --glm --1 --out $OUTPUT
