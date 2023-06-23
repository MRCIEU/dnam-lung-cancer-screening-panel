#!/bin/bash

#SBATCH --job-name=extract-snps
#SBATCH --account=sscm009461
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --array=1-22

CHR=chr$((SLURM_ARRAY_TASK_ID))

BASE_DIR=$SLURM_SUBMIT_DIR

GENO_DIR=$BASE_DIR/1000G
CHR_FILE=$GENO_DIR/ALL.$CHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

OUTPUT=$BASE_DIR/genotypes-$CHR.vcf.gz

module load libs/bcftools/1.8

srun bcftools view -v snps -T $BASE_DIR/snps.txt $CHR_FILE | gzip -c > $OUTPUT

