#!/bin/bash

#SBATCH --job-name=extract-genotypes
#SBATCH --account=sscm009461
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=0:30:00
#SBATCH --mem=8G
#SBATCH --array=1-22

CHR=chr$((SLURM_ARRAY_TASK_ID))

SITES_FILENAME=$1
GENO_DIR=$2
OUTPUT_DIR=$3

cd $SLURM_SUBMIT_DIR

CHR_FILE=$GENO_DIR/ALL.$CHR.shapeit2_integrated_v1a.GRCh38.20181129.phased.vcf.gz

OUTPUT=$OUTPUT_DIR/genotypes-$CHR.vcf.gz

module load libs/bcftools/1.8

srun bcftools view -v snps -T $SITES_FILENAME $CHR_FILE | gzip -c > $OUTPUT

