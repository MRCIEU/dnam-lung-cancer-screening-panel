#!/bin/bash

#SBATCH --job-name=filter-sites
#SBATCH --account=sscm009461
#SBATCH --partition=mrcieu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task=1
#SBATCH --time=00:60:00
#SBATCH --mem=16G
#SBATCH --array=0-1
  ## 0-30 31 ancestries

readarray -t ANCESTRIES < ancestries.txt

ANCESTRY=${ANCESTRIES[$SLURM_ARRAY_TASK_ID]}

cd $SLURM_SUBMIT_DIR

Rscript src/filter-sites.r $ANCESTRY $1 $2 $3 $4 $5
