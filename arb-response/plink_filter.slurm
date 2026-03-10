#!/bin/bash

#SBATCH --job-name=plink-filter
#SBATCH --output=plink-filter.log

#SBATCH --ntasks=5
#SBATCH --time=5:00:00

#SBATCH --mem=36G
#SBATCH --no-requeue

#SBATCH --export=NONE
unset SLURM_EXPORT_ENV

export OMP_NUM_THREADS=1

module load plink/1.90

plink --memory 20000 --threads 4 --bfile autosomes --extract snps_to_keep.txt --geno 1 --make-bed --out filtered_for_arb_response


