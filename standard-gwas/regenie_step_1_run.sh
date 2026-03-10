#!/usr/bin/env bash
#SBATCH --job-name=rg_step1_rows
#SBATCH --array=11-368
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=20
#SBATCH --mem=64G
#SBATCH --no-requeue
#SBATCH --output=logs/step1_%A_%a.out
#SBATCH --error=logs/step1_%A_%a.err

set -euo pipefail
mkdir -p logs

CONDA_ENV="regenie412"
BED_STEP1="<BED_STEP1>"
OUTROOT="<OUTROOT>"
STEP1DIR="${OUTROOT}/step1"
TMPROOT="${OUTROOT}/regenie_tmp"
PAIRS="${OUTROOT}/pairs_all.tsv"
mkdir -p "${STEP1DIR}" "${TMPROOT}"

module load conda
conda activate "${CONDA_ENV}"
REG_BIN="$(which regenie)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

LINE="$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "${PAIRS}")"
IFS=$'\t' read -r RUN_ID SCOPE TRAIT PHENOTYPE COVARS_IN_MODEL PHENO_PATH COVAR_PATH INTERACTION_VAR <<< "${LINE}"

TRAIT_FLAG="--qt"
[[ "${TRAIT}" == "CVD" || "${TRAIT}" =~ ^CLASS_ ]] && TRAIT_FLAG="--bt"

OUTDIR="${STEP1DIR}/${RUN_ID}"
mkdir -p "${OUTDIR}"

"${REG_BIN}" \
  --step 1 \
  --bed "${BED_STEP1}" \
  --phenoFile "${PHENO_PATH}" --phenoCol "${PHENOTYPE}" \
  --covarFile "${COVAR_PATH}" \
  ${TRAIT_FLAG} \
  --bsize 1000 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --lowmem --lowmem-prefix "${TMPROOT}/${RUN_ID}" \
  --out "${OUTDIR}/${RUN_ID}_step1"
