#!/usr/bin/env bash
#SBATCH --job-name=rg_step2_chr
#SBATCH --array=66-6512
#SBATCH --time=3-00:00:00
#SBATCH --cpus-per-task=8
#SBATCH --mem=32G
#SBATCH --no-requeue
#SBATCH --output=logs/step2_chr_%A_%a.out
#SBATCH --error=logs/step2_chr_%A_%a.err

set -euo pipefail
mkdir -p logs

CONDA_ENV="regenie412"
BED_STEP2="<BED_STEP2>"

OUTROOT="<OUTROOT>"
STEP1DIR="${OUTROOT}/step1"
STEP2DIR="${OUTROOT}/step2"
PAIRS="${OUTROOT}/pairs_all.tsv"

mkdir -p "${STEP2DIR}"

module load conda
conda activate "${CONDA_ENV}"
REG_BIN="$(which regenie)"
export OMP_NUM_THREADS="${SLURM_CPUS_PER_TASK}"

IDX0=$((SLURM_ARRAY_TASK_ID - 1))
PAIR_IDX=$((IDX0 / 22 + 1))
CHR=$((IDX0 % 22 + 1))

LINE="$(sed -n "$((PAIR_IDX + 1))p" "${PAIRS}")"

IFS=$'\t' read -r RUN_ID SCOPE TRAIT PHENOTYPE COVARS_IN_MODEL PHENO_PATH COVAR_PATH INTERACTION_VAR <<< "${LINE}"

TRAIT_FLAG="--qt"
FIRTH_OPTS=""
if [[ "${TRAIT}" == "CVD" || "${TRAIT}" =~ ^CLASS_ ]]; then
  TRAIT_FLAG="--bt"
  FIRTH_OPTS="--firth --approx --pThresh 0.000001"
fi

PRED_LIST="${STEP1DIR}/${RUN_ID}/${RUN_ID}_step1_pred.list"

OUTDIR="${STEP2DIR}/${RUN_ID}/chr${CHR}"
mkdir -p "${OUTDIR}"

"${REG_BIN}" \
  --step 2 \
  --bed "${BED_STEP2}" \
  --chr "${CHR}" \
  --phenoFile "${PHENO_PATH}" --phenoCol "${PHENOTYPE}" \
  --covarFile "${COVAR_PATH}" \
  ${TRAIT_FLAG} \
  ${FIRTH_OPTS} \
  --pred "${PRED_LIST}" \
  --bsize 1000 \
  --threads "${SLURM_CPUS_PER_TASK}" \
  --out "${OUTDIR}/${RUN_ID}_chr${CHR}_step2"
