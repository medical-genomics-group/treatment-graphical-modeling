#!/bin/bash
#SBATCH --job-name=ci_bp_%a
#SBATCH --output=./log/ci_%J_%a.out
#SBATCH --error=./log/ci_%J_%a.err
#SBATCH --time=00:30:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=1-100

set -euo pipefail
module load python

# ── Paths: fill these in before running ──────────────────────────────────────
CI_GWAS_ROOT="path/to/ci-gwas"
MXM_ROOT="path/to/mxm"          # directory containing per-chromosome LD binary files
DATA_ROOT="path/to/bedfiles"    # directory with PLINK bedfiles
CIGWAS_BIN="path/to/ci-gwas.py"
# ─────────────────────────────────────────────────────────────────────────────

ALPHA_EXP=3
MAX_LEVEL=3
MAX_LEVEL_TWO=14
NUM_SAMPLES=450000
DEPTH=1

ID=${SLURM_ARRAY_TASK_ID}


ARGS=$(sed -n "${ID}p" args.txt)
read -r CIX BLOCK_IX GEN_SETUP <<< "${ARGS}"

ZERO_BASED_BLOCK_IX=$((BLOCK_IX-1))
ALPHA=$(awk "BEGIN{print 10 ^ -$ALPHA_EXP}")

EXPDIR=${CI_GWAS_ROOT}/mxp/main/${GEN_SETUP}
PXPDIR=${CI_GWAS_ROOT}/pxp/${GEN_SETUP}
OUTROOT=${CI_GWAS_ROOT}/out/${GEN_SETUP}
RUNID=e${ALPHA_EXP}_l${MAX_LEVEL}_d${DEPTH}
OUTDIR=${OUTROOT}/${RUNID}

mkdir -p "${OUTDIR}"

BLOCK_FILE=${CI_GWAS_ROOT}/submit/blockfiles/c${CIX}_m11000.blocks
MXM=${MXM_ROOT}/c${CIX}_b${BLOCK_IX}.ld.bin
MXP=${EXPDIR}/c${CIX}_cors.tsv
MXPSE=${EXPDIR}/c${CIX}_sds.tsv
PXP=${PXPDIR}/cor_pxp.tsv
PXPSE=${PXPDIR}/ses_pxp.tsv
TIMEINDEX=${CI_GWAS_ROOT}/time_ind/${GEN_SETUP}_time_index.txt

echo "HOST=$(hostname)"
echo "TASK=${ID}  CIX=${CIX}  BLOCK_IX=${BLOCK_IX}  ZB=${ZERO_BASED_BLOCK_IX}  GEN_SETUP=${GEN_SETUP}"
echo "MXM=${MXM}"
echo "BLOCK_FILE=${BLOCK_FILE}"
echo "OUTDIR=${OUTDIR}"

python preflight_check.py \
  --gen-setup "${GEN_SETUP}" \
  --pxp "${PXP}" \
  --mxp-dir "${EXPDIR}" \
  --time-index "${TIMEINDEX}" \
  --out-time-table "${OUTDIR}/preflight_table.tsv"

srun ${CIGWAS_BIN} cuskss \
  --mxm "${MXM}" \
  --mxp "${MXP}" \
  --pxp "${PXP}" \
  --mxp-se "${MXPSE}" \
  --pxp-se "${PXPSE}" \
  --num-samples "${NUM_SAMPLES}" \
  --block-index "${ZERO_BASED_BLOCK_IX}" \
  --blockfile "${BLOCK_FILE}" \
  --alpha "${ALPHA}" \
  --max-level-one "${MAX_LEVEL}" \
  --max-level-two "${MAX_LEVEL_TWO}" \
  --max-depth "${DEPTH}" \
  --outdir "${OUTDIR}" \
  --time-index "${TIMEINDEX}"

