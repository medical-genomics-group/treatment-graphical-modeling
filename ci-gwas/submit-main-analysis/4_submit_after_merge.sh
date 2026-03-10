#!/bin/bash
#SBATCH --job-name=ci_merged
#SBATCH --output=log/ci_merged_%A_%a.out
#SBATCH --error=log/ci_merged_%A_%a.err
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --array=0-100

set -euo pipefail
module load python/3.10.6

# ── Paths: fill these in before running ──────────────────────────────────────
CI_GWAS_ROOT="path/to/ci-gwas"
CIGWAS_BIN="path/to/ci-gwas.py"
# ─────────────────────────────────────────────────────────────────────────────

GEN_SETUPS=(
  "ever_cvd_with_drugs"
  "sbp_pre_post_1to60_rel_up_to_2nd_no_cvd"
  "dbp_pre_post_1to60_rel_up_to_2nd_no_cvd"
  "sbp_pre_post_1to60_rel_up_to_3rd_no_cvd"
  "dbp_pre_post_1to60_rel_up_to_3rd_no_cvd"
)

GEN_SETUP="${GEN_SETUPS[$SLURM_ARRAY_TASK_ID]}"

ALPHA_EXP=3
MAX_LEVEL=3
MAX_LEVEL_TWO=4
NUM_SAMPLES=450000
DEPTH=1

ALPHA=$(awk "BEGIN{print 10 ^ -$ALPHA_EXP}")

OUTROOT=${CI_GWAS_ROOT}/out/${GEN_SETUP}
RUNID=e${ALPHA_EXP}_l${MAX_LEVEL}_d${DEPTH}
CUSKDIR=${OUTROOT}/${RUNID}

MXM=${CUSKDIR}/merged_blocks.mxm.ld.bin
IXS=${CUSKDIR}/merged_blocks.ixs

MXPDIR=${CI_GWAS_ROOT}/mxp/main/${GEN_SETUP}
MXP=${MXPDIR}/merged_mxp_cors.tsv
MXPSE=${MXPDIR}/merged_mxp_sds.tsv

PXPDIR=${CI_GWAS_ROOT}/pxp/${GEN_SETUP}
PXP=${PXPDIR}/cor_pxp.tsv
PXPSE=${PXPDIR}/ses_pxp.tsv

TIMEINDEX=${CI_GWAS_ROOT}/time_ind/${GEN_SETUP}_time_index.txt

mkdir -p "${CUSKDIR}"

echo "HOST=$(hostname)"
echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "GEN_SETUP=${GEN_SETUP}"
echo "CUSKDIR=${CUSKDIR}"
echo "MXM=${MXM}"
echo "IXS=${IXS}"
echo "MXP=${MXP}"
echo "MXPSE=${MXPSE}"
echo "PXP=${PXP}"
echo "PXPSE=${PXPSE}"
echo "TIMEINDEX=${TIMEINDEX}"

srun ${CIGWAS_BIN} cuskss \
    --mxm "${MXM}" \
    --mxp "${MXP}" \
    --pxp "${PXP}" \
    --mxp-se "${MXPSE}" \
    --pxp-se "${PXPSE}" \
    --num-samples "${NUM_SAMPLES}" \
    --marker-indices "${IXS}" \
    --alpha "${ALPHA}" \
    --max-level-one "${MAX_LEVEL}" \
    --max-level-two "${MAX_LEVEL_TWO}" \
    --max-depth "${DEPTH}" \
    --time-index "${TIMEINDEX}" \
    --outdir "${CUSKDIR}"
