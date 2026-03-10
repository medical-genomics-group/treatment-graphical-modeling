#!/bin/bash
#SBATCH --job-name=merge_blocks
#SBATCH --output=./log/merge_blocks_%A_%a.out
#SBATCH --error=./log/merge_blocks_%A_%a.err
#SBATCH --time=00:20:00
#SBATCH --mem=8G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
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
DEPTH=1

OUTROOT=${CI_GWAS_ROOT}/out/${GEN_SETUP}
RUNID=e${ALPHA_EXP}_l${MAX_LEVEL}_d${DEPTH}
OUTDIR=${OUTROOT}/${RUNID}

BLOCK_FILE=${CI_GWAS_ROOT}/submit/blockfiles/autosomes_m11000.blocks

echo "TASK=${SLURM_ARRAY_TASK_ID} GEN_SETUP=${GEN_SETUP}"
echo "OUTDIR=${OUTDIR}"
echo "BLOCK_FILE=${BLOCK_FILE}"
echo "HOST=$(hostname)"

srun ${CIGWAS_BIN} \
    merge-block-outputs "${OUTDIR}" "${BLOCK_FILE}"