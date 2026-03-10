#!/bin/bash
#SBATCH --job-name=merge_blocks
#SBATCH --output=log/mxm_merge_%A_%a.out
#SBATCH --error=log/mxm_merge_%A_%a.err
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
CUSK_DIR=${OUTROOT}/${RUNID}
PXP=${CI_GWAS_ROOT}/pxp/${GEN_SETUP}/cor_pxp.tsv

echo "CUSK_DIR=${CUSK_DIR}"
echo "PXP=${PXP}"
echo "HOST=$(hostname)"

srun calc_mxm_merged.py "${CUSK_DIR}" "${PXP}"
