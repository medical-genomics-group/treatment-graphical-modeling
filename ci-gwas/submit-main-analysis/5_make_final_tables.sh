#!/bin/bash
#SBATCH --job-name=final_table
#SBATCH --output=log/final_table_%A_%a.out
#SBATCH --error=log/final_table_%A_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=10G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --array=0-100

set -euo pipefail
module load python/3.10.6

mkdir -p log

# ── Paths: fill these in before running ──────────────────────────────────────
CI_GWAS_ROOT="path/to/ci-gwas"
# ─────────────────────────────────────────────────────────────────────────────

GEN_SETUPS=(
    "sbp_pre_post_1to60_no_cvd"
    "dbp_pre_post_1to60_no_cvd"
    "bp_pre_post_1to60_no_cvd"
    "sbp_pre_post_1to60_age5_with_statins_no_cvd"
    "dbp_pre_post_1to60_age5_with_statins_no_cvd"
    "sbp_pre_post_1to60_age5_no_statins_no_cvd"
    "dbp_pre_post_1to60_age5_no_statins_no_cvd"
    "sbp_pre_post_1to60_no_statins_no_cvd"
    "dbp_pre_post_1to60_no_statins_no_cvd"
    "sbp_pre_post_1to60_age4_with_statins_no_cvd"
    "dbp_pre_post_1to60_age4_with_statins_no_cvd"
    "sbp_pre_post_1to60_age5_shifted_with_statins_no_cvd"
    "dbp_pre_post_1to60_age5_shifted_with_statins_no_cvd"
    "sbp_pre_post_1to30_no_cvd"
    "dbp_pre_post_1to30_no_cvd"
    "ldl_pre_post_1to30_no_cvd"
    "ldl_pre_post_1to60_no_cvd"
    "ldl_pre_post_1to60_age5_no_cvd"
    "ldl_pre_post_1to60_statins_only_no_cvd"
    "sbp_pre_post_1to60_with_prior_cvd"
    "dbp_pre_post_1to60_with_prior_cvd"
    "sbp_pre_post_1to60_age5_with_prior_cvd"
    "dbp_pre_post_1to60_age5_with_prior_cvd"
    "sbp_pre_post_1to60_with_ever_cvd"
    "dbp_pre_post_1to60_with_ever_cvd"
    "prior_cvd"
    "ever_cvd"
    "ever_cvd_with_drugs"
    "sbp_pre_post_1to60_with_copying_no_cvd"
    "dbp_pre_post_1to60_with_copying_no_cvd"
    "sbp_pre_post_1to60_on_drugs_no_cvd"
    "dbp_pre_post_1to60_on_drugs_no_cvd"
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

echo "SLURM_ARRAY_TASK_ID=${SLURM_ARRAY_TASK_ID}"
echo "GEN_SETUP=${GEN_SETUP}"
echo "CUSK_DIR=${CUSK_DIR}"
echo "PXP=${PXP}"
echo "HOST=$(hostname)"

srun create_table_cuskss.py "${CUSK_DIR}" "${PXP}"
