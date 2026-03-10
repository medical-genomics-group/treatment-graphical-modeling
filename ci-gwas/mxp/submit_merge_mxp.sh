#!/bin/bash
#SBATCH --job-name=merge_mxp
#SBATCH --output=./mxp_merge_%A_%a.out
#SBATCH --error=./mxp_merge_%A_%a.err
#SBATCH --time=24:00:00
#SBATCH --mem=64G
#SBATCH --cpus-per-task=1
#SBATCH --array=0-100

module load python

# ── Paths: fill these in before running ──────────────────────────────────────
MXP_ROOT="path/to/mxp/pruned"
BIM="path/to/bedfiles/autosomes.bim"
# ─────────────────────────────────────────────────────────────────────────────

srun python merge_mxp.py \
  --mxp-root "$MXP_ROOT" \
  --bim "$BIM" \
  --only dbp_pre_post_1to60_no_cvd \
  --only sbp_pre_post_1to60_no_cvd \
  --only dbp_pre_post_1to60_age5_with_statins_no_cvd \
  --only sbp_pre_post_1to60_age5_with_statins_no_cvd