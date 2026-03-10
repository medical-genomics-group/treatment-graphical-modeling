#!/bin/bash
#SBATCH -a 1-100
#SBATCH --job-name=mxp%a
#SBATCH --output=./log/%A_%a.out
#SBATCH --error=./log/%A_%a.err
#SBATCH --time=48:00:00
#SBATCH --mem=64G

ID=$SLURM_ARRAY_TASK_ID
ARGS=$(sed -n "${ID}p" args-for-rerun.txt)
echo "[$ID] $ARGS"
read PHENODIR RUN CHR OUTDIR <<< "$ARGS"

module load python

mkdir -p "$OUTDIR"
srun python3 calc_mxp.py "$PHENODIR" "$RUN" "$CHR" "$OUTDIR"
