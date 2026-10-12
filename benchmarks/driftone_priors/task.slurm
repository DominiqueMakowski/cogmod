#!/bin/bash
# One array task of the cogmod prior check: row SLURM_ARRAY_TASK_ID of
# cells.csv, compiled and sampled in its own scratch directory.
idx="${SLURM_ARRAY_TASK_ID:-1}"
work="${PC_SCRATCH}/work/${SLURM_ARRAY_JOB_ID:-0}_${idx}"
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S'): cell ${idx}"
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${PC_EB_MODULES}"
module load "${PC_R_MODULE}"
export R_LIBS="${PC_LIB}:${PC_PROD_LIB}${R_LIBS:+:$R_LIBS}"
rm -rf "$work" && mkdir -p "$work" || exit 1
export PC_WORK="$work"
Rscript "${PC_DIR}/refit_cell.R" "$idx"
status=$?
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
rm -rf "$work"
exit $status
