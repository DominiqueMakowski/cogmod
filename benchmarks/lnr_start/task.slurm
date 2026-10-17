#!/bin/bash
# One array task of the LNR start check: row SLURM_ARRAY_TASK_ID of
# cells.csv, compiled and sampled in its own scratch directory.
idx="${SLURM_ARRAY_TASK_ID:-1}"
work="${LS_SCRATCH}/work/${SLURM_ARRAY_JOB_ID:-0}_${idx}"
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S'): cell ${idx}"
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${LS_EB_MODULES}"
module load "${LS_R_MODULE}"
export R_LIBS="${LS_LIB}:${LS_PROD_LIB}${R_LIBS:+:$R_LIBS}"
rm -rf "$work" && mkdir -p "$work" || exit 1
export LS_WORK="$work"
Rscript "${LS_DIR}/fit_cell.R" "$idx"
status=$?
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
# Status 3 = a chain is incomplete (fit_cell.R): CmdStan's files are kept.
[ "$status" -eq 0 ] && rm -rf "$work"
exit $status
