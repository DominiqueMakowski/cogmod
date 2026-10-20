#!/usr/bin/env bash
# run.sh -- the LNR start check on Artemis (README.md). Run prep.R first.
#
#   ./run.sh push     push/ -> the cluster (one tar pipe)
#   ./run.sh smoke    the first cell, synchronously with srun
#   ./run.sh submit   the whole array, LS_MAX cells at a time (default 24)
#   ./run.sh queue    this user's queue
#   ./run.sh log      the tail of every task's output so far
#   ./run.sh pull     results back into results/, for summarise.R
#
# Each cell is one fit of 4 chains on 4 CPUs and compiles its own model on
# its own node. The library and the exclusions are the driftone prior check's
# (benchmarks/driftone_priors/run.sh): the stanli experiment's library (brms
# 2.23.1, as on the laptop) ahead of the production one, which is never
# written to; artemis-rtx-01 excluded. The GlobalProtect VPN must be up, and
# sshd rate-limits bursts of logins: one connection per command.
#
# --hint=nomultithread: Artemis counts hyperthreads as CPUs, so a bare
# --cpus-per-task=4 is 2 physical cores for 4 chains plus R (array
# 11596401). With whole cores a task takes 8 CPUs off the account's cap
# instead of 4. The chains lost in that array were not the cores' doing
# but a processx id collision between cells on one node (fit_cell.R).
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALIAS="${LS_HPC_ALIAS:-artemis}"
REMOTE="${LS_REMOTE:-/mnt/lustre/users/psych/dmm56/cogmod_lnrstart}"
SCRATCH="${LS_SCRATCH:-/mnt/lustre/scratch/psych/dmm56/cogmod_lnrstart}"
LIB="${LS_LIB:-/mnt/lustre/users/psych/dmm56/cogmod_stanli/R_libs}"
PROD_LIB="${LS_PROD_LIB:-/mnt/lustre/users/psych/dmm56/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${LS_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB="${LS_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
MAX="${LS_MAX:-24}"
SSH=ssh
[ -x /c/Windows/System32/OpenSSH/ssh.exe ] && SSH=/c/Windows/System32/OpenSSH/ssh.exe
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@"; }
EXPORTS="LS_DIR=${REMOTE},LS_SCRATCH=${SCRATCH},LS_EB_MODULES=${EB},LS_R_MODULE=${R_MODULE},LS_LIB=${LIB},LS_PROD_LIB=${PROD_LIB}"
n_cells() { tail -n +2 "$HERE/push/cells.csv" | wc -l | tr -d ' '; }

case "${1:-}" in
  push)
    for f in "$HERE"/push/*.R "$HERE"/push/*.slurm; do sed -i 's/\r$//' "$f"; done
    tar -C "$HERE/push" -cf - . | remote "mkdir -p $REMOTE $SCRATCH && rm -rf $REMOTE/tree && tar -C $REMOTE -xf - && chmod +x $REMOTE/task.slurm" ;;
  smoke)
    remote "cd $REMOTE && srun --partition=short --time=01:00:00 --cpus-per-task=4 --hint=nomultithread --mem=8G --exclude=artemis-rtx-01 bash -c 'export ${EXPORTS//,/ } SLURM_ARRAY_TASK_ID=1 SLURM_ARRAY_JOB_ID=smoke; bash task.slurm'" ;;
  submit)
    remote "cd $REMOTE && sbatch --array=1-$(n_cells)%${MAX} --partition=general --cpus-per-task=4 --hint=nomultithread --mem=8G --time=03:00:00 --exclude=artemis-rtx-01 --job-name=cogmod_lnrstart --output=$SCRATCH/ls_%A_%a.out --error=$SCRATCH/ls_%A_%a.err --chdir=$REMOTE --export=ALL,$EXPORTS task.slurm" ;;
  queue)
    remote "squeue -u \$USER -o '%.10i %.9P %.16j %.2t %.10M %.6D %R'" ;;
  log)
    remote "cd $SCRATCH && for f in ls_*.out; do echo \"== \$f\"; tail -n 2 \$f; done" ;;
  pull)
    mkdir -p "$HERE/results"
    remote "cd $REMOTE/results && tar -cf - ." | tar -C "$HERE/results" -xf - ;;
  *) sed -n '2,9p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
