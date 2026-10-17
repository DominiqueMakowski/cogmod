#!/usr/bin/env bash
# run.sh -- the driftone prior refits on Artemis (README.md). Run prep.R first.
#
#   ./run.sh push     push/ -> the cluster (one tar pipe)
#   ./run.sh smoke    the last cell (bisa, arm B), synchronously with srun
#   ./run.sh submit   the whole array, PC_MAX cells at a time (default 32)
#   ./run.sh pull     results back into results/cells/, for collect.R
#
# Each cell is one fit of 4 chains on 4 CPUs and compiles its own model on
# its own node, so nothing compiled is shared across node types. The library
# is the stanli experiment's (brms 2.23.1, as on the laptop, so both platforms
# sample the same Stan program) ahead of the production one, which is never
# written to. artemis-rtx-01 is excluded: it ran CmdStan 5x slow on
# 2026-10-12. The GlobalProtect VPN must be up, and sshd rate-limits bursts of
# logins: one connection per command.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ALIAS="${PC_HPC_ALIAS:-artemis}"
REMOTE="${PC_REMOTE:-/mnt/lustre/users/psych/dmm56/cogmod_priorcheck}"
SCRATCH="${PC_SCRATCH:-/mnt/lustre/scratch/psych/dmm56/cogmod_priorcheck}"
LIB="${PC_LIB:-/mnt/lustre/users/psych/dmm56/cogmod_stanli/R_libs}"
PROD_LIB="${PC_PROD_LIB:-/mnt/lustre/users/psych/dmm56/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${PC_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB="${PC_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
MAX="${PC_MAX:-32}"
SSH=ssh
[ -x /c/Windows/System32/OpenSSH/ssh.exe ] && SSH=/c/Windows/System32/OpenSSH/ssh.exe
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@"; }
EXPORTS="PC_DIR=${REMOTE},PC_SCRATCH=${SCRATCH},PC_EB_MODULES=${EB},PC_R_MODULE=${R_MODULE},PC_LIB=${LIB},PC_PROD_LIB=${PROD_LIB}"
n_cells() { tail -n +2 "$HERE/push/cells.csv" | wc -l | tr -d ' '; }

case "${1:-}" in
  push)
    for f in "$HERE"/push/*.R "$HERE"/push/*.slurm; do sed -i 's/\r$//' "$f"; done
    tar -C "$HERE/push" -cf - . | remote "mkdir -p $REMOTE $SCRATCH && rm -rf $REMOTE/tree && tar -C $REMOTE -xf - && chmod +x $REMOTE/task.slurm" ;;
  smoke)
    remote "cd $REMOTE && srun --partition=short --time=00:40:00 --cpus-per-task=4 --mem=8G --exclude=artemis-rtx-01 bash -c 'export ${EXPORTS//,/ } SLURM_ARRAY_TASK_ID=$(n_cells) SLURM_ARRAY_JOB_ID=smoke; bash task.slurm'" ;;
  submit)
    remote "cd $REMOTE && sbatch --array=1-$(n_cells)%${MAX} --partition=general --cpus-per-task=4 --mem=8G --time=04:00:00 --exclude=artemis-rtx-01 --job-name=cogmod_priorcheck --output=$SCRATCH/pc_%A_%a.out --error=$SCRATCH/pc_%A_%a.err --chdir=$REMOTE --export=ALL,$EXPORTS task.slurm" ;;
  pull)
    mkdir -p "$HERE/results/cells"
    remote "cd $REMOTE/results && tar -cf - ." | tar -C "$HERE/results/cells" -xf - ;;
  *) sed -n '2,7p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
