#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run.sh -- drive the LBA screen on Artemis (Sussex HPC) over SSH.
#
#   ./run.sh build              build the source tarball of this tree (local)
#   ./run.sh push               tarball + scripts -> the cluster (one tar pipe)
#   ./run.sh install            install the tarball into the screen's own R library
#   ./run.sh submit lba|rdm [..]  one array job, one task per chain (extra args -> sbatch)
#   ./run.sh queue              squeue for this user
#   ./run.sh progress           how many chains of each experiment have written a row
#   ./run.sh log [pattern]      tail the newest .out/.err (optionally matching)
#   ./run.sh pull               copy the results back into results/
#   ./run.sh cancel <jobid>     scancel
#   ./run.sh sh '<cmd>'         run a command on the login node
#
# Copied from benchmarks/inits_ablation/run.sh, with its own directories so
# that neither the production fits nor the ablation's results are touched.
# The GlobalProtect VPN must be up.
#
# The array is throttled (%SCREEN_MAX concurrent tasks, default 32, so 128
# CPUs): the account has a 550-CPU cap across every partition, and the
# September ablation's unthrottled array held all of it.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/../.." && pwd)"
HERE_W="$(cygpath -m "$HERE" 2>/dev/null || echo "$HERE")"
PKG_W="$(cygpath -m "$PKG" 2>/dev/null || echo "$PKG")"

ALIAS="${SCREEN_HPC_ALIAS:-artemis}"
USER_="${SCREEN_HPC_USER:-dmm56}"
GROUP_="${SCREEN_HPC_GROUP:-psych}"
REMOTE="${SCREEN_REMOTE:-/mnt/lustre/users/${GROUP_}/${USER_}/cogmod_lba_screen}"
SCRATCH="${SCREEN_SCRATCH:-/mnt/lustre/scratch/${GROUP_}/${USER_}/cogmod_lba_screen}"
LIB="${REMOTE}/R_libs"
PROD_LIB="${SCREEN_PROD_LIB:-/mnt/lustre/users/${GROUP_}/${USER_}/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${SCREEN_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB_MODULES="${SCREEN_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
MAX="${SCREEN_MAX:-32}"
RSCRIPT_LOCAL="${RSCRIPT:-/c/Program Files/R/R-4.5.3/bin/Rscript.exe}"
R_LOCAL="${R_EXE:-/c/Program Files/R/R-4.5.3/bin/R.exe}"

if [ -x /c/Windows/System32/OpenSSH/ssh.exe ]; then
  SSH=/c/Windows/System32/OpenSSH/ssh.exe
else
  SSH=ssh
fi
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
info() { printf '\033[36m==>\033[0m %s\n' "$*"; }
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@" 2> >(grep -v 'module: command not found' >&2); }

n_cells() {
  "$RSCRIPT_LOCAL" -e "setwd('$HERE_W'); source('cells.R'); cat(nrow(screen_cells('$1')))"
}

cmd_build() {
  info "building the source tarball of $PKG"
  rm -f "$HERE"/cogmod_*.tar.gz
  (cd "$HERE" && "$R_LOCAL" CMD build --no-build-vignettes --no-manual "$PKG_W" 2>&1 | tail -5)
  ls -la "$HERE"/cogmod_*.tar.gz
}

cmd_push() {
  info "pushing to ${REMOTE}"
  local tmp; tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN
  for f in "$HERE"/*.R "$HERE"/*.slurm; do sed 's/\r$//' "$f" > "$tmp/$(basename "$f")"; done
  cp "$HERE"/cogmod_*.tar.gz "$tmp/"
  tar -C "$tmp" -cf - . | remote "mkdir -p '${REMOTE}' '${SCRATCH}' && tar -C '${REMOTE}' -xf - && chmod +x '${REMOTE}'/*.slurm && ls -la '${REMOTE}'"
}

cmd_install() {
  info "installing this tree's cogmod into ${LIB} (compute node, a few minutes)"
  remote "cd '${REMOTE}' && srun --partition=short --time=00:40:00 --cpus-per-task=2 --mem=8G --job-name=cogmod_screen_install bash -c '
    . /etc/profile.d/lmod.sh; module use \"${EB_MODULES}\"; module load \"${R_MODULE}\"
    export SCREEN_LIB=\"${LIB}\"; export R_LIBS=\"${LIB}:${PROD_LIB}\"
    Rscript install.R'"
}

cmd_submit() {
  local exp="${1:?usage: run.sh submit lba|rdm [sbatch args]}"; shift || true
  local n; n="$(n_cells "$exp")"
  info "submitting ${exp}: ${n} chains, 4 CPUs each, at most ${MAX} at once"
  remote "mkdir -p '${SCRATCH}' '${REMOTE}/results' && cd '${REMOTE}' && sbatch \
    --array=1-${n}%${MAX} --cpus-per-task=4 --mem=16G \
    --job-name='cogmod_lba_screen_${exp}' \
    --output='${SCRATCH}/screen_${exp}_%A_%a.out' --error='${SCRATCH}/screen_${exp}_%A_%a.err' \
    --chdir='${REMOTE}' \
    --export=ALL,SCREEN_EXPERIMENT='${exp}',SCREEN_DIR='${REMOTE}',SCREEN_OUT='${REMOTE}/results',SCREEN_LIB='${LIB}',SCREEN_PROD_LIB='${PROD_LIB}',SCREEN_R_MODULE='${R_MODULE}',SCREEN_EB_MODULES='${EB_MODULES}' \
    $* fit.slurm"
}

cmd_queue() { remote "squeue -u '${USER_}' -o '%.12i %.12P %.24j %.8T %.10M %.10L %.5C %R' | grep -E 'JOBID|cogmod_lba_screen' || echo '(no screen jobs queued)'"; }

cmd_progress() {
  remote "cd '${REMOTE}/results' 2>/dev/null || { echo 'no results yet'; exit 0; }
    for e in lba rdm; do n=\$(ls \${e}_*.csv 2>/dev/null | wc -l); echo \"\$e: \$n chains\"; done
    echo '--- errors in rows:'; grep -l 'error:' *.csv 2>/dev/null | head; echo '--- failed tasks (.err with Error):'; grep -l '^Error\|Execution halted' '${SCRATCH}'/*.err 2>/dev/null | wc -l"
}

cmd_log() {
  local key="${1:-}"
  remote "cd '${SCRATCH}' 2>/dev/null || exit 1; for f in \$(ls -t *${key}*.out *${key}*.err 2>/dev/null | head -4); do echo \"===== \$f =====\"; tail -n 25 \"\$f\"; done"
}

cmd_pull() {
  mkdir -p "$HERE/results"
  info "pulling results into ${HERE}/results"
  remote "cd '${REMOTE}/results' && tar -cf - *.csv *.rds" | tar -C "$HERE/results" -xf -
  ls "$HERE/results" | wc -l
}

cmd_cancel() { remote "scancel '${1:?jobid}' && echo cancelled"; }
cmd_sh() { remote "${1:?cmd}"; }

case "${1:-}" in
  build)    shift; cmd_build "$@" ;;
  push)     shift; cmd_push "$@" ;;
  install)  shift; cmd_install "$@" ;;
  submit)   shift; cmd_submit "$@" ;;
  queue)    shift; cmd_queue "$@" ;;
  progress) shift; cmd_progress "$@" ;;
  log)      shift; cmd_log "$@" ;;
  pull)     shift; cmd_pull "$@" ;;
  cancel)   shift; cmd_cancel "$@" ;;
  sh)       shift; cmd_sh "$@" ;;
  *) sed -n '2,22p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
