#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run.sh -- drive the inits ablation on Artemis (Sussex HPC) over SSH.
#
#   ./run.sh build            build the source tarball of this tree (local)
#   ./run.sh push             tarball + scripts + data -> the cluster (one tar pipe)
#   ./run.sh install          install the tarball into the ablation's own R library
#   ./run.sh submit A|B [..]  submit one array job for an experiment (extra args -> sbatch)
#   ./run.sh queue            squeue for this user
#   ./run.sh progress         how many cells of each experiment have written a row
#   ./run.sh log [pattern]    tail the newest .out/.err (optionally matching)
#   ./run.sh pull             copy the result rows back into results/
#   ./run.sh cancel <jobid>   scancel
#   ./run.sh sh '<cmd>'       run a command on the login node
#
# Conventions follow analysis/server/hpc of the Illusion Game project (SSH
# alias `artemis`, the CmdStanR module, one SSH connection per command
# because sshd rate-limits bursts). The GlobalProtect VPN must be up.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/../.." && pwd)"
# Windows R wants Windows paths; Git Bash gives POSIX ones.
HERE_W="$(cygpath -m "$HERE" 2>/dev/null || echo "$HERE")"
PKG_W="$(cygpath -m "$PKG" 2>/dev/null || echo "$PKG")"

ALIAS="${ABL_HPC_ALIAS:-artemis}"
USER_="${ABL_HPC_USER:-dmm56}"
GROUP_="${ABL_HPC_GROUP:-psych}"
REMOTE="${ABL_REMOTE:-/mnt/lustre/users/${GROUP_}/${USER_}/cogmod_inits_ablation}"
SCRATCH="${ABL_SCRATCH:-/mnt/lustre/scratch/${GROUP_}/${USER_}/cogmod_inits_ablation}"
LIB="${REMOTE}/R_libs"
PROD_LIB="${ABL_PROD_LIB:-/mnt/lustre/users/${GROUP_}/${USER_}/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${ABL_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB_MODULES="${ABL_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
RSCRIPT_LOCAL="${RSCRIPT:-/c/Program Files/R/R-4.5.3/bin/Rscript.exe}"
R_LOCAL="${R_EXE:-/c/Program Files/R/R-4.5.3/bin/R.exe}"

if [ -x /c/Windows/System32/OpenSSH/ssh.exe ]; then
  SSH=/c/Windows/System32/OpenSSH/ssh.exe; SCP=/c/Windows/System32/OpenSSH/scp.exe
else
  SSH=ssh; SCP=scp
fi
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
info() { printf '\033[36m==>\033[0m %s\n' "$*"; }
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@" 2> >(grep -v 'module: command not found' >&2); }

# Number of cells per experiment, read from cells.R so the array is never
# wider or narrower than the design.
n_cells() {
  "$RSCRIPT_LOCAL" -e "source('$HERE_W/cells.R'); cat(nrow(ablation_cells('$1')))"
}

cmd_build() {
  info "exporting speed_acc (rtdists) to data/speed_acc.csv"
  mkdir -p "$HERE/data"
  "$RSCRIPT_LOCAL" -e "data(speed_acc, package = 'rtdists'); d <- data.frame(id = as.integer(as.character(speed_acc\$id)), condition = as.character(speed_acc\$condition), rt = speed_acc\$rt); write.csv(d, '$HERE_W/data/speed_acc.csv', row.names = FALSE); cat(nrow(d), 'rows,', sum(d\$rt <= 2), 'at or under 2 s\n')"
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
  mkdir -p "$tmp/data" && cp "$HERE"/data/*.csv "$tmp/data/"
  tar -C "$tmp" -cf - . | remote "mkdir -p '${REMOTE}' '${SCRATCH}' && tar -C '${REMOTE}' -xf - && chmod +x '${REMOTE}'/*.slurm && ls -la '${REMOTE}' '${REMOTE}/data'"
}

cmd_install() {
  info "installing this tree's cogmod into ${LIB} (compute node, a few minutes)"
  remote "cd '${REMOTE}' && srun --partition=short --time=00:40:00 --cpus-per-task=2 --mem=8G --job-name=cogmod_inits_install bash -c '
    . /etc/profile.d/lmod.sh; module use \"${EB_MODULES}\"; module load \"${R_MODULE}\"
    export ABL_LIB=\"${LIB}\"; export R_LIBS=\"${LIB}:${PROD_LIB}\"
    Rscript install.R'"
}

cmd_submit() {
  local exp="${1:?usage: run.sh submit A|B [sbatch args]}"; shift || true
  local n; n="$(n_cells "$exp")"
  local cpus mem
  case "$exp" in
    A) cpus=8;  mem=16G ;;
    B) cpus=16; mem=32G ;;
    *) die "unknown experiment $exp" ;;
  esac
  info "submitting experiment ${exp}: ${n} cells, ${cpus} CPUs each"
  remote "mkdir -p '${SCRATCH}' '${REMOTE}/results' && cd '${REMOTE}' && sbatch \
    --array=1-${n} --cpus-per-task=${cpus} --mem=${mem} \
    --job-name='cogmod_inits_${exp}' \
    --output='${SCRATCH}/inits_${exp}_%A_%a.out' --error='${SCRATCH}/inits_${exp}_%A_%a.err' \
    --chdir='${REMOTE}' \
    --export=ALL,ABL_EXPERIMENT='${exp}',ABL_DIR='${REMOTE}',ABL_OUT='${REMOTE}/results',ABL_DATA='${REMOTE}/data',ABL_LIB='${LIB}',ABL_PROD_LIB='${PROD_LIB}',ABL_R_MODULE='${R_MODULE}',ABL_EB_MODULES='${EB_MODULES}' \
    $* fit.slurm"
}

cmd_queue() { remote "squeue -u '${USER_}' -o '%.12i %.12P %.20j %.8T %.10M %.10L %.5C %R' | grep -E 'JOBID|cogmod_inits' || echo '(no ablation jobs queued)'"; }

cmd_progress() {
  remote "cd '${REMOTE}/results' 2>/dev/null || { echo 'no results yet'; exit 0; }
    for e in A B; do n=\$(ls \${e}_*.csv 2>/dev/null | wc -l); echo \"\$e: \$n rows\"; done
    echo '--- errors in rows:'; grep -l 'error:' *.csv 2>/dev/null | head; echo '--- failed tasks (.err with Error):'; grep -l '^Error\|Execution halted' '${SCRATCH}'/*.err 2>/dev/null | wc -l"
}

cmd_log() {
  local key="${1:-}"
  remote "cd '${SCRATCH}' 2>/dev/null || exit 1; for f in \$(ls -t *${key}*.out *${key}*.err 2>/dev/null | head -4); do echo \"===== \$f =====\"; tail -n 25 \"\$f\"; done"
}

cmd_pull() {
  mkdir -p "$HERE/results"
  info "pulling result rows into ${HERE}/results"
  remote "cd '${REMOTE}/results' && tar -cf - *.csv" | tar -C "$HERE/results" -xf -
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
  *) sed -n '2,20p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
