#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run.sh -- run the stanli benchmark (benchmarks/stanli) on Artemis over SSH.
#
#   ./run.sh push               the scripts, the minimal package tree and the
#                               stanli tarball -> the cluster (one tar pipe)
#   ./run.sh install            brms 2.23.1 (GitHub) + stanli $STANLI_VERSION (default
#                               0.19.1; + runtime) into the experiment's own R library,
#                               on a compute node. STANLI_TAG=0.19.1 on submit/pull keeps
#                               a rerun's results apart (results/<kind>_<tag>_<index>/)
#   ./run.sh submit fit [n]     n seeds of bench.R's 2x2 fit, 4 CPUs each (default 16)
#   ./run.sh submit fitpar [n]  n seeds with stanli's chains in parallel (one process
#                               per chain; stanli's threads from random inits), 4 CPUs
#   ./run.sh submit grad [n]    n copies of grad + bisect + repro, 2 CPUs each (default 4)
#                               n may be a range ("10-12") to rerun those indices;
#                               extra args go to sbatch
#
# artemis-general-02 (a VM, CPU "AMD EPYC-Genoa Processor") killed stanli with
# SIGILL (exit 132) inside stanli_model() on 2026-10-03, every time, while the
# EPYC 9334/9355 nodes ran it; reruns there want --exclude=artemis-general-02.
#   ./run.sh queue              squeue for this user
#   ./run.sh progress           how many tasks of each kind have finished
#   ./run.sh log [pattern]      tail the newest .out/.err (optionally matching)
#   ./run.sh pull               copy results/ back into ../results/hpc/
#   ./run.sh cancel <jobid>     scancel
#   ./run.sh sh '<cmd>'         run a command on the login node
#
# Copied from benchmarks/lba_screen/run.sh, with its own directories, so that
# neither the production fits nor the other experiments are touched. brms is
# installed at the local version (2.23.1 from GitHub, e71e9d7, against the
# cluster's 2.21.0) so that both platforms sample the same Stan program. The GlobalProtect VPN must
# be up. Arrays are throttled (%STANLI_MAX, default 16): the account has a
# 550-CPU cap across every partition.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/../../.." && pwd)"

ALIAS="${STANLI_HPC_ALIAS:-artemis}"
USER_="${STANLI_HPC_USER:-dmm56}"
GROUP_="${STANLI_HPC_GROUP:-psych}"
REMOTE="${STANLI_REMOTE:-/mnt/lustre/users/${GROUP_}/${USER_}/cogmod_stanli}"
SCRATCH="${STANLI_SCRATCH:-/mnt/lustre/scratch/${GROUP_}/${USER_}/cogmod_stanli}"
LIB="${REMOTE}/R_libs"
PROD_LIB="${STANLI_PROD_LIB:-/mnt/lustre/users/${GROUP_}/${USER_}/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${STANLI_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB_MODULES="${STANLI_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
MAX="${STANLI_MAX:-16}"
STANLI_VERSION="${STANLI_VERSION:-0.19.1}"
# Tasks write to results/<kind>[_<tag>]_<index>/ and bench.R grad writes its
# files with --tag, so a rerun on a new stanli keeps the old results.
TAG="${STANLI_TAG:-}"

if [ -x /c/Windows/System32/OpenSSH/ssh.exe ]; then
  SSH=/c/Windows/System32/OpenSSH/ssh.exe
else
  SSH=ssh
fi
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
info() { printf '\033[36m==>\033[0m %s\n' "$*"; }
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@" 2> >(grep -v 'module: command not found' >&2); }
EXPORTS="STANLI_DIR='${REMOTE}',STANLI_SCRATCH='${SCRATCH}',STANLI_LIB='${LIB}',STANLI_PROD_LIB='${PROD_LIB}',STANLI_R_MODULE='${R_MODULE}',STANLI_EB_MODULES='${EB_MODULES}',STANLI_TAG='${TAG}'"

# What load_all() and the benchmark need, nothing else: the package's R code
# and metadata, gradient_programs.R, and benchmarks/stanli's scripts and Stan
# files (not its results, which each task writes afresh).
cmd_push() {
  info "pushing to ${REMOTE}"
  local tmp; tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN
  mkdir -p "$tmp/tree/benchmarks/stanli"
  (cd "$PKG" && cp -r DESCRIPTION NAMESPACE R data inst "$tmp/tree/")
  sed 's/\r$//' "$PKG/benchmarks/gradient_programs.R" > "$tmp/tree/benchmarks/gradient_programs.R"
  for f in "$PKG"/benchmarks/stanli/*.R "$PKG"/benchmarks/stanli/*.stan; do
    sed 's/\r$//' "$f" > "$tmp/tree/benchmarks/stanli/$(basename "$f")"
  done
  for f in "$HERE"/*.R "$HERE"/*.slurm; do sed 's/\r$//' "$f" > "$tmp/$(basename "$f")"; done
  local tgz="stanli_${STANLI_VERSION}.tar.gz"
  curl -sSfL -o "$tmp/$tgz" "https://github.com/seantalts/stanli/releases/download/v${STANLI_VERSION}/${tgz}"
  # install.R takes the one tarball it finds, so the previous version's goes.
  tar -C "$tmp" -cf - . | remote "mkdir -p '${REMOTE}' '${SCRATCH}' && rm -rf '${REMOTE}/tree' '${REMOTE}'/stanli_*.tar.gz && tar -C '${REMOTE}' -xf - && chmod +x '${REMOTE}'/*.slurm && ls -la '${REMOTE}'"
}

# The R module's Makeconf has -march=native, the library is on Lustre, and
# srun lands on whatever node is free: on 2026-10-03 that was an EPYC 9355,
# so stanli's small R bridge (stanli.so) got an AVX-512 instruction
# (vcvttsd2usi, the seed cast in stanli_r_model_new) and every task on a
# node without AVX-512 - the EPYC-Genoa VM, the EPYC 7513 nodes - died with
# SIGILL in stanli_model(). A Makevars with a fixed -march, put first via
# R_MAKEVARS_USER, makes whatever this installs run on every node.
# x86-64-v3 (AVX2) is what the VM and the 7513s have; nothing compiled here
# is on a hot path. The file is written on the login side so that no quoting
# happens inside the srun command.
cmd_install() {
  info "installing brms + stanli into ${LIB} (compute node, a few minutes)"
  printf 'CFLAGS = -O2 -march=x86-64-v3\nCXXFLAGS = -O2 -march=x86-64-v3\nCXX11FLAGS = -O2 -march=x86-64-v3\nCXX14FLAGS = -O2 -march=x86-64-v3\nCXX17FLAGS = -O2 -march=x86-64-v3\n' \
    | remote "cat > '${REMOTE}/Makevars.portable'"
  remote "cd '${REMOTE}' && srun --partition=short --time=00:45:00 --cpus-per-task=4 --mem=8G --job-name=cogmod_stanli_install bash -c '
    . /etc/profile.d/lmod.sh; module use \"${EB_MODULES}\"; module load \"${R_MODULE}\"
    export STANLI_LIB=\"${LIB}\" STANLI_DIR=\"${REMOTE}\" MAKEFLAGS=-j4; export R_LIBS=\"${LIB}:${PROD_LIB}\"
    export R_MAKEVARS_USER=\"${REMOTE}/Makevars.portable\"
    Rscript install.R'"
}

cmd_submit() {
  local kind="${1:?usage: run.sh submit fit|grad [n] [sbatch args]}"; shift || true
  local n cpus time
  case "$kind" in
    fit|fitpar) n="${1:-16}"; cpus=4; time="02:00:00" ;;
    grad) n="${1:-4}";  cpus=2; time="01:30:00" ;;
    *) die "unknown kind: $kind" ;;
  esac
  shift || true
  # n is a count (1..n) or a range of indices ("10-12") to rerun those tasks.
  local range="$n"; case "$n" in *-*) ;; *) range="1-${n}" ;; esac
  info "submitting ${kind}: tasks ${range}, ${cpus} CPUs each, at most ${MAX} at once"
  remote "mkdir -p '${SCRATCH}' '${REMOTE}/results' && cd '${REMOTE}' && sbatch \
    --array=${range}%${MAX} --cpus-per-task=${cpus} --mem=16G --time=${time} \
    --job-name='cogmod_stanli_${kind}' \
    --output='${SCRATCH}/${kind}_%A_%a.out' --error='${SCRATCH}/${kind}_%A_%a.err' \
    --chdir='${REMOTE}' --export=ALL,STANLI_KIND='${kind}',${EXPORTS} \
    $* task.slurm"
}

cmd_queue() { remote "squeue -u '${USER_}' -o '%.14i %.10P %.24j %.8T %.10M %.10L %.5C %R' | grep -E 'JOBID|cogmod_stanli' || echo '(no stanli jobs queued)'"; }

cmd_progress() {
  remote "cd '${REMOTE}/results' 2>/dev/null || { echo 'no results yet'; exit 0; }
    echo \"fit : \$(ls fit_*/fit_summary.csv 2>/dev/null | wc -l) seeds done\"
    echo \"fitpar: \$(ls fitpar_*/fit_summary_par.csv 2>/dev/null | wc -l) seeds done\"
    echo \"grad: \$(ls grad_*/time.csv 2>/dev/null | wc -l) done\"
    echo '--- failed tasks (.err with Error):'; grep -l '^Error\|Execution halted' '${SCRATCH}'/*.err 2>/dev/null | wc -l"
}

cmd_log() {
  local key="${1:-}"
  remote "cd '${SCRATCH}' 2>/dev/null || exit 1; for f in \$(ls -t *${key}*.out *${key}*.err 2>/dev/null | head -4); do echo \"===== \$f =====\"; tail -n 25 \"\$f\"; done"
}

cmd_pull() {
  local dest="$HERE/../results/hpc"
  mkdir -p "$dest"
  info "pulling results into ${dest}"
  remote "cd '${REMOTE}/results' && tar -cf - --exclude='*.log' --exclude='*.rds' ." | tar -C "$dest" -xf -
  ls "$dest" | wc -l
}

cmd_cancel() { remote "scancel '${1:?jobid}' && echo cancelled"; }
cmd_sh() { remote "${1:?cmd}"; }

case "${1:-}" in
  push)     shift; cmd_push "$@" ;;
  install)  shift; cmd_install "$@" ;;
  submit)   shift; cmd_submit "$@" ;;
  queue)    shift; cmd_queue "$@" ;;
  progress) shift; cmd_progress "$@" ;;
  log)      shift; cmd_log "$@" ;;
  pull)     shift; cmd_pull "$@" ;;
  cancel)   shift; cmd_cancel "$@" ;;
  sh)       shift; cmd_sh "$@" ;;
  *) sed -n '2,/^# ----/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
