#!/usr/bin/env bash
# ---------------------------------------------------------------------------
# run.sh -- the nutpie benchmark (benchmarks/nutpieR/bench.R) on Artemis over SSH.
#
#   ./run.sh push            the scripts, the minimal package tree, and stanli
#                            $NP_STANLI_VERSION's (default 0.26.0) source and Linux
#                            runtime, sha256-checked -> the cluster (one tar pipe)
#   ./run.sh install         Rust (rustup, into the experiment's directory), then
#                            brms 2.23.1, stanli and nutpieR (from source) into the
#                            experiment's own R library and the runtime into
#                            runtime_v<version>/, on a compute node; checks every piece
#   ./run.sh smoke           one short task per model in NP_MODELS (--iter 100, the
#                            default arms), waited on, before submitting
#   ./run.sh submit [n]      n array tasks of bench.R, 4 whole cores each, one
#                            (model, seed) per task: task k is model k mod M of
#                            NP_MODELS (colon-separated, default "lnr"), seed
#                            (k - 1) div M + 1. Default n = 18 seeds per model (every
#                            arm order 6 times over 3 arms); n may be a range
#                            ("10-12") to rerun those indices.
#                            NP_ARGS goes to bench.R (e.g. "--arms cmdstanr:stanli:nutpie:nutpie_lowrank",
#                            colon-separated: sbatch --export splits at commas),
#                            NP_TAG keeps a run's results apart
#                            (results/fit_<tag>_<model>_<seed>/), NP_TIME replaces
#                            the 01:30:00 limit; extra args go to sbatch
#   ./run.sh queue           squeue for this user
#   ./run.sh progress        how many seeds of each model have finished, and failures
#   ./run.sh log [pattern]   tail the newest .out/.err (optionally matching)
#   ./run.sh pull            copy results/ back into ../results/hpc/
#   ./run.sh cancel <jobid>  scancel
#   ./run.sh sh '<cmd>'      run a command on the login node
#
# Modelled on benchmarks/stanli/hpc/run.sh, with its own directories and its
# own R library: brms 2.23.1 from GitHub (as on the laptop; the R module has
# 2.21.0), stanli and nutpieR. R_LIBS puts it ahead of production, which is
# never written to. (The stanli experiment's brms 2.23.1 library is gone
# from the cluster; cogmod_stanli/ now holds a stanli 0.25.0 one.)
#
# Every task gets 4 whole cores (--hint=nomultithread: Artemis counts
# hyperthreads as CPUs, so a bare --cpus-per-task=4 is two cores, which
# hurt stanli more than CmdStan; benchmarks/stanli/README.md, "Two chains
# per core"), so a task takes 8 CPUs off the account's 550-CPU cap.
# artemis-rtx-01 (slow CmdStan) and the general-01/-02 VMs are excluded, as
# in the stanli runs. The GlobalProtect VPN must be up.
# ---------------------------------------------------------------------------
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PKG="$(cd "$HERE/../../.." && pwd)"

ALIAS="${NP_HPC_ALIAS:-artemis}"
USER_="${NP_HPC_USER:-dmm56}"
GROUP_="${NP_HPC_GROUP:-psych}"
REMOTE="${NP_REMOTE:-/mnt/lustre/users/${GROUP_}/${USER_}/cogmod_nutpie}"
SCRATCH="${NP_SCRATCH:-/mnt/lustre/scratch/${GROUP_}/${USER_}/cogmod_nutpie}"
LIB="${REMOTE}/R_libs"
PROD_LIB="${NP_PROD_LIB:-/mnt/lustre/users/${GROUP_}/${USER_}/cluster_R_libs/x86_64-pc-linux-gnu-library/4.3}"
R_MODULE="${NP_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
EB_MODULES="${NP_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
# A module with rustc >= 1.85 in place of rustup (install.sh says why not by default).
RUST_MODULE="${NP_RUST_MODULE:-}"
MAX="${NP_MAX:-20}"
STANLI_VERSION="${NP_STANLI_VERSION:-0.26.0}"
TAG="${NP_TAG:-}"
MODELS="${NP_MODELS:-lnr}"
N_MODELS="$(tr ':' '\n' <<< "$MODELS" | grep -c .)"
EXCLUDE="${NP_EXCLUDE:-artemis-rtx-01,artemis-general-01,artemis-general-02}"

if [ -x /c/Windows/System32/OpenSSH/ssh.exe ]; then
  SSH=/c/Windows/System32/OpenSSH/ssh.exe
else
  SSH=ssh
fi
die() { printf '\033[31merror:\033[0m %s\n' "$*" >&2; exit 1; }
info() { printf '\033[36m==>\033[0m %s\n' "$*"; }
remote() { "$SSH" -o BatchMode=yes "$ALIAS" "$@" 2> >(grep -v 'module: command not found' >&2); }
# The environment a task gets: exports TAG ARGS.
exports() {
  printf "%s" "NP_DIR='${REMOTE}',NP_SCRATCH='${SCRATCH}',NP_LIB='${LIB}'," \
    "NP_PROD_LIB='${PROD_LIB}',NP_R_MODULE='${R_MODULE}',NP_EB_MODULES='${EB_MODULES}'," \
    "NP_RUST_MODULE='${RUST_MODULE}',NP_STANLI_VERSION='${STANLI_VERSION}',NP_MODELS='${MODELS}'," \
    "NP_TAG='$1',NP_ARGS='$2'"
}
EXPORTS="$(exports "$TAG" "${NP_ARGS:-}")"

# A release asset of stanli into $1, checked against the sha256 GitHub lists
# for it.
fetch_checked() {
  local dir="$1" name="$2" want got
  curl -sSfL -o "$dir/$name" "https://github.com/seantalts/stanli/releases/download/v${STANLI_VERSION}/${name}"
  # Each asset's "digest" comes some 25 lines after its "name" (the uploader
  # block sits between), so keep those two keys only, then pair them up.
  want="$(curl -sSfL "https://api.github.com/repos/seantalts/stanli/releases/tags/v${STANLI_VERSION}" \
    | grep -E '"(name|digest)": ' | grep -A1 "\"name\": \"${name}\"" \
    | sed -n 's/.*"digest": "sha256:\([0-9a-f]*\)".*/\1/p' | head -1)"
  got="$(sha256sum "$dir/$name" | cut -d' ' -f1)"
  [ -n "$want" ] && [ "$want" = "$got" ] || die "${name}: sha256 ${got} does not match GitHub's '${want}'"
  info "${name}: sha256 ok"
}

# What load_all() and the benchmark need, nothing else: the package's R code
# and metadata, gradient_programs.R (gp_args(), gp_load()) and bench.R. The
# runtime is not fetched by stanli_install() on the cluster: that deletes
# every other release's cached runtime, and IGComputational runs stanli
# 0.25.0 from its own copy (cogmod_stanli/src_0.25.0/install_stanli_0.25.0.sh
# does the same as here). install.sh unpacks it into runtime_v<version>/,
# which STANLI_RUNTIME points at; nothing goes into ~/.cache.
cmd_push() {
  info "pushing to ${REMOTE}"
  local tmp; tmp="$(mktemp -d)"; trap 'rm -rf "$tmp"' RETURN
  mkdir -p "$tmp/tree/benchmarks/nutpieR"
  (cd "$PKG" && cp -r DESCRIPTION NAMESPACE R data inst "$tmp/tree/")
  sed 's/\r$//' "$PKG/benchmarks/gradient_programs.R" > "$tmp/tree/benchmarks/gradient_programs.R"
  sed 's/\r$//' "$PKG/benchmarks/nutpieR/bench.R" > "$tmp/tree/benchmarks/nutpieR/bench.R"
  for f in "$HERE"/*.R "$HERE"/*.slurm "$HERE"/install.sh; do sed 's/\r$//' "$f" > "$tmp/$(basename "$f")"; done
  fetch_checked "$tmp" "stanli_${STANLI_VERSION}.tar.gz"
  fetch_checked "$tmp" "stanli-runtime-linux-x86_64.tar.gz"
  # install.R takes the one tarball it finds, so a previous version's goes.
  tar -C "$tmp" -cf - . | remote "mkdir -p '${REMOTE}' '${SCRATCH}' && rm -rf '${REMOTE}/tree' '${REMOTE}'/stanli_*.tar.gz '${REMOTE}'/stanli-runtime-*.tar.gz && tar -C '${REMOTE}' -xf - && chmod +x '${REMOTE}'/*.slurm '${REMOTE}'/install.sh && ls -la '${REMOTE}'"
}

# stanli's R bridge and nutpieR's C entry point are compiled with the R
# module's Makeconf, which has -march=native: built on an AVX-512 node, they
# die with SIGILL on the nodes without it (benchmarks/stanli/hpc/run.sh,
# cmd_install). A Makevars with a fixed -march, first via R_MAKEVARS_USER,
# makes the library run on every node. Rust builds for baseline x86-64
# unless told otherwise, and nutpieR does not. Written on the login side so
# that no quoting happens inside the srun command.
cmd_install() {
  info "installing Rust, stanli ${STANLI_VERSION} and nutpieR into ${REMOTE} (compute node, 10-20 minutes)"
  printf 'CFLAGS = -O2 -march=x86-64-v3\nCXXFLAGS = -O2 -march=x86-64-v3\nCXX11FLAGS = -O2 -march=x86-64-v3\nCXX14FLAGS = -O2 -march=x86-64-v3\nCXX17FLAGS = -O2 -march=x86-64-v3\n' \
    | remote "cat > '${REMOTE}/Makevars.portable'"
  remote "cd '${REMOTE}' && srun --partition=short --time=01:00:00 --cpus-per-task=4 --mem=8G \
    --job-name=cogmod_nutpie_install --export=ALL,${EXPORTS} bash install.sh"
}

# One task per model (seed 1), in parallel; sbatch --wait returns when all
# have finished, then each task's status line and the tail of its fit.log.
cmd_smoke() {
  info "smoke: ${MODELS}, --iter 100 (short partition, waited on)"
  local ex; ex="$(exports smoke '--iter 100')"
  remote "mkdir -p '${SCRATCH}' '${REMOTE}/results' && cd '${REMOTE}' && sbatch --wait \
    --array=1-${N_MODELS} --partition=short --time=01:00:00 --cpus-per-task=4 --hint=nomultithread \
    --mem=8G --exclude='${EXCLUDE}' --job-name=cogmod_nutpie_smoke \
    --output='${SCRATCH}/smoke_%A_%a.out' --error='${SCRATCH}/smoke_%A_%a.err' \
    --chdir='${REMOTE}' --export=ALL,${ex} task.slurm; \
    for m in \$(tr ':' ' ' <<< '${MODELS}'); do d=results/fit_smoke_\${m}_1; echo \"===== \$m: \$(tail -n 1 \$d/node.txt)\"; \
      grep -v '^Chain\|^\$' \$d/fit.log | tail -n 12; done"
}

cmd_submit() {
  local n="${1:-$((18 * N_MODELS))}"; shift || true
  # n is a count (1..n) or a range of indices ("10-12") to rerun those tasks.
  local range="$n"; case "$n" in *-*) ;; *) range="1-${n}" ;; esac
  info "submitting tasks ${range} over ${MODELS}, 4 whole cores each, at most ${MAX} at once"
  remote "mkdir -p '${SCRATCH}' '${REMOTE}/results' && cd '${REMOTE}' && sbatch \
    --array=${range}%${MAX} --partition=general --cpus-per-task=4 --hint=nomultithread --mem=8G \
    --time=${NP_TIME:-01:30:00} --exclude='${EXCLUDE}' --job-name=cogmod_nutpie \
    --output='${SCRATCH}/fit_%A_%a.out' --error='${SCRATCH}/fit_%A_%a.err' \
    --chdir='${REMOTE}' --export=ALL,${EXPORTS} $* task.slurm"
}

cmd_queue() { remote "squeue -u '${USER_}' -o '%.14i %.10P %.24j %.8T %.10M %.10L %.5C %R' | grep -E 'JOBID|cogmod_nutpie' || echo '(no nutpie jobs queued)'"; }

cmd_progress() {
  remote "cd '${REMOTE}/results' 2>/dev/null || { echo 'no results yet'; exit 0; }
    for m in \$(tr ':' ' ' <<< '${MODELS}'); do
      echo \"\$m${TAG:+ $TAG}: \$(ls fit_${TAG:+${TAG}_}\${m}_[0-9]*/fit_*.csv 2>/dev/null | wc -l) seeds done\"; done
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
  remote "cd '${REMOTE}/results' && tar -cf - --exclude='*.rds' ." | tar -C "$dest" -xf -
  ls "$dest" | wc -l
}

cmd_cancel() { remote "scancel '${1:?jobid}' && echo cancelled"; }
cmd_sh() { remote "${1:?cmd}"; }

case "${1:-}" in
  push)     shift; cmd_push "$@" ;;
  install)  shift; cmd_install "$@" ;;
  smoke)    shift; cmd_smoke "$@" ;;
  submit)   shift; cmd_submit "$@" ;;
  queue)    shift; cmd_queue "$@" ;;
  progress) shift; cmd_progress "$@" ;;
  log)      shift; cmd_log "$@" ;;
  pull)     shift; cmd_pull "$@" ;;
  cancel)   shift; cmd_cancel "$@" ;;
  sh)       shift; cmd_sh "$@" ;;
  *) sed -n '2,/^# ----/p' "${BASH_SOURCE[0]}" | sed 's/^# \{0,1\}//' ;;
esac
