#!/bin/bash
# Run on a compute node by `run.sh install` (srun), in the experiment's
# directory: the R module, a Rust toolchain for nutpieR, then install.R.
#
# nutpieR has no binary for R 4.3 on Linux, so it builds from source, which
# needs cargo and rustc >= 1.85. EasyBuild has Rust 1.86 and 1.88, but on
# GCCcore 13.3 and 14.3, against the R module's 12.3 (foss-2023a): loading
# one next to R invites Lmod to swap GCCcore under R. So rustup installs a
# minimal stable toolchain into the experiment's directory instead (about
# 0.5 GB, once; nothing in $HOME), unless NP_RUST_MODULE names a module. nutpieR's Makevars points CARGO_HOME at a
# temporary directory for the build itself; the rustup proxies on PATH find
# the toolchain through RUSTUP_HOME, which it leaves alone. Compute nodes
# reach crates.io and GitHub, as the stanli install already relied on.
set -uo pipefail
cd "${NP_DIR:?}" || exit 1
[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${NP_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${NP_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"

if [ -n "${NP_RUST_MODULE:-}" ]; then
  module load "$NP_RUST_MODULE"
else
  export RUSTUP_HOME="${NP_DIR}/rust/rustup" CARGO_HOME="${NP_DIR}/rust/cargo"
  if [ ! -x "${CARGO_HOME}/bin/cargo" ]; then
    curl --proto '=https' --tlsv1.2 -sSf https://sh.rustup.rs \
      | sh -s -- -y --profile minimal --default-toolchain stable --no-modify-path || exit 1
  fi
  export PATH="${CARGO_HOME}/bin:${PATH}"
fi
echo "rustc: $(rustc --version 2>&1)"
echo "cargo: $(cargo --version 2>&1)"

# The runtime run.sh pushed (sha256-checked), unpacked once; never
# stanli_install() (run.sh, cmd_push, says why).
rt="${NP_DIR}/runtime_v${NP_STANLI_VERSION:?}"
if [ ! -f "$rt/libstanli.so" ]; then
  mkdir -p "$rt" && tar -xzf stanli-runtime-linux-x86_64.tar.gz -C "$rt" || exit 1
fi
export STANLI_RUNTIME="$rt/libstanli.so"
echo "runtime: $(ls -l "$STANLI_RUNTIME")"

export MAKEFLAGS=-j4 R_MAKEVARS_USER="${NP_DIR}/Makevars.portable"
export R_LIBS="${NP_LIB}:${NP_PROD_LIB}"
Rscript install.R
