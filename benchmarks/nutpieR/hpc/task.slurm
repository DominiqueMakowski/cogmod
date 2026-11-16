#!/bin/bash
#SBATCH --job-name=cogmod_nutpie
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=8G
#SBATCH --partition=general
#SBATCH --output=nutpie_%A_%a.out
#SBATCH --error=nutpie_%A_%a.err

# One array task of the nutpie benchmark: bench.R for one (model, seed), with
# NP_ARGS passed on. Task k takes model k mod M of NP_MODELS (M of them,
# colon-separated: sbatch --export splits at commas) and seed
# (k - 1) div M + 1, so tasks 1-M are one seed of every model. Everything
# else that varies is on the sbatch command line from run.sh.
#
# Each task works in its own copy of the pushed tree on scratch, so no two
# tasks share an output directory or a compiled executable, then copies what
# it wrote to results/fit[_<tag>]_<seed>/ with a note of the node. The
# experiment's R library (brms 2.23.1, stanli 0.26.0, nutpieR) goes first on
# R_LIBS, then production, which is never written to.

idx="${SLURM_ARRAY_TASK_ID:-1}"
job="${SLURM_ARRAY_JOB_ID:-0}"
tag="${NP_TAG:-}"
IFS=: read -ra models <<< "${NP_MODELS:-lnr}"
model="${models[$(( (idx - 1) % ${#models[@]} ))]}"
seed=$(( (idx - 1) / ${#models[@]} + 1 ))
work="${NP_SCRATCH}/work/fit_${job}_${idx}"
dest="${NP_DIR}/results/fit${tag:+_$tag}_${model}_${seed}"
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S'): ${model} seed ${seed} -> ${dest}"

[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${NP_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${NP_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
export R_LIBS="${NP_LIB}:${NP_PROD_LIB}${R_LIBS:+:$R_LIBS}"
# The runtime install.sh unpacked, never ~/.cache (run.sh, cmd_push).
export STANLI_RUNTIME="${NP_DIR}/runtime_v${NP_STANLI_VERSION}/libstanli.so"
echo "Rscript: $(command -v Rscript)"

rm -rf "$work" && mkdir -p "$work" "$dest" && cp -r "${NP_DIR}/tree/." "$work/" || exit 1
cd "$work" || exit 1
{ echo "node: $(hostname)"; echo "cpus: $(nproc)"; grep -m1 'model name' /proc/cpuinfo
  echo "job: ${job}_${idx} started $(date '+%Y-%m-%d %H:%M:%S')"; } > "$dest/node.txt"

# NP_ARGS unquoted on purpose: it is a list of options.
Rscript benchmarks/nutpieR/bench.R --model "$model" --seed "$seed" ${tag:+--tag $tag} ${NP_ARGS:-} \
  --out "$work/out" > "$dest/fit.log" 2>&1
status=$?
cp benchmarks/nutpieR/results/* "$dest/" 2>/dev/null
echo "finished $(date '+%Y-%m-%d %H:%M:%S') status $status" >> "$dest/node.txt"
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
rm -rf "$work"
exit $status
