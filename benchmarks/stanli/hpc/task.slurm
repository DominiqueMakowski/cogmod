#!/bin/bash
#SBATCH --job-name=cogmod_stanli
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --output=stanli_%A_%a.out
#SBATCH --error=stanli_%A_%a.err

# One array task of the stanli benchmark; everything that varies is on the
# sbatch command line from run.sh.
#
#   STANLI_KIND=fit   bench.R's 2x2 fit at --seed = the array index
#                     (STANLI_ARMS, if set, replaces the default arms)
#   STANLI_KIND=grad  bench.R grad, bisect.R, repro_issue.R and
#                     repro_parallel_init.R: the per-gradient checks and
#                     timings, single-threaded, then the parallel-chains check
#
# Each task works in its own copy of the pushed tree on scratch, so no two
# tasks share an output directory or a compiled executable, and then copies
# what it wrote to results/<kind>_<index>/, with a note of the node. The
# experiment's R library goes first on R_LIBS so that its brms 2.23.1 and
# stanli shadow the production library, which is never written to.

idx="${SLURM_ARRAY_TASK_ID:-1}"
kind="${STANLI_KIND:-fit}"
tag="${STANLI_TAG:-}"
work="${STANLI_SCRATCH}/work/${kind}_${SLURM_ARRAY_JOB_ID:-0}_${idx}"
dest="${STANLI_DIR}/results/${kind}${tag:+_$tag}_${idx}"
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S'): ${kind} ${idx} -> ${dest}"

[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${STANLI_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${STANLI_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
export R_LIBS="${STANLI_LIB}:${STANLI_PROD_LIB}${R_LIBS:+:$R_LIBS}"
echo "Rscript: $(command -v Rscript)"

rm -rf "$work" && mkdir -p "$work" "$dest" && cp -r "${STANLI_DIR}/tree/." "$work/" || exit 1
cd "$work" || exit 1
{ echo "node: $(hostname)"; echo "cpus: $(nproc)"; grep -m1 'model name' /proc/cpuinfo
  echo "job: ${SLURM_ARRAY_JOB_ID:-}_${idx} started $(date '+%Y-%m-%d %H:%M:%S')"; } > "$dest/node.txt"

R=benchmarks/stanli
Rscript $R/bench.R emit > "$dest/emit.log" 2>&1
status=$?
if [ $status -eq 0 ] && [ "$kind" = fit ]; then
  Rscript $R/bench.R fit --seed "$idx" ${tag:+--tag $tag} ${STANLI_ARMS:+--arms $STANLI_ARMS} \
    > "$dest/fit.log" 2>&1; status=$?
elif [ $status -eq 0 ] && [ "$kind" = grad ]; then
  Rscript $R/bench.R grad --reps 21 ${tag:+--tag $tag} > "$dest/grad.log" 2>&1; status=$?
  Rscript $R/bisect.R ${tag:+--tag $tag} > "$dest/bisect.log" 2>&1 || status=$?
  Rscript $R/repro_issue.R > "$dest/repro_issue.txt" 2>&1 || status=$?
  # Chain times are each chain's own wall clock, so four chains sharing two
  # cores slow alike and wall/sum stays near 1/4 when they overlap; 1.0 is
  # serial.
  Rscript $R/repro_parallel_init.R > "$dest/repro_parallel_init.txt" 2>&1 || status=$?
fi
cp $R/results/* "$dest/" 2>/dev/null
echo "finished $(date '+%Y-%m-%d %H:%M:%S') status $status" >> "$dest/node.txt"
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
rm -rf "$work"
exit $status
