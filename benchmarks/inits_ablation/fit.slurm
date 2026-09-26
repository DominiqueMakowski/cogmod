#!/bin/bash
#SBATCH --job-name=cogmod_inits
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --output=inits_%A_%a.out
#SBATCH --error=inits_%A_%a.err

# One array task = one cell of the inits ablation (see cells.R). Everything
# that varies is on the sbatch command line from run.sh: --array, --partition,
# --cpus-per-task, --mem, --output/--error and the ABL_* variables.
#
# Module lines follow analysis/server/fit.slurm of the Illusion Game project,
# which is the configuration measured to work on Artemis (its AGENT.md 3.1,
# 3.2). The ablation's own R library goes FIRST on R_LIBS so that its cogmod
# (this tree, with the data-aware layer) shadows the production one, which is
# left untouched.

cd "${ABL_DIR:-$SLURM_SUBMIT_DIR}" || exit 1
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S')"
echo "Experiment ${ABL_EXPERIMENT:-A} cell ${SLURM_ARRAY_TASK_ID:-1} -> ${ABL_OUT}"

[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${ABL_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${ABL_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
export R_LIBS="${ABL_LIB}:${ABL_PROD_LIB}${R_LIBS:+:$R_LIBS}"
echo "Rscript: $(command -v Rscript)"

Rscript fit_cell.R
status=$?
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
