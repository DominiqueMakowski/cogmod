#!/bin/bash
#SBATCH --job-name=cogmod_lba_screen
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --partition=general
#SBATCH --output=screen_%A_%a.out
#SBATCH --error=screen_%A_%a.err

# One array task = one chain of the LBA screen (see cells.R). Everything that
# varies is on the sbatch command line from run.sh.
#
# Module lines follow analysis/server/fit.slurm of the Illusion Game project.
# The screen's own R library goes FIRST on R_LIBS so that its cogmod (this
# tree) shadows the production one; fit_cell.R removes it again for the
# 0.3.3 arm, which wants the production cogmod. The production library is
# never written to.

cd "${SCREEN_DIR:-$SLURM_SUBMIT_DIR}" || exit 1
echo "Running on $(hostname) at $(date '+%Y-%m-%d %H:%M:%S')"
echo "Experiment ${SCREEN_EXPERIMENT:-lba} cell ${SLURM_ARRAY_TASK_ID:-1} -> ${SCREEN_OUT}"

[ -n "${LMOD_CMD:-}" ] || [ ! -f /etc/profile.d/lmod.sh ] || . /etc/profile.d/lmod.sh
module use "${SCREEN_EB_MODULES:-/mnt/shared/easybuild/modules/all}"
module load "${SCREEN_R_MODULE:-CmdStanR/0.7.1-foss-2023a-R-4.3.2}"
export R_LIBS="${SCREEN_LIB}:${SCREEN_PROD_LIB}${R_LIBS:+:$R_LIBS}"
echo "Rscript: $(command -v Rscript)"

Rscript fit_cell.R
status=$?
echo "FINISHED at $(date '+%Y-%m-%d %H:%M:%S') with status $status"
exit $status
