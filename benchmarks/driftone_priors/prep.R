# Build push/ for the cluster: the vignette data as plain data frames (no
# brmsfit, so nothing on the cluster needs rstan to read them), the cell
# table, the package tree being tested, and the two scripts. Run from the
# package root:
#
#   Rscript benchmarks/driftone_priors/prep.R
#
# Arm A is whatever this tree's priors are; B and C patch the registries
# in-session (refit_cell.R), so to rerun the 2026-10-12 comparison as it was,
# push a tree of 0.3.4 (git archive 5382c79), whose priors arm A measured.
here <- "benchmarks/driftone_priors"
P <- file.path(here, "push")
unlink(P, recursive = TRUE)
dir.create(file.path(P, "data"), recursive = TRUE)
dir.create(file.path(P, "tree"))

# The decision_making article's data (two participants of rtdists' speed_acc,
# errors kept) and the rt_models article's (three participants, correct
# responses only), exactly as the cached fits had them.
dm <- as.data.frame(readRDS("vignettes/models/m_rdm.rds")$data)[, c("RT", "Error", "Condition")]
rt <- as.data.frame(readRDS("vignettes/models/m_wald4.rds")$data)[, c("RT", "Condition")]
# The stress case: 70% of the trials, which leaves the Accuracy condition 71
# errors - the ~75 on which a flat `driftone` prior ran the RDM to the floor
set.seed(2026)
sub <- dm[sort(sample(nrow(dm), round(0.7 * nrow(dm)))), ]
saveRDS(dm, file.path(P, "data", "dm.rds"), version = 2)
saveRDS(sub, file.path(P, "data", "sub.rds"), version = 2)
saveRDS(rt, file.path(P, "data", "rt.rds"), version = 2)

# Most informative first, since the array runs in row order: the three race
# fits under all three arms, then the threshold and outlier checks under A, B.
cells <- rbind(
  expand.grid(seed = 1:4, arm = c("A", "B", "C"), model = "rdm_sub", stringsAsFactors = FALSE),
  expand.grid(seed = 1:4, arm = c("A", "B", "C"), model = "rdm", stringsAsFactors = FALSE),
  expand.grid(seed = 1:3, arm = c("A", "B", "C"), model = "lba2", stringsAsFactors = FALSE),
  expand.grid(seed = 1:2, arm = c("A", "B"),
              model = c("wald4", "exwald", "ddm", "ddm_bias", "ddm_pout", "bisa"),
              stringsAsFactors = FALSE)
)[, c("model", "arm", "seed")]
utils::write.csv(cells, file.path(P, "cells.csv"), row.names = FALSE)

file.copy(c("DESCRIPTION", "NAMESPACE"), file.path(P, "tree"))
for (d in c("R", "data", "inst")) file.copy(d, file.path(P, "tree"), recursive = TRUE)
file.copy(file.path(here, c("refit_cell.R", "task.slurm")), P)
cat(nrow(cells), "cells written to", P, "\n")
