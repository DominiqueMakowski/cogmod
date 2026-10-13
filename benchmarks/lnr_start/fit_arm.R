# Step 2 of the LNR start check (README.md): one arm's cells on the laptop,
# in one process. The per-cell work is fit_cell.R's; this only loops. The
# two arms can run side by side (4 chains each): CmdStan's output goes to a
# directory of this arm's own, so each cleans up only after itself.
#
#   Rscript benchmarks/lnr_start/fit_arm.R A [seeds] [models]
#   Rscript benchmarks/lnr_start/fit_arm.R B [seeds] [models]
#
# `seeds` is a comma-separated list, 1,2,3,4 by default; `models` likewise,
# in the order to run them, every model in spec.rds by default. Cells whose
# chain diagnostics already exist are skipped, so an interrupted arm resumes.
args <- commandArgs(trailingOnly = TRUE)
arm <- args[1]
stopifnot(arm %in% c("A", "B"))
seeds <- if (length(args) > 1) as.integer(strsplit(args[2], ",")[[1]]) else 1:4
models <- if (length(args) > 2) strsplit(args[3], ",")[[1]] else NULL
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(brms))
here <- "benchmarks/lnr_start"
source(file.path(here, "fit_cell.R"))
W <- file.path(here, "work")
R <- file.path(here, "results")
dir.create(R, showWarnings = FALSE)
out_dir <- file.path(W, arm)
dir.create(out_dir, showWarnings = FALSE)
spec <- readRDS(file.path(W, "spec.rds"))
if (is.null(models)) models <- names(spec)
stopifnot(all(models %in% names(spec)))

for (name in models) {
  for (seed in seeds) {
    fit_cell(name, arm, seed, spec, stan_dir = W, out_dir = out_dir, resdir = R)
  }
}
cat(format(Sys.time(), "%H:%M:%S"), arm, "arm finished\n")
