# Stack the per-cell CSVs that refit_cell.R writes (pulled into results/cells/
# by run.sh pull) into four files: chains.csv, pars.csv, pp.csv and
# priors.csv, each with model/arm/seed columns. Run from the package root.
#
#   Rscript benchmarks/driftone_priors/collect.R [cells dir] [out dir]
args <- commandArgs(trailingOnly = TRUE)
here <- "benchmarks/driftone_priors"
cells <- if (length(args) >= 1) args[1] else file.path(here, "results", "cells")
out <- if (length(args) >= 2) args[2] else file.path(here, "results")
dir.create(out, showWarnings = FALSE, recursive = TRUE)

stack <- function(what) {
  files <- list.files(cells, pattern = sprintf("__%s\\.csv$", what), full.names = TRUE)
  do.call(rbind, lapply(files, function(f) {
    parts <- strsplit(sub(sprintf("__%s\\.csv$", what), "", basename(f)), "__")[[1]]
    x <- utils::read.csv(f, check.names = FALSE)
    if (!nrow(x)) return(NULL)
    cbind(model = parts[1], arm = parts[2], seed = as.integer(sub("^s", "", parts[3])), x)
  }))
}
for (what in c("chains", "pars", "pp", "priors")) {
  x <- stack(what)
  if (is.null(x)) next
  # one priors table per model and arm: every seed wrote the same one
  if (what == "priors") x <- unique(x[, setdiff(names(x), "seed")])
  utils::write.csv(x, file.path(out, paste0(what, ".csv")), row.names = FALSE)
  cat(what, ":", nrow(x), "rows\n")
}
