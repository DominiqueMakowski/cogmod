# Step 2 of the LNR start check (README.md), part 1: the data, the formulas
# and the compiled programs, built once so that the two arms' processes
# (fit_arm.R) share them and never compile the same file concurrently. Run
# from the package root:
#
#   Rscript benchmarks/lnr_start/prep.R
#
# Three models, each sampled under arm A (nuone starts at 0.7, the 0.3.5
# value) and arm B (nuone starts at 0):
#
#   vignette       the decision_making article's LNR on its own data (two
#                  participants of rtdists' speed_acc, 4,620 trials, 7.3%
#                  errors): RT | dec(Error) ~ Condition, nuone ~ Condition,
#                  ndt ~ Condition, sigmabias = 0
#   vignette_bias  the same with sigmabias ~ 1, the LBA-with-LogNormal-drifts
#                  variant, since both are in use
#   sim_hier       simulated hierarchical LNR data: 20 participants x 120
#                  trials, a Condition effect on both rates, participant
#                  intercepts on both rates and on ndt, sigmabias = 0.
#                  Generating error rates 7.8% (Accuracy) and 20.6% (Speed)
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(brms))
here <- "benchmarks/lnr_start"
W <- file.path(here, "work")
dir.create(W, showWarnings = FALSE, recursive = TRUE)

dm <- as.data.frame(readRDS("vignettes/models/m_lnr.rds")$data)
dm <- dm[, c("RT", "Error", "Condition")]

# Simulated hierarchical data. The population values follow the vignette's
# posterior (mu 1.16, nuone -1.06 in Accuracy, sigmaone about twice
# sigmazero) with a larger Condition effect on the error accumulator, so that
# the Speed condition has one error in five.
set.seed(2026)
P <- 20L
per <- 60L
u <- stats::rnorm(P, 0, 0.25)
w <- stats::rnorm(P, 0, 0.30)
e <- stats::rnorm(P, 0, 0.10)
sim <- do.call(rbind, lapply(seq_len(P), function(p) {
  do.call(rbind, lapply(c(0, 1), function(s) {
    r <- rcogmod_lnr(per, nuzero = 1.0 + 0.4 * s + u[p],
                     nuone = -0.9 + 1.2 * s + w[p],
                     sigmazero = 0.6, sigmaone = 1.2,
                     ndt = 0.25 * exp(e[p]), sigmabias = 0, poutlier = 0.02)
    data.frame(Participant = factor(p, levels = seq_len(P)),
               Condition = c("Accuracy", "Speed")[s + 1],
               RT = r$rt, Error = r$response)
  }))
}))
rownames(sim) <- NULL
cat(sprintf("sim_hier: %d rows, errors %.1f%% (Accuracy %.1f%%, Speed %.1f%%)\n",
            nrow(sim), 100 * mean(sim$Error),
            100 * mean(sim$Error[sim$Condition == "Accuracy"]),
            100 * mean(sim$Error[sim$Condition == "Speed"])))

spec <- list(
  vignette = list(
    data = dm,
    f = bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1,
           sigmaone ~ 1, sigmabias = 0, ndt ~ Condition,
           family = cogmod_lnr())
  ),
  vignette_bias = list(
    data = dm,
    f = bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1,
           sigmaone ~ 1, sigmabias ~ 1, ndt ~ Condition,
           family = cogmod_lnr())
  ),
  sim_hier = list(
    data = sim,
    f = bf(RT | dec(Error) ~ Condition + (1 | Participant),
           nuone ~ Condition + (1 | Participant), sigmazero ~ 1,
           sigmaone ~ 1, sigmabias = 0, ndt ~ 1 + (1 | Participant),
           family = cogmod_lnr())
  )
)

for (name in names(spec)) {
  sp <- spec[[name]]
  pri <- suppressMessages(suppressWarnings(cogmod_priors(sp$f, sp$data)))
  sv <- cogmod_stanvars(sp$f)
  code <- suppressWarnings(brms::make_stancode(sp$f, data = sp$data,
                                               prior = pri, stanvars = sv))
  file <- cmdstanr::write_stan_file(code, dir = W, basename = name)
  cat(name, ": compiling", basename(file), "\n")
  mod <- cmdstanr::cmdstan_model(file, quiet = TRUE)
  spec[[name]]$stan_file <- file
  spec[[name]]$exe_file <- mod$exe_file()
}
saveRDS(spec, file.path(W, "spec.rds"))
cat("wrote", file.path(W, "spec.rds"), "\n")

# push/ for the cluster (run.sh): the same data and formulas without the
# local program paths (each cell compiles its own), the cell table, the
# package tree, and the two scripts. Most informative first, since the array
# runs in row order: the hierarchical model, then the vignette's, then the
# sigmabias variant; A and B of a seed adjacent so that pairs complete
# together.
P <- file.path(here, "push")
unlink(P, recursive = TRUE)
dir.create(file.path(P, "tree"), recursive = TRUE)
dir.create(file.path(P, "results"))
for (name in names(spec)) spec[[name]][c("stan_file", "exe_file")] <- NULL
saveRDS(spec, file.path(P, "spec.rds"), version = 2)
cells <- expand.grid(arm = c("A", "B"), seed = 1:4,
                     model = c("sim_hier", "vignette", "vignette_bias"),
                     stringsAsFactors = FALSE)[, c("model", "arm", "seed")]
utils::write.csv(cells, file.path(P, "cells.csv"), row.names = FALSE)
file.copy(c("DESCRIPTION", "NAMESPACE"), file.path(P, "tree"))
for (d in c("R", "data", "inst")) file.copy(d, file.path(P, "tree"), recursive = TRUE)
file.copy(file.path(here, c("fit_cell.R", "task.slurm")), P)
cat(nrow(cells), "cells written to", P, "\n")
