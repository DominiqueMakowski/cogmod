# Benchmark: can a fit on a subset of participants warm-start the full fit?
#
# Companion to the "Improving Sampling Efficiency and Performance" article
# (vignettes/articles/performance.qmd), "Warm Starts" section, and the
# integration check of cogmod_warmstart(). Not part of the R package; run it
# from the repository root:
#
#   Rscript benchmarks/warm_start_subset.R
#
# or for another family / size, from an R session:
#
#   bench_overrides <- list(family = "ddm"); source("benchmarks/warm_start_subset.R")
#
# The question: a pilot fit on some of the participants has adapted a metric
# and a step size, and has posterior means for every parameter. The full data
# set has more participants and so more parameters (one standardized random
# effect per participant per term), so Stan refuses the pilot's metric as is -
# the lengths differ. cogmod_warmstart() maps it across by parameter name:
# population-level entries are carried over, a pilot participant's `z` entries
# follow it to its position in the full model, and each new participant gets
# the average of the pilot's entries for the same effect. Initial values come
# from the pilot's posterior means (new participants start at z = 0).
#
# Model: the family's "mixed" formula of helpers.R (participant random
# intercepts on two or three parameters), fitted with brm(). Runs on the full
# data, all with 500 sampling iterations:
#   reference        cogmod inits, warmup 500
#   cold             cogmod inits, warmup 100
#   warm inits       pilot posterior means as inits, default metric, warmup 100
#   warm             + pilot metric and step size (cogmod_warmstart()), warmup 100
#
# Output: benchmarks/results/warm_start_subset_<family>.csv (one row per run,
# the pilot fit included: wall time, min bulk ESS, ESS per second,
# divergences, max Rhat).

source("benchmarks/helpers.R")

defaults <- list(
  family = "lnr",
  n_participants = 8,
  n_subset = 4,
  n_trials = 150,
  chains = 4,
  cores = 4,
  out_dir = file.path("benchmarks", "results")
)
settings <- if (exists("bench_overrides")) modifyList(defaults, bench_overrides) else defaults
list2env(settings, envir = environment())
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
out_file <- file.path(out_dir, sprintf("warm_start_subset_%s.csv", family))

df <- bench_data(n_participants, n_trials)
sub_ids <- levels(df$Participant)[seq_len(n_subset)]
df_sub <- droplevels(df[df$Participant %in% sub_ids, ])
cat(sprintf("Family: %s. Pilot subset: %d participants, %d trials\n", family, n_subset, nrow(df_sub)))

runs <- list()
record <- function(label, m, warmup, n_part) {
  wall <- max(rowSums(rstan::get_elapsed_time(m$fit)))
  er <- ess_rhat(m)
  row <- data.frame(family = family, run = label, participants = n_part, warmup = warmup,
                    wall_s = round(wall, 1),
                    ess_bulk_min = round(unname(er["ess_bulk_min"])),
                    ess_per_s = round(unname(er["ess_bulk_min"]) / wall, 2),
                    divergent = sum(nuts_params(m, pars = "divergent__")$Value),
                    rhat_max = round(unname(er["rhat_max"]), 3))
  print(row, row.names = FALSE)
  runs[[length(runs) + 1]] <<- row
  write.csv(do.call(rbind, runs), out_file, row.names = FALSE)
}

f <- make_formula(family, "mixed")
prior <- cogmod_priors(f, df)
stanvars <- cogmod_stanvars(f)

# ---- The pilot fit -----------------------------------------------------------------------

cat("\n==== Pilot fit on the subset ====\n")
m_pilot <- suppressMessages(
  brm(f, data = df_sub, prior = prior, stanvars = stanvars, init = cogmod_inits(f, df_sub),
      backend = "cmdstanr", chains = chains, cores = cores,
      iter = 1000, warmup = 500, seed = 1, refresh = 0, silent = 2)
)
record("pilot (subset)", m_pilot, 500, n_subset)

ws <- cogmod_warmstart(m_pilot, data = df)   # the pilot's model, on all the data
print(ws)

# ---- Runs on the full data ---------------------------------------------------------
# update() reuses the compiled program: the number of participants is data.

refit <- function(init, warmup, seed, ...) {
  suppressMessages(
    update(m_pilot, newdata = df, init = init, chains = chains, cores = cores,
           iter = warmup + 500, warmup = warmup, seed = seed, refresh = 0, silent = 2, ...)
  )
}

cat("\n==== Full data ====\n")
m_ref <- refit(cogmod_inits(f, df), 500, 1)
record("reference", m_ref, 500, n_participants)

m_cold <- refit(cogmod_inits(f, df), 100, 2)
record("cold, short warmup", m_cold, 100, n_participants)

m_init <- refit(ws$init, 100, 2)
record("warm inits only", m_init, 100, n_participants)

m_warm <- refit(ws$init, 100, 2, inv_metric = ws$inv_metric, step_size = ws$step_size)
record("warm inits + metric + step size", m_warm, 100, n_participants)

# ---- How far was the pilot's adaptation from the full fit's? ---------------------------

ws_ref <- cogmod_warmstart(m_ref)
pop <- is.na(ws$table$group)
cat("\nInverse metric (posterior variances), pilot mapped vs full reference, population-level:\n")
cmp <- rbind(pilot = ws$inv_metric[pop], reference = ws_ref$inv_metric[pop])
colnames(cmp) <- ws$table$parameter[pop]
print(round(cmp, 4))
cat(sprintf("z entries, mean: pilot %.3f, reference %.3f\n",
            mean(ws$inv_metric[!pop]), mean(ws_ref$inv_metric[!pop])))
cat(sprintf("step size: pilot %.3f, reference %.3f\n", ws$step_size, ws_ref$step_size))

cat("\nPopulation-level estimates (should agree):\n")
print(round(cbind(pilot = fixef(m_pilot)[, 1], reference = fixef(m_ref)[, 1], cold = fixef(m_cold)[, 1],
                  warm_inits = fixef(m_init)[, 1], warm = fixef(m_warm)[, 1]), 3))

cat("\n==== Summary ====\n")
print(do.call(rbind, runs), row.names = FALSE)
cat(sprintf("\nResults written to %s\n", normalizePath(out_file)))
