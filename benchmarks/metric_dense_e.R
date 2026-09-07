# Benchmark: do `metric = "dense_e"` and the Stan compiler optimizations pay off
# for evidence accumulation models?
#
# Companion to the "Improving Sampling Efficiency and Performance" article
# (vignettes/articles/performance.qmd), whose findings come from this script.
# Not part of the R package; run it from the repository root:
#
#   Rscript benchmarks/metric_dense_e.R
#
# or `source()` it from an interactive session. Settings live in `defaults`
# below; override any of them without editing the file by defining a list
# named `bench_overrides` before sourcing, e.g.
#
#   bench_overrides <- list(families = "ddm", structures = "fixed", iter = 400)
#   source("benchmarks/metric_dense_e.R")
#
# What it does. For every family (DDM, LBA, LNR, RDM), model structure
# (population-level only vs. participant random intercepts) and compiler
# setting it compiles one Stan program, then samples it once per metric and
# seed, holding data, priors, inits, chains and iterations fixed. Two things
# vary:
#
#   metric   `diag_e` (Stan's default) learns one variance per parameter during
#            warmup; `dense_e` learns the full covariance, so it can absorb the
#            posterior correlations (boundary vs ndt, drift vs boundary, ...)
#            that make these posteriors slow to sample. Changes how many
#            gradient evaluations an effective draw costs.
#   optim    `default` compiles as brms does out of the box; `all` adds stanc's
#            O1 optimizations (memory layout of automatic differentiation) and
#            CmdStan's STAN_CPP_OPTIMS and STAN_NO_RANGE_CHECKS flags. Changes
#            how long one gradient evaluation takes, nothing else.
#
# What it records, per fit (one row of benchmarks/results/metric_dense_e.csv):
#   compile_s       seconds to compile this (family, structure, optim) program
#   wall_s          longest chain, warmup + sampling, in seconds
#   ess_bulk_min    smallest bulk ESS over all parameters (posterior::ess_bulk)
#   ess_tail_min    idem, tail ESS
#   ess_per_s       ess_bulk_min / wall_s: the headline figure
#   ess_per_grad    ess_bulk_min / post-warmup leapfrog steps: how well the
#                   metric fits the geometry, independent of hardware
#   grad_per_s      post-warmup leapfrog steps per core-second: how fast one
#                   gradient evaluation is, which is what the compiler flags act on
#   leapfrog_mean   mean leapfrog steps per post-warmup iteration
#   treedepth_mean, stepsize_mean, divergent, rhat_max, n_pars
#   max_abs_cor     largest |posterior correlation| between population-level
#                   parameters - the thing dense_e is supposed to exploit
#
# The summary at the end reports every (metric, optim) combination relative to
# the baseline (diag_e, default) for each family x structure. ess_per_s_ratio
# above 1 means the combination was the better choice on this machine for this
# model. For the metric, read ess_per_grad_ratio and leapfrog_ratio; for the
# compiler flags, grad_per_s_ratio.
#
# Runtime. Each (family, structure, optim) compiles once (minutes on Windows)
# and samples length(metrics) * length(seeds) times. With the defaults below
# expect a couple of hours; the LBA is the slow one. Resist shortening
# `warmup` to make it faster: the covariance estimate is the whole point of
# dense_e, and with a couple of hundred warmup iterations it is noise, so a
# short run is biased against dense_e and says nothing about a real fit.

source("benchmarks/helpers.R")

# ---- Settings -----------------------------------------------------------------

defaults <- list(
  families = c("ddm", "lba2", "lnr", "rdm"),
  structures = c("fixed", "mixed"),   # "fixed": no random effects; "mixed": (1 | Participant)
  metrics = c("diag_e", "dense_e"),
  optims = c("default", "all"),       # compiler settings, see `optim_args` below
  seeds = 1,                          # more seeds = less noisy ratios, proportionally longer
  n_participants = 6,                 # participants drawn from rtdists::speed_acc
  n_trials = 150,                     # trials per participant, balanced across conditions
  chains = 4,
  cores = 4,
  iter = 1000,
  warmup = 500,
  out_dir = file.path("benchmarks", "results"),
  save_fits = FALSE                   # TRUE writes each brmsfit to out_dir (large)
)
settings <- if (exists("bench_overrides")) modifyList(defaults, bench_overrides) else defaults
list2env(settings, envir = environment())
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# Passed to cmdstanr::cmdstan_model() through brm(stan_model_args = ). "O1" is
# a stanc flag; the other two are CmdStan makefile flags. STAN_NO_RANGE_CHECKS
# drops the bounds checks on array/vector indexing, which is safe once the
# program is known to run without indexing errors, as these are.
optim_args <- list(
  default = list(),
  all = list(
    stanc_options = list("O1"),
    cpp_options = list(STAN_CPP_OPTIMS = TRUE, STAN_NO_RANGE_CHECKS = TRUE)
  )
)

df <- bench_data(n_participants, n_trials)

# ---- Measurement ----------------------------------------------------------------

chain_times <- function(fit) {
  # Stan's own per-chain timings (columns warmup, sample), which survive the CSV
  # to stanfit conversion brms does for cmdstanr. NULL if unavailable.
  t <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)
  if (is.null(t) || !all(is.finite(t))) NULL else t
}

measure <- function(fit, fallback_wall) {
  vars <- variables(fit)
  vars <- vars[!vars %in% c("lp__", "lprior") & !startsWith(vars, "z_")]
  draws <- as_draws_array(fit, variable = vars)
  s <- summarise_draws(draws, "rhat", "ess_bulk", "ess_tail")

  # Population-level parameters: fixed effects and random-effect SDs. Their
  # pairwise correlations are what a dense metric can absorb and a diagonal
  # one cannot.
  pop <- vars[startsWith(vars, "b_") | startsWith(vars, "sd_")]
  cors <- cor(as_draws_matrix(fit, variable = pop))
  max_abs_cor <- max(abs(cors[upper.tri(cors)]))

  np <- nuts_params(fit)  # post-warmup only
  by_par <- function(p) np$Value[np$Parameter == p]
  n_leapfrog <- sum(by_par("n_leapfrog__"))

  # Chains run in parallel, so wall time is the slowest chain. Gradient
  # throughput is per core-second: all post-warmup leapfrog steps over the
  # summed sampling time of the chains.
  t <- chain_times(fit)
  wall_s <- if (is.null(t)) fallback_wall else max(rowSums(t))
  sample_s <- if (is.null(t)) NA_real_ else sum(t[, "sample"])
  ess_min <- min(s$ess_bulk, na.rm = TRUE)

  data.frame(
    wall_s = wall_s,
    ess_bulk_min = ess_min,
    ess_tail_min = min(s$ess_tail, na.rm = TRUE),
    ess_per_s = ess_min / wall_s,
    ess_per_grad = ess_min / n_leapfrog,
    grad_per_s = n_leapfrog / sample_s,
    leapfrog_mean = mean(by_par("n_leapfrog__")),
    treedepth_mean = mean(by_par("treedepth__")),
    stepsize_mean = mean(by_par("stepsize__")),
    divergent = sum(by_par("divergent__")),
    rhat_max = max(s$rhat, na.rm = TRUE),
    # Dimension of the space Stan samples in: the b_*_Intercept columns are
    # recomputed from the centred `Intercept*` parameters, not sampled twice.
    n_pars = sum(!grepl("^b_(.*_)?Intercept$", vars)),
    max_abs_cor = max_abs_cor
  )
}

# ---- Run ------------------------------------------------------------------------

results_file <- file.path(out_dir, "metric_dense_e.csv")
results <- list()

for (family in families) {
  for (structure in structures) {
    tag <- paste(family, structure, sep = "_")
    f <- make_formula(family, structure)

    for (optim in optims) {
      cat(sprintf("\n==== %s | optim = %s: compiling ====\n", tag, optim))
      # One executable directory per compiler setting, so that cmdstanr never
      # mistakes a program compiled under one setting for the other.
      exe_dir <- file.path(normalizePath(out_dir), "exe", optim)
      dir.create(exe_dir, showWarnings = FALSE, recursive = TRUE)

      # chains = 0 compiles without a real run (brms still takes a single
      # transition to build the object, whose diagnostics are meaningless;
      # hence the suppress*). Every fit below reuses this executable.
      compile_time <- system.time(
        compiled <- tryCatch(
          suppressMessages(suppressWarnings(
            brm(f, data = df,
                prior = cogmod_priors(f, df),
                init = cogmod_inits(f, df),
                stanvars = cogmod_stanvars(f),
                backend = "cmdstanr", chains = 0, silent = 2,
                stan_model_args = c(optim_args[[optim]], list(dir = exe_dir)))
          )),
          error = function(e) { message("compile failed: ", conditionMessage(e)); NULL }
        )
      )
      if (is.null(compiled)) next
      cat(sprintf("     compiled in %.0f s\n", compile_time["elapsed"]))

      for (seed in seeds) {
        for (metric in metrics) {
          cat(sprintf("---- %s | optim = %s | metric = %s | seed = %d\n", tag, optim, metric, seed))
          timing <- system.time(
            fit <- tryCatch(
              brm(fit = compiled,
                  init = cogmod_inits(f, df),
                  chains = chains, cores = cores, iter = iter, warmup = warmup,
                  seed = seed, metric = metric, silent = 2, refresh = 0),
              error = function(e) { message("sampling failed: ", conditionMessage(e)); NULL }
            )
          )
          if (is.null(fit)) next

          row <- cbind(
            data.frame(family = family, structure = structure, optim = optim,
                       metric = metric, seed = seed, n_obs = nrow(df),
                       n_participants = n_participants, chains = chains,
                       iter = iter, warmup = warmup,
                       compile_s = unname(compile_time["elapsed"])),
            measure(fit, unname(timing["elapsed"]))
          )
          print(row[, c("wall_s", "ess_bulk_min", "ess_per_s", "ess_per_grad", "grad_per_s",
                        "leapfrog_mean", "divergent", "rhat_max", "max_abs_cor")],
                row.names = FALSE, digits = 3)
          results[[length(results) + 1]] <- row
          write.csv(do.call(rbind, results), results_file, row.names = FALSE)

          if (isTRUE(save_fits)) {
            saveRDS(fit, file.path(out_dir, sprintf("%s_%s_%s_seed%d.rds", tag, optim, metric, seed)))
          }
          rm(fit); gc(verbose = FALSE)
        }
      }
    }
  }
}

# ---- Summary ------------------------------------------------------------------------

results <- do.call(rbind, results)
has_baseline <- !is.null(results) &&
  any(results$metric == "diag_e" & results$optim == "default")
if (!has_baseline) {
  cat("\nNothing to compare (need the diag_e / default baseline).\n")
} else {
  # Mean over seeds, then every (metric, optim) cell relative to the baseline
  # cell (diag_e, default) of the same family x structure.
  agg <- aggregate(
    cbind(compile_s, wall_s, ess_per_s, ess_per_grad, grad_per_s, leapfrog_mean, divergent) ~
      family + structure + optim + metric,
    data = results, FUN = mean
  )
  base <- agg[agg$metric == "diag_e" & agg$optim == "default", ]
  cmp <- merge(agg, base, by = c("family", "structure"), suffixes = c("", "_base"))
  first <- match(paste(cmp$family, cmp$structure, "diag_e", "default"),
                 paste(results$family, results$structure, results$metric, results$optim))
  summary_tbl <- data.frame(
    family = cmp$family,
    structure = cmp$structure,
    metric = cmp$metric,
    optim = cmp$optim,
    n_pars = results$n_pars[first],
    max_abs_cor = round(results$max_abs_cor[first], 2),
    compile_s = round(cmp$compile_s),
    wall_s = round(cmp$wall_s),
    wall_ratio = round(cmp$wall_s / cmp$wall_s_base, 2),
    grad_per_s_ratio = round(cmp$grad_per_s / cmp$grad_per_s_base, 2),
    leapfrog_ratio = round(cmp$leapfrog_mean / cmp$leapfrog_mean_base, 2),
    ess_per_grad_ratio = round(cmp$ess_per_grad / cmp$ess_per_grad_base, 2),
    ess_per_s_ratio = round(cmp$ess_per_s / cmp$ess_per_s_base, 2),
    divergent = cmp$divergent
  )
  summary_tbl <- summary_tbl[order(summary_tbl$family, summary_tbl$structure,
                                   summary_tbl$metric != "diag_e", summary_tbl$optim != "default"), ]
  cat("\n==== relative to the baseline (metric = diag_e, optim = default) ====\n")
  cat("ess_per_s_ratio > 1: more effective draws per second (the number that matters).\n",
      "grad_per_s_ratio > 1: each gradient evaluation got cheaper (what the compiler\n",
      "flags do). leapfrog_ratio < 1 / ess_per_grad_ratio > 1: fewer gradient\n",
      "evaluations per effective draw, i.e. the metric fits the geometry better.\n\n", sep = "")
  print(summary_tbl, row.names = FALSE)
  write.csv(summary_tbl, file.path(out_dir, "metric_dense_e_summary.csv"), row.names = FALSE)

  if (requireNamespace("ggplot2", quietly = TRUE)) {
    library(ggplot2)
    p <- ggplot(results, aes(x = family, y = ess_per_s, fill = metric)) +
      geom_col(position = position_dodge(width = 0.7), width = 0.65) +
      facet_grid(structure ~ optim, labeller = label_both) +
      scale_fill_manual(values = c(diag_e = "#3F51B5", dense_e = "#F4511E")) +
      labs(x = NULL, y = "min bulk ESS per second (higher is better)", fill = "metric") +
      theme_minimal()
    ggsave(file.path(out_dir, "metric_dense_e.png"), p, width = 8, height = 6, dpi = 150)
  }
  cat(sprintf("\nResults written to %s\n", normalizePath(out_dir)))
}
