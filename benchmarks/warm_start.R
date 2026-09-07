# Benchmark: what does a warm start buy?
#
# Companion to the "Improving Sampling Efficiency and Performance" article
# (vignettes/articles/performance.qmd), whose findings on approximations and
# warm starts come from this script. Not part of the R package; run it from the
# repository root, after metric_dense_e.R or on its own:
#
#   Rscript benchmarks/warm_start.R
#
# Two demonstrations, on the same data as metric_dense_e.R:
#
# A. Reusing the adapted metric and step size of a previous fit. brms keeps
#    them in attr(fit$fit, "metadata"); a refit of the same model that starts
#    from them can run a much shorter warmup. Model: LNR, population-level only.
#    Runs: reference (warmup 500), warm restart (warmup 100, stored metric and
#    step size), and a cold run with the same short warmup, to show what the
#    stored quantities are worth.
#
# B. Initializing MCMC from Pathfinder, at the cmdstanr level: brms writes the
#    program and the data, Pathfinder supplies the initial values and, from
#    its draws mapped to the unconstrained space, a starting dense metric.
#    Model: DDM with participant random intercepts. Runs: reference (cogmod
#    inits, diag_e, warmup 500), warm start with the metric left to adapt
#    (warmup 300), and warm start with the metric held fixed (warmup 100,
#    init_buffer = warmup so that no adaptation window runs). Laplace and
#    Pathfinder are also compared with MCMC on the population-level
#    parameters, and the warm-started run is wrapped back into a brmsfit to
#    check that post-processing works.
#
# Output: benchmarks/results/warm_start.csv (one row per run: wall time
# including the approximation where one was used, min bulk ESS, ESS per
# second, divergences, max Rhat) and warm_start_approx.csv (posterior means
# and SDs from MCMC, Pathfinder and Laplace).

source("benchmarks/helpers.R")

defaults <- list(
  n_participants = 6,
  n_trials = 150,
  chains = 4,
  cores = 4,
  out_dir = file.path("benchmarks", "results")
)
settings <- if (exists("bench_overrides")) modifyList(defaults, bench_overrides) else defaults
list2env(settings, envir = environment())
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

df <- bench_data(n_participants, n_trials)
runs <- list()

record <- function(part, label, wall_s, x, divergent, warmup, approx_s = 0) {
  er <- ess_rhat(x)
  row <- data.frame(part = part, run = label, warmup = warmup,
                    approx_s = round(approx_s, 1), wall_s = round(wall_s + approx_s, 1),
                    ess_bulk_min = round(unname(er["ess_bulk_min"])),
                    ess_per_s = round(unname(er["ess_bulk_min"]) / (wall_s + approx_s), 2),
                    divergent = divergent, rhat_max = round(unname(er["rhat_max"]), 3))
  print(row, row.names = FALSE)
  runs[[length(runs) + 1]] <<- row
  write.csv(do.call(rbind, runs), file.path(out_dir, "warm_start.csv"), row.names = FALSE)
}

# ---- A. Metric and step size from a previous fit --------------------------------------

cat("\n==== A. Reusing a previous fit's metric (LNR, population-level) ====\n")
f <- make_formula("lnr", "fixed")
brms_wall <- function(m) max(rowSums(rstan::get_elapsed_time(m$fit)))
brms_div <- function(m) sum(nuts_params(m, pars = "divergent__")$Value)

m_ref <- suppressMessages(
  brm(f, data = df, prior = cogmod_priors(f, df), init = cogmod_inits(f, df),
      stanvars = cogmod_stanvars(f), backend = "cmdstanr", chains = chains, cores = cores,
      iter = 1000, warmup = 500, seed = 1, refresh = 0, silent = 2)
)
record("A", "reference", brms_wall(m_ref), m_ref, brms_div(m_ref), warmup = 500)

md <- attr(m_ref$fit, "metadata")   # $inv_metric and $step_size, one entry per chain
m_warm <- suppressMessages(
  brm(fit = m_ref, init = cogmod_inits(f, df), chains = chains, cores = cores,
      iter = 600, warmup = 100, seed = 2, refresh = 0, silent = 2,
      inv_metric = md$inv_metric[[1]], step_size = md$step_size[[1]])
)
record("A", "warm restart", brms_wall(m_warm), m_warm, brms_div(m_warm), warmup = 100)

m_cold <- suppressMessages(
  brm(fit = m_ref, init = cogmod_inits(f, df), chains = chains, cores = cores,
      iter = 600, warmup = 100, seed = 2, refresh = 0, silent = 2)
)
record("A", "cold, short warmup", brms_wall(m_cold), m_cold, brms_div(m_cold), warmup = 100)

cat("\nPopulation-level estimates (should agree):\n")
print(rbind(reference = fixef(m_ref)[, 1], warm = fixef(m_warm)[, 1], cold = fixef(m_cold)[, 1]), digits = 3)
rm(m_ref, m_warm, m_cold); gc(verbose = FALSE)

# ---- B. Pathfinder as a warm start -----------------------------------------------------

cat("\n==== B. Pathfinder warm start (DDM, participant random intercepts) ====\n")
f <- make_formula("ddm", "mixed")
prior <- cogmod_priors(f, df)
stanvars <- cogmod_stanvars(f)
init <- cogmod_inits(f, df)

# 1. brms writes the Stan program and the data; cmdstanr compiles it, with the
#    model methods needed to map draws to the unconstrained space.
scode <- stancode(f, data = df, prior = prior, stanvars = stanvars, backend = "cmdstanr")
sdata <- standata(f, data = df, prior = prior, stanvars = stanvars)
mod <- cmdstan_model(write_stan_file(scode), compile_model_methods = TRUE)

# 2. Approximations
t_pf <- system.time(
  pf <- mod$pathfinder(data = sdata, init = init, num_paths = 8, draws = 1000, seed = 1, refresh = 0)
)["elapsed"]
t_lp <- system.time(
  lp <- tryCatch(mod$laplace(data = sdata, init = init, draws = 1000, seed = 1, refresh = 0),
                 error = function(e) { message("Laplace failed: ", conditionMessage(e)); NULL })
)["elapsed"]
cat(sprintf("Pathfinder: %.1f s   Laplace: %.1f s\n", t_pf, t_lp))

# 3. Starting metric: covariance of the Pathfinder draws on the unconstrained scale
inv_metric <- cov(pf$unconstrain_draws(format = "draws_matrix"))

cmd_wall <- function(fit) max(fit$time()$chains$total)
cmd_div <- function(fit) sum(fit$diagnostic_summary(quiet = TRUE)$num_divergent)
cmd_draws <- function(fit) {
  vars <- fit$metadata()$stan_variables
  fit$draws(variables = vars[!vars %in% c("lp__", "lprior") & !startsWith(vars, "z_")])
}

# 4. Runs
fit_ref <- mod$sample(data = sdata, init = init, chains = chains, parallel_chains = cores,
                      iter_warmup = 500, iter_sampling = 500, seed = 1, refresh = 0)
record("B", "reference", cmd_wall(fit_ref), cmd_draws(fit_ref), cmd_div(fit_ref), warmup = 500)

fit_warm <- mod$sample(data = sdata, init = pf, chains = chains, parallel_chains = cores,
                       iter_warmup = 300, iter_sampling = 500, seed = 1, refresh = 0,
                       metric = "dense_e", inv_metric = inv_metric)
record("B", "warm, metric adapted", cmd_wall(fit_warm), cmd_draws(fit_warm), cmd_div(fit_warm),
       warmup = 300, approx_s = t_pf)

fit_fixed <- mod$sample(data = sdata, init = pf, chains = chains, parallel_chains = cores,
                        iter_warmup = 100, iter_sampling = 500, seed = 1, refresh = 0,
                        metric = "dense_e", inv_metric = inv_metric,
                        init_buffer = 100, term_buffer = 0, window = 0)
record("B", "warm, metric fixed", cmd_wall(fit_fixed), cmd_draws(fit_fixed), cmd_div(fit_fixed),
       warmup = 100, approx_s = t_pf)

# 5. Approximations against MCMC on the population-level parameters
vars <- grep("^(b|Intercept|sd_)", fit_ref$metadata()$stan_variables, value = TRUE)
approx <- data.frame(
  variable = fit_ref$summary(vars)$variable,
  mcmc_mean = fit_ref$summary(vars)$mean, mcmc_sd = fit_ref$summary(vars)$sd,
  pathfinder_mean = pf$summary(vars)$mean, pathfinder_sd = pf$summary(vars)$sd
)
if (!is.null(lp)) {
  approx$laplace_mean <- lp$summary(vars)$mean
  approx$laplace_sd <- lp$summary(vars)$sd
}
cat("\nApproximations vs MCMC (Pathfinder Pareto k is printed above by cmdstanr):\n")
print(approx, digits = 3, row.names = FALSE)
write.csv(approx, file.path(out_dir, "warm_start_approx.csv"), row.names = FALSE)

# 6. Back into brms, to check that post-processing works on the wrapped fit
m <- brm(f, data = df, prior = prior, stanvars = stanvars, backend = "cmdstanr", empty = TRUE)
m$fit <- read_csv_as_stanfit(fit_warm$output_files(),
                             variables = fit_warm$metadata()$stan_variables, model = mod)
m <- rename_pars(m)
cat("\nWrapped brmsfit: fixef() rows =", nrow(fixef(m)),
    "; posterior_predict() dims =", paste(dim(posterior_predict(m, ndraws = 5)), collapse = " x "), "\n")

cat("\n==== Summary ====\n")
print(do.call(rbind, runs), row.names = FALSE)
cat(sprintf("\nResults written to %s\n", normalizePath(out_dir)))
