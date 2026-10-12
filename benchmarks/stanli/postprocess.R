# Which post-processing tools work on a stanli fit taken into brms through
# brms's "mock" backend (the route in the performance article)?
#
#   Rscript benchmarks/stanli/postprocess.R [--out DIR] [--tag NAME]
#
# Run after `bench.R emit` (reads DIR/inputs.rds, default
# benchmarks/results/stanli). Fits the benchmark's LNR and a reduced LNR
# (`nuone ~ 1`, for loo_compare()) through stanli, the same two through
# brm(backend = "cmdstanr") as the control, and runs each tool on both, so
# that a failure can be put on the backend or on the family. Then an RT-only
# lognormal through stanli, for the tools that need posterior_epred(), which
# the LNR does not have by design. 4 chains x (300 + 300); about 15 minutes
# on Windows, most of it the cmdstanr compilations and the two rstan
# recompilations brms makes for moment matching and bridge sampling. Writes
# results/postprocess[_TAG].csv (tool, stanli, cmdstanr).

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()
opt <- gp_args(list(out = "benchmarks/results/stanli", tag = ""))
gp_load(".")
suppressPackageStartupMessages({ library(brms); library(stanli) })
options(mc.cores = 4)
inp <- readRDS(file.path(opt$out, "inputs.rds"))
df <- inp$df
f1 <- inp$f
f0 <- brms::bf(RT | dec(Error) ~ Condition, nuone ~ 1, sigmazero ~ 1,
               sigmaone ~ 1, sigmabias = 0, ndt ~ Condition, family = cogmod_lnr())
N_IT <- 300

# The performance article's four steps. save_pars(all = TRUE) because
# bridge_sampler() asks for it.
fit_stanli <- function(f) {
  prior <- cogmod_priors(f, df); stanvars <- cogmod_stanvars(f)
  scode <- stancode(f, data = df, prior = prior, stanvars = stanvars)
  sdata <- standata(f, data = df, prior = prior, stanvars = stanvars)
  model <- stanli_model(code = scode, data = sdata)
  inits <- cogmod_inits(f, df)
  init <- do.call(rbind, lapply(1:4, function(chain) unconstrain(model, inits(chain))))
  fit <- sample_model(model, chains = 4, parallel_chains = 4, warmup = N_IT, samples = N_IT,
                      init = init, refresh = 0, seed = 1)
  brm(f, data = df, prior = prior, stanvars = stanvars, save_pars = save_pars(all = TRUE),
      backend = "mock", mock_fit = as_stanfit(fit))
}
fit_cmdstanr <- function(f) brm(
  f, data = df, prior = cogmod_priors(f, df), stanvars = cogmod_stanvars(f),
  init = cogmod_inits(f, df), save_pars = save_pars(all = TRUE), backend = "cmdstanr",
  chains = 4, cores = 4, iter = 2 * N_IT, warmup = N_IT, refresh = 0, silent = 2, seed = 1)
set.seed(1)
fits <- list(stanli = list(fit_stanli(f1), fit_stanli(f0)),
             cmdstanr = list(fit_cmdstanr(f1), fit_cmdstanr(f0)))

# One line per result: a shape, a few numbers, or the error.
short <- function(x) {
  if (inherits(x, "ggplot")) { ggplot2::ggplot_build(x); return("ggplot built") }
  if (inherits(x, "brms_conditional_effects")) return(paste("effects:", paste(names(x), collapse = ",")))
  if (is.matrix(x) || is.array(x)) return(paste("dim", paste(dim(x), collapse = "x")))
  if (is.data.frame(x)) return(paste("rows", nrow(x)))
  substr(paste(format(signif(unlist(x), 5)), collapse = " "), 1, 60)
}
attempt <- function(fun, ...) tryCatch(
  withCallingHandlers(short(fun(...)), warning = function(w) invokeRestart("muffleWarning")),
  error = function(e) paste("ERROR:", gsub("\\s+", " ", substr(conditionMessage(e), 1, 100))))
rows <- list()
run <- function(tool, fun, engines = names(fits)) {
  r <- vapply(names(fits), function(e) if (e %in% engines) attempt(fun, fits[[e]][[1]], fits[[e]][[2]]) else "not run", "")
  cat(sprintf("\n## %s\n  stanli  : %s\n  cmdstanr: %s\n", tool, r[["stanli"]], r[["cmdstanr"]]))
  rows[[length(rows) + 1]] <<- data.frame(tool = tool, stanli = r[["stanli"]], cmdstanr = r[["cmdstanr"]])
}

run("fixef", function(m, m0) fixef(m)[, "Estimate"])
run("loo", function(m, m0) loo(m)$estimates["elpd_loo", ])
run("loo_compare (full vs reduced)", function(m, m0) loo_compare(loo(m), loo(m0))[2, 1:2])
run("waic", function(m, m0) waic(m)$estimates["elpd_waic", ])
run("loo, moment matching on the 3 worst Pareto k", function(m, m0) {
  k <- loo(m)$diagnostics$pareto_k
  loo(m, moment_match = TRUE, k_threshold = sort(k, decreasing = TRUE)[3] - 1e-9)$estimates["elpd_loo", ]
})
run("bridge_sampler (log marginal likelihood)", function(m, m0) bridgesampling::bridge_sampler(m, silent = TRUE)$logml)
run("log_lik", function(m, m0) log_lik(m, ndraws = 20))
run("posterior_predict", function(m, m0) posterior_predict(m, ndraws = 20))
run("posterior_linpred", function(m, m0) posterior_linpred(m, ndraws = 20))
run("posterior_epred", function(m, m0) posterior_epred(m, ndraws = 20))
run("conditional_effects", function(m, m0) conditional_effects(m))
run("pp_check", function(m, m0) pp_check(m, ndraws = 10))
run("hypothesis", function(m, m0) hypothesis(m, "ConditionSpeed > 0")$hypothesis$Post.Prob)
run("nuts_params, rhat, neff_ratio", function(m, m0) c(nrow(nuts_params(m)), max(rhat(m), na.rm = TRUE)))
run("priorsense::powerscale_sensitivity", function(m, m0) as.data.frame(priorsense::powerscale_sensitivity(m)))
run("bayestestR::describe_posterior", function(m, m0) bayestestR::describe_posterior(m, test = "pd"))
run("emmeans", function(m, m0) as.data.frame(emmeans::emmeans(m, ~ Condition)))
run("modelbased::estimate_means", function(m, m0) as.data.frame(modelbased::estimate_means(m, by = "Condition")))
run("modelbased::estimate_means, link scale", function(m, m0)
  as.data.frame(modelbased::estimate_means(m, by = "Condition", predict = "link"))$Median)
run("marginaleffects::avg_predictions, link scale", function(m, m0)
  marginaleffects::avg_predictions(m, by = "Condition", type = "link")$estimate)
run("saveRDS / readRDS, then loo", function(m, m0) {
  p <- tempfile(fileext = ".rds"); saveRDS(m, p); loo(readRDS(p))$estimates["elpd_loo", 1]
})
run("update(iter = 100)", function(m, m0) class(update(m, iter = 100, refresh = 0)), engines = "stanli")

# A family with posterior_epred(): the tools built on expectations.
lg <- fit_stanli(bf(RT ~ Condition, family = cogmod_lognormal()))
for (tool in c("posterior_epred", "conditional_effects", "pp_check", "bayes_R2",
               "modelbased::estimate_contrasts", "marginaleffects::avg_comparisons")) {
  fun <- switch(tool,
    posterior_epred = function(m) posterior_epred(m, ndraws = 20),
    conditional_effects = function(m) conditional_effects(m),
    pp_check = function(m) pp_check(m, ndraws = 10),
    bayes_R2 = function(m) bayes_R2(m)[, "Estimate"],
    `modelbased::estimate_contrasts` = function(m) as.data.frame(modelbased::estimate_contrasts(m, contrast = "Condition"))$Median,
    `marginaleffects::avg_comparisons` = function(m) marginaleffects::avg_comparisons(m, variables = "Condition")$estimate)
  r <- attempt(fun, lg)
  cat(sprintf("\n## lognormal: %s\n  stanli  : %s\n", tool, r))
  rows[[length(rows) + 1]] <- data.frame(tool = paste("lognormal:", tool), stanli = r, cmdstanr = "not run")
}
res <- do.call(rbind, rows)
utils::write.csv(res, paste0("benchmarks/stanli/results/postprocess", if (nzchar(opt$tag)) paste0("_", opt$tag), ".csv"),
                 row.names = FALSE)
