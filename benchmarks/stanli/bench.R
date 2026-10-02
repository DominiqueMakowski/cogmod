# stanli vs CmdStan on the LNR: does it run, and is it any faster?
# =================================================================
#
# Question (2026-10-02): paul-buerkner/brms#1911 proposes stanli - an
# interpreter for Stan programs, no C++ toolchain - as a brms backend, and
# shows a "light" route that works with brms as released: brm(empty = TRUE),
# sample its Stan code and data with stanli, as_stanfit(), rename_pars(). Does
# that work for a cogmod family, and how does it compare with brms's cmdstanr
# backend? README.md has the answer; this script reproduces it.
#
# Model: the decision_making vignette's LNR on its data (speed_acc,
# participants 1-3, RT <= 2 s; 4620 trials), as in benchmarks/lnr_vectorize.
# stanli refuses the program brms writes for it (`orig`) over one call,
# log_mix(), so it runs three variants, each checked against `orig`:
#
#   lse  `orig` with log_mix() written out as log_sum_exp(); the package's
#        density otherwise, branches and all
#   bf   a branch-free copy of the likelihood (lnr_branchfree_functions.stan,
#        sigmabias = 0 only) that drops the far-tail series of log Phi
#   sel  the same with the package's numerics kept, every branch written as
#        a select of clamped arms (lnr_select_functions.stan, sigmabias = 0)
#
# Run from the package root, in steps; each reads what the previous wrote.
#
#   Rscript benchmarks/stanli/bench.R emit [--out DIR]
#   Rscript benchmarks/stanli/bench.R grad [--out DIR] [--reps 21]
#   Rscript benchmarks/stanli/bench.R fit  [--out DIR]
#
#   emit  `orig`, the three variants, data and an init
#   grad  log density and gradient of `orig`, `bf` and `sel` in CmdStan and of the
#         three variants in stanli at 6 points (results/check.csv); the same
#         at a point and a response that put log Phi far into its tail
#         (results/tail.csv); then the cost of one gradient, alternating
#         blocks (results/time.csv)
#   fit   4 chains x 500 + 500: brm(backend = "cmdstanr") on `orig` against
#         the light route on `sel` (results/fit_chains.csv, fit_pars.csv)
#
# DIR (generated programs, executables, fits) defaults to
# benchmarks/results/stanli. Needs stanli (seantalts.r-universe.dev, then
# stanli::stanli_install() once), rtdists (the data), cmdstanr and CmdStan.
# `grad` compiles three programs with model methods, a minute or two each on
# Windows; `fit` compiles brms's own. Run nothing else on the machine while
# `grad` times or `fit` samples.
#
# How the cost is measured: R loops calling grad_log_prob() (cmdstanr model
# methods) and stanli's log_prob_grad() over 200 points, a different one each
# call, so nothing can be cached; programs alternate block by block, as in
# benchmarks/gradient_cost.R, and the headline is the ratio of medians within
# one run. Both go through R and a C++ boundary per call; at ~1-4 ms per
# gradient that overhead does not register. CPU time is recorded beside
# elapsed time, to show each uses one core.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile(), gp_methods()

mode <- commandArgs(trailingOnly = TRUE)[1]
steps <- c("emit", "grad", "fit")
if (is.na(mode) || !mode %in% steps) {
  stop("usage: Rscript benchmarks/stanli/bench.R ", paste(steps, collapse = "|"),
       " [--options]", call. = FALSE)
}
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})
opt <- gp_args(list(out = "benchmarks/results/stanli", reps = 21L))
out <- opt$out
res_dir <- "benchmarks/stanli/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
path <- function(...) file.path(out, paste0(...))

gp_load(".")
suppressPackageStartupMessages(library(stanli))

# ---- Data and model ----------------------------------------------------------

# benchmarks/lnr_vectorize/bench.R's lnr_data() and its `vignette` model.
lnr_data <- function() {
  data(speed_acc, package = "rtdists", envir = environment())
  df <- data.frame(
    Participant = as.integer(as.character(speed_acc$id)),
    Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[
      as.character(speed_acc$condition)]),
    RT = speed_acc$rt,
    Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat))
  )
  df[df$Participant %in% c(1, 2, 3) & df$RT <= 2, ]
}
lnr_formula <- function() {
  brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1,
           sigmaone ~ 1, sigmabias = 0, ndt ~ Condition, family = cogmod_lnr())
}

# The brms program with its functions block swapped for one of ours, and
# sigmabias dropped from the likelihood call. Parameters, priors and
# generated quantities are untouched, so rename_pars() and every brms
# post-processing function see the program they expect.
swap_functions <- function(code, file) {
  x <- strsplit(code, "\n")[[1]]
  f0 <- grep("^functions [{]", x)
  d0 <- grep("^data [{]", x)
  y <- c(x[seq_len(f0 - 1)], readLines(file), x[d0:length(x)])
  call0 <- "sigmaone[n], sigmabias, ndt[n]"
  if (sum(grepl(call0, y, fixed = TRUE)) != 1) stop("likelihood call not found", call. = FALSE)
  paste(sub(call0, "sigmaone[n], ndt[n]", y, fixed = TRUE), collapse = "\n")
}

# The brms program with log_mix() - the one function stanli 0.18.1 refuses
# inside a branch on a parameter - written out. Nothing else changes.
logsumexp <- function(code) {
  old <- "return log_mix(poutlier, lp_out, lp_dec);"
  if (!grepl(old, code, fixed = TRUE)) stop("log_mix() call not found", call. = FALSE)
  sub(old, "return log_sum_exp(log(poutlier) + lp_out, log1m(poutlier) + lp_dec);",
      code, fixed = TRUE)
}

.VARIANTS <- list(
  lse = logsumexp,
  bf = function(code) swap_functions(code, "benchmarks/stanli/lnr_branchfree_functions.stan"),
  sel = function(code) swap_functions(code, "benchmarks/stanli/lnr_select_functions.stan")
)

inputs <- function() readRDS(path("inputs.rds"))
code_of <- function(v) paste(readLines(path("lnr_", v, ".stan")), collapse = "\n")
stanli_of <- function(v, data) stanli_model(code = code_of(v), data = data, threads_per_chain = 1)
init_json <- function(init, file) {
  # write_stan_json() writes a length-1 R vector as a scalar, which CmdStan
  # rejects for the vector[1] coefficients (`b`, `b_nuone`, `b_ndt`) here; a
  # 1-d array keeps the brackets.
  vec <- grepl("^b(_|$)", names(init))
  init[vec] <- lapply(init[vec], function(v) array(v, dim = length(v)))
  cmdstanr::write_stan_json(init, file)
}
relerr <- function(g, ref) max(abs(g - ref) / pmax(1, abs(ref)))

# ---- emit --------------------------------------------------------------------
if (mode == "emit") {
  df <- lnr_data()
  f <- lnr_formula()
  prior <- suppressMessages(cogmod_priors(f, df))
  sv <- cogmod_stanvars(f)
  code <- as.character(brms::make_stancode(f, data = df, prior = prior, stanvars = sv))
  sdat <- lapply(unclass(brms::make_standata(f, data = df, prior = prior, stanvars = sv)), identity)
  writeLines(code, path("lnr_orig.stan"))
  for (v in names(.VARIANTS)) writeLines(.VARIANTS[[v]](code), path("lnr_", v, ".stan"))
  cmdstanr::write_stan_json(sdat, path("lnr.data.json"))
  init <- cogmod_inits(f, df, jitter = 0)(1)
  init_json(init, path("lnr.init.json"))
  saveRDS(list(df = df, f = f, prior = prior, sv = sv, sdat = sdat, init = init), path("inputs.rds"))
  cat("N =", sdat$N, "- wrote", out, "\n")
  e <- tryCatch({ stanli_model(code = code, data = sdat); "built" },
                error = function(e) conditionMessage(e))
  cat("stanli on the original program:", e, "\n")
}

# ---- grad --------------------------------------------------------------------
if (mode == "grad") {
  inp <- inputs()
  cat("compiling with model methods\n")
  mo <- gp_compile(path("lnr_orig.stan"))
  fo <- gp_methods(mo, path("lnr.data.json"), path("lnr.init.json"))
  fb <- gp_methods(gp_compile(path("lnr_bf.stan")), path("lnr.data.json"), path("lnr.init.json"))
  fsel <- gp_methods(gp_compile(path("lnr_sel.stan")), path("lnr.data.json"), path("lnr.init.json"))
  t0 <- Sys.time()
  ms <- lapply(setNames(names(.VARIANTS), names(.VARIANTS)), stanli_of, data = inp$sdat)
  cat("stanli build, three variants:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

  # Check: the init and five points around it. stanli's lp carries the
  # Jacobian and the constants, as CmdStan's log_prob(jacobian = TRUE) does.
  q0 <- fo$up0
  set.seed(1)
  Q <- rbind(q0, t(replicate(5, q0 + stats::rnorm(length(q0), 0, 0.15))))
  check <- do.call(rbind, lapply(seq_len(nrow(Q)), function(i) {
    q <- Q[i, ]
    go <- fo$fit$grad_log_prob(q)
    gb <- fb$fit$grad_log_prob(q)
    gs <- fsel$fit$grad_log_prob(q)
    s <- lapply(ms, log_prob_grad, q = q)
    data.frame(point = i, lp_orig = attr(go, "log_prob"), lp_cmdstan_bf = attr(gb, "log_prob"),
               lp_cmdstan_sel = attr(gs, "log_prob"), grad_cmdstan_sel = relerr(gs, go),
               lp_stanli_lse = s$lse$lp, lp_stanli_bf = s$bf$lp, lp_stanli_sel = s$sel$lp,
               grad_cmdstan_bf = relerr(gb, go), grad_stanli_lse = relerr(s$lse$grad, go),
               grad_stanli_bf = relerr(s$bf$grad, go), grad_stanli_sel = relerr(s$sel$grad, go))
  }))
  print(check, digits = 10)
  utils::write.csv(check, file.path(res_dir, "check.csv"), row.names = FALSE)

  # The far tail, which cogmod_log_Phi()'s x < -25 series is there for. Both
  # sigmas at softplus(-3) = 0.049 and the first trial replaced by a 12 s
  # response (gradient_check.R's slowest): the loser's survival then sits
  # at x ~ -60, past the x = -38 where the erfc route underflows.
  sd12 <- inp$sdat
  sd12$Y[1] <- 12
  cmdstanr::write_stan_json(sd12, path("lnr12.data.json"))
  init12 <- inp$init
  init12$Intercept_sigmazero <- init12$Intercept_sigmaone <- -3
  init_json(init12, path("lnr12.init.json"))
  fo12 <- gp_methods(mo, path("lnr12.data.json"), path("lnr12.init.json"))
  q12 <- fo12$up0
  go <- fo12$fit$grad_log_prob(q12)
  s_l <- log1p(exp(-3))
  nu_l <- if (sd12$dec[1] == 0) init12$Intercept_nuone else init12$Intercept
  cat(sprintf("tail check: loser's x on the 12 s trial = %.1f\n",
              (-nu_l - log(12 - exp(init12$Intercept_ndt))) / s_l))
  tail <- do.call(rbind, lapply(names(.VARIANTS), function(v) {
    g <- log_prob_grad(stanli_of(v, sd12), q12)
    data.frame(variant = v, lp_orig = attr(go, "log_prob"), lp = g$lp,
               grad_finite = all(is.finite(g$grad)),
               grad_relerr = if (all(is.finite(g$grad))) relerr(g$grad, go) else NA_real_)
  }))
  print(tail, digits = 10)
  utils::write.csv(tail, file.path(res_dir, "tail.csv"), row.names = FALSE)

  # Cost.
  K <- 200
  set.seed(3)
  P <- lapply(seq_len(K), function(k) q0 + stats::rnorm(length(q0), 0, 0.1))
  timers <- list(
    cmdstan_orig = function() for (k in seq_len(K)) fo$fit$grad_log_prob(P[[k]]),
    cmdstan_bf = function() for (k in seq_len(K)) fb$fit$grad_log_prob(P[[k]]),
    cmdstan_sel = function() for (k in seq_len(K)) fsel$fit$grad_log_prob(P[[k]]),
    stanli_lse = function() for (k in seq_len(K)) log_prob_grad(ms$lse, P[[k]]),
    stanli_bf = function() for (k in seq_len(K)) log_prob_grad(ms$bf, P[[k]]),
    stanli_sel = function() for (k in seq_len(K)) log_prob_grad(ms$sel, P[[k]])
  )
  for (nm in names(timers)) timers[[nm]]()  # warm up
  rows <- list()
  for (r in seq_len(opt$reps)) {
    for (nm in sample(names(timers))) {
      st <- system.time(timers[[nm]]())
      rows[[length(rows) + 1]] <- data.frame(rep = r, program = nm,
        us_per_grad = st[["elapsed"]] / K * 1e6,
        cpu_over_elapsed = (st[["user.self"]] + st[["sys.self"]]) / st[["elapsed"]])
    }
  }
  tm <- do.call(rbind, rows)
  utils::write.csv(tm, file.path(res_dir, "time.csv"), row.names = FALSE)
  med <- tapply(tm$us_per_grad, tm$program, stats::median)
  print(round(med))
  cat("ratio to cmdstan_orig:\n")
  print(round(med / med[["cmdstan_orig"]], 2))
  cat("CPU / elapsed (median):\n")
  print(round(tapply(tm$cpu_over_elapsed, tm$program, stats::median), 2))
}

# ---- fit ---------------------------------------------------------------------
if (mode == "fit") {
  inp <- inputs()
  f <- inp$f; df <- inp$df; prior <- inp$prior; sv <- inp$sv
  set.seed(11)
  inits <- cogmod_inits(f, df)
  init_list <- lapply(1:4, function(i) inits(i))

  cat("brm, cmdstanr backend\n")
  t0 <- Sys.time()
  fit_c <- brms::brm(f, data = df, prior = prior, stanvars = sv, backend = "cmdstanr",
                     chains = 4, cores = 4, warmup = 500, iter = 1000, init = init_list,
                     seed = 11, refresh = 0, silent = 2)
  wall_c <- as.numeric(Sys.time() - t0, units = "secs")

  # The light route from brms#1911, on `sel`. sample_model() takes no list of
  # per-chain inits in stanli 0.18.1 ("'list' object cannot be coerced to
  # type 'double'"), so its chains start from its own random inits (radius 2).
  cat("stanli, light route\n")
  t0 <- Sys.time()
  dummy <- brms::brm(f, data = df, prior = prior, stanvars = sv, empty = TRUE)
  m <- stanli_model(code = .VARIANTS$sel(as.character(brms::stancode(dummy))),
                    data = brms::standata(dummy))
  fs <- sample_model(m, chains = 4, seed = 11, warmup = 500, samples = 500,
                     parallel_chains = 4, refresh = 0)
  dummy$fit <- as_stanfit(fs)
  fit_s <- brms:::rename_pars(dummy)
  wall_s <- as.numeric(Sys.time() - t0, units = "secs")

  fits <- list(cmdstanr = fit_c, stanli = fit_s)
  saveRDS(fits, path("fits.rds"))
  chains <- do.call(rbind, lapply(names(fits), function(b) {
    fit <- fits[[b]]
    el <- rstan::get_elapsed_time(fit$fit)
    np <- brms::nuts_params(fit)
    per <- function(p) as.numeric(tapply(np$Value[np$Parameter == p], np$Chain[np$Parameter == p], sum))
    data.frame(backend = b, chain = seq_len(nrow(el)), warmup_s = el[, "warmup"],
               sample_s = el[, "sample"], leapfrog = per("n_leapfrog__"),
               divergent = per("divergent__"),
               wall_total_s = if (b == "cmdstanr") wall_c else wall_s)
  }))
  pars <- do.call(rbind, lapply(names(fits), function(b) {
    d <- posterior::summarise_draws(
      brms::as_draws_df(fits[[b]], variable = "^b_|^poutlier", regex = TRUE),
      "mean", "sd", "rhat", "ess_bulk", "ess_tail")
    data.frame(backend = b, as.data.frame(d))
  }))
  utils::write.csv(chains, file.path(res_dir, "fit_chains.csv"), row.names = FALSE)
  utils::write.csv(pars, file.path(res_dir, "fit_pars.csv"), row.names = FALSE)
  print(chains)
  print(pars, digits = 4)
  print(brms::loo_compare(brms::loo(fit_c), brms::loo(fit_s)))
}
