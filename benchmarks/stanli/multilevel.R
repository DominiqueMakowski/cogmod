# stanli vs CmdStan on multilevel decision-making models
# ======================================================
#
# Question (2026-10-12): bench.R measured one model, the vignette's LNR with
# population-level effects only, on 4620 trials. What speed-up can one expect
# on the models people actually fit to choice + RT data: several families,
# participants as random effects, and few trials each?
#
# Seven models, each fitted to data simulated from its own RNG - outlier-free
# (poutlier = 0), 10 participants x 2 conditions (Accuracy, Speed) x 50
# trials, every participant with their own drift, boundary and non-decision
# time - so that every model is correctly specified:
#
#   lnr    cogmod_lnr(), sigmabias = 0 (the vignette's LNR)
#   lnr6   cogmod_lnr(), sigmabias estimated (all six parameters)
#   ddm4   cogmod_ddm(), sigmadrift = sigmabias = sigmandt = 0
#   ddm5   cogmod_ddm(), sigmadrift estimated
#   ddm7   cogmod_ddm(), sigmadrift, sigmabias and sigmandt estimated
#   lba    cogmod_lba2(), sigmazero = 1
#   rdm    cogmod_rdm()
#
# In each, the parameter that carries the choice gets `~ Condition +
# (Condition | Participant)`, the threshold (or the second accumulator) and
# `ndt` a random intercept, the rest are population-level intercepts: the
# vignette's formulas with participants added. Population values are close
# to the vignette's fits, with more errors (Accuracy ~5-10%, Speed ~15-25%),
# so that the error side is identified.
#
# Run from the package root:
#
#   Rscript benchmarks/stanli/multilevel.R emit --model lnr [--out DIR] [--trials 50]
#   Rscript benchmarks/stanli/multilevel.R fit  --model lnr --seed 1 [--out DIR] [--iter 500] [--tag NAME]
#   Rscript benchmarks/stanli/multilevel.R grad --model lnr --seed 1 [--out DIR] [--draws 40]
#                                              [--reps 15] [--tag NAME]
#
#   emit  simulate the data (fixed seed, so every `fit` of a model sees the
#         same data), write the brms program and inputs, check that stanli
#         builds it
#   fit   4 chains x (ITER warmup + ITER sampling) from the same
#         cogmod_inits() and seed, once through cmdstanr and once through
#         stanli, arms in an order that alternates with the seed; compile and
#         build timed apart from sampling. Keeps the cmdstanr CSVs for `grad`
#   grad  at DRAWS posterior draws of that seed's cmdstanr fit, CmdStan's
#         and stanli's log density and gradient compared, then the cost of
#         one gradient, single-threaded, in REPS alternating blocks
#
# Writes benchmarks/stanli/results/ml_{fit,chains,pars,check,time}_<model>_<seed>[_<tag>].csv.
#
# Both arms of `fit` end the same way, so they pay the same overhead: the
# draws as a stanfit (brms's own CSV reader for CmdStan, as brm(backend =
# "cmdstanr") does; as_stanfit() for stanli) through brms:::rename_pars() on
# an empty brm(). brm(backend = "mock", mock_fit = ...) gives the same
# brmsfit (README, "0.20.0"). The cmdstanr arm samples with cmdstanr
# directly, as brm() does, so that its CSVs can be kept.
#
# `grad` is a process of its own because it has to be: on Linux, loading
# CmdStan's model methods (an Rcpp shared object) into a process in which
# rstan has read a fit segfaults in dyn.load() - the same crash AGENTS.md
# records for `_expose()` after an example's brm(). The first cluster smoke
# run hit it in every task (2026-10-12); Windows does not.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()

mode <- commandArgs(trailingOnly = TRUE)[1]
steps <- c("emit", "fit", "grad")
if (is.na(mode) || !mode %in% steps) {
  stop("usage: Rscript benchmarks/stanli/multilevel.R ", paste(steps, collapse = "|"),
       " --model NAME [--options]", call. = FALSE)
}
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})
opt <- gp_args(list(out = "benchmarks/results/stanli_ml", model = "", seed = 1L,
                    iter = 500L, draws = 40L, reps = 15L, tag = "", trials = 50L,
                    participants = 10L))
res_dir <- "benchmarks/stanli/results"
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)

gp_load(".")
suppressPackageStartupMessages(library(stanli))

# ---- Models ------------------------------------------------------------------

# Per model: the formula, and the truth the data are drawn from. `truth` is
# per dpar on its link scale: `b` = (intercept, Speed effect), `sd` = the
# participants' SD of (intercept, Speed effect), drawn independently; `fixed`
# a value on the natural scale, for a dpar fixed in the formula. `rng` maps
# the dpars, named as brms names them (`mu` is the first one), to the
# family's r<family>().
P <- "Participant"
.ML <- list(
  lnr = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      nuone ~ Condition + (1 | Participant), sigmazero ~ 1, sigmaone ~ 1,
      sigmabias = 0, ndt ~ 1 + (1 | Participant), family = cogmod_lnr()),
    truth = list(mu = list(b = c(1, 0.4), sd = c(0.3, 0.2)),
                 nuone = list(b = c(-0.7, 0.8), sd = 0.3),
                 sigmazero = list(b = -0.2), sigmaone = list(b = 0.25),
                 sigmabias = list(fixed = 0),
                 ndt = list(b = log(0.25), sd = 0.15)),
    rng = function(n, p) rcogmod_lnr(n, nuzero = p$mu, nuone = p$nuone, sigmazero = p$sigmazero,
                                     sigmaone = p$sigmaone, sigmabias = p$sigmabias, ndt = p$ndt)),
  lnr6 = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      nuone ~ Condition + (1 | Participant), sigmazero ~ 1, sigmaone ~ 1,
      sigmabias ~ 1, ndt ~ 1 + (1 | Participant), family = cogmod_lnr()),
    truth = list(mu = list(b = c(1, 0.4), sd = c(0.3, 0.2)),
                 nuone = list(b = c(-0.7, 0.8), sd = 0.3),
                 sigmazero = list(b = -0.2), sigmaone = list(b = 0.25),
                 sigmabias = list(b = -0.43),  # softplus: 0.5
                 ndt = list(b = log(0.25), sd = 0.15)),
    rng = function(n, p) rcogmod_lnr(n, nuzero = p$mu, nuone = p$nuone, sigmazero = p$sigmazero,
                                     sigmaone = p$sigmaone, sigmabias = p$sigmabias, ndt = p$ndt)),
  ddm4 = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      boundary ~ Condition + (1 | Participant), bias ~ 1, ndt ~ 1 + (1 | Participant),
      sigmadrift = 0, sigmabias = 0, sigmandt = 0, family = cogmod_ddm()),
    truth = list(mu = list(b = c(-1.4, 0), sd = c(0.4, 0.2)),
                 boundary = list(b = c(1.38, -0.6), sd = 0.2),  # softplus: 1.6, 1.15
                 bias = list(b = 0),
                 sigmadrift = list(fixed = 0), sigmabias = list(fixed = 0), sigmandt = list(fixed = 0),
                 ndt = list(b = log(0.3), sd = 0.15)),
    rng = function(n, p) rcogmod_ddm(n, drift = p$mu, boundary = p$boundary, bias = p$bias,
                                     ndt = p$ndt, sigmadrift = p$sigmadrift,
                                     sigmabias = p$sigmabias, sigmandt = p$sigmandt)),
  ddm5 = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      boundary ~ Condition + (1 | Participant), bias ~ 1, ndt ~ 1 + (1 | Participant),
      sigmadrift ~ 1, sigmabias = 0, sigmandt = 0, family = cogmod_ddm()),
    truth = list(mu = list(b = c(-1.4, 0), sd = c(0.4, 0.2)),
                 boundary = list(b = c(1.38, -0.6), sd = 0.2),
                 bias = list(b = 0),
                 sigmadrift = list(b = 0.25),  # softplus: 0.83
                 sigmabias = list(fixed = 0), sigmandt = list(fixed = 0),
                 ndt = list(b = log(0.3), sd = 0.15)),
    rng = function(n, p) rcogmod_ddm(n, drift = p$mu, boundary = p$boundary, bias = p$bias,
                                     ndt = p$ndt, sigmadrift = p$sigmadrift,
                                     sigmabias = p$sigmabias, sigmandt = p$sigmandt)),
  ddm7 = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      boundary ~ Condition + (1 | Participant), bias ~ 1, ndt ~ 1 + (1 | Participant),
      sigmadrift ~ 1, sigmabias ~ 1, sigmandt ~ 1, family = cogmod_ddm()),
    truth = list(mu = list(b = c(-1.4, 0), sd = c(0.4, 0.2)),
                 boundary = list(b = c(1.38, -0.6), sd = 0.2),
                 bias = list(b = 0),
                 sigmadrift = list(b = 0.25),
                 sigmabias = list(b = qlogis(0.3)),
                 sigmandt = list(b = log(0.1)),
                 ndt = list(b = log(0.3), sd = 0.15)),
    rng = function(n, p) rcogmod_ddm(n, drift = p$mu, boundary = p$boundary, bias = p$bias,
                                     ndt = p$ndt, sigmadrift = p$sigmadrift,
                                     sigmabias = p$sigmabias, sigmandt = p$sigmandt)),
  lba = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      driftone ~ Condition + (1 | Participant), sigmazero = 1, sigmaone ~ 1,
      sigmabias ~ 1, boundary ~ 1 + (1 | Participant), ndt ~ 1 + (1 | Participant),
      family = cogmod_lba2()),
    truth = list(mu = list(b = c(3, 0.3), sd = c(0.4, 0.2)),
                 driftone = list(b = c(0.8, 0.8), sd = 0.3),
                 sigmazero = list(fixed = 1), sigmaone = list(b = 0.54),  # softplus: 1
                 sigmabias = list(b = -0.43),                             # softplus: 0.5
                 boundary = list(b = -0.71, sd = 0.2),                    # softplus: 0.4
                 ndt = list(b = log(0.2), sd = 0.15)),
    rng = function(n, p) rcogmod_lba2(n, driftzero = p$mu, driftone = p$driftone,
                                      sigmazero = p$sigmazero, sigmaone = p$sigmaone,
                                      sigmabias = p$sigmabias, boundary = p$boundary, ndt = p$ndt)),
  rdm = list(
    formula = function() brms::bf(
      RT | dec(Error) ~ Condition + (Condition | Participant),
      driftone ~ Condition + (1 | Participant), sigmabias ~ 1,
      boundary ~ 1 + (1 | Participant), ndt ~ 1 + (1 | Participant), family = cogmod_rdm()),
    truth = list(mu = list(b = c(2.95, 0.3), sd = c(0.3, 0.2)),   # softplus: 3
                 driftone = list(b = c(-0.03, 0.6), sd = 0.3),    # softplus: 0.68
                 sigmabias = list(b = -1.05),                     # softplus: 0.3
                 boundary = list(b = 0.25, sd = 0.2),             # softplus: 0.83
                 ndt = list(b = log(0.25), sd = 0.15)),
    rng = function(n, p) rcogmod_rdm(n, vzero = p$mu, vone = p$driftone, boundary = p$boundary,
                                     bias = p$sigmabias, ndt = p$ndt))
)
if (!opt$model %in% names(.ML)) {
  stop("--model must be one of ", paste(names(.ML), collapse = ", "), call. = FALSE)
}
spec <- .ML[[opt$model]]
path <- function(...) file.path(opt$out, paste0(opt$model, ...))

linkinv <- function(link, eta) switch(link, identity = eta, log = exp(eta),
                                      softplus = log1p(exp(eta)), logit = stats::plogis(eta),
                                      stop("link ", link, call. = FALSE))

# The family's links, by dpar, off the formula's family object.
family_links <- function(f) {
  fam <- f$family
  setNames(vapply(fam$dpars, function(d) {
    if (d == "mu") fam$link else fam[[paste0("link_", d)]]
  }, character(1)), fam$dpars)
}

simulate <- function(spec, n_ppt, n_trials, seed = 2026) {
  set.seed(seed)
  d <- expand.grid(trial = seq_len(n_trials), Condition = c("Accuracy", "Speed"),
                   Participant = seq_len(n_ppt), stringsAsFactors = FALSE)
  speed <- as.numeric(d$Condition == "Speed")
  links <- family_links(spec$formula())
  pars <- lapply(setNames(names(spec$truth), names(spec$truth)), function(dp) {
    t <- spec$truth[[dp]]
    if (!is.null(t$fixed)) return(rep(t$fixed, nrow(d)))
    b <- c(t$b, 0)[1:2]
    s <- c(t$sd, 0, 0)[1:2]
    u0 <- stats::rnorm(n_ppt, 0, s[1])
    u1 <- stats::rnorm(n_ppt, 0, s[2])
    linkinv(links[[dp]], b[1] + u0[d$Participant] + (b[2] + u1[d$Participant]) * speed)
  })
  sim <- spec$rng(nrow(d), pars)
  data.frame(Participant = d$Participant, Condition = d$Condition, RT = sim$rt,
             Error = as.integer(sim$response))
}

# ---- emit --------------------------------------------------------------------
if (mode == "emit") {
  df <- simulate(spec, opt$participants, opt$trials)
  f <- spec$formula()
  prior <- suppressMessages(cogmod_priors(f, df))
  sv <- cogmod_stanvars(f)
  code <- as.character(brms::make_stancode(f, data = df, prior = prior, stanvars = sv))
  sdat <- lapply(unclass(brms::make_standata(f, data = df, prior = prior, stanvars = sv)), identity)
  writeLines(code, path(".stan"))
  saveRDS(list(df = df, f = f, prior = prior, sv = sv, sdat = sdat), path("_inputs.rds"))
  cat(opt$model, ": N =", nrow(df), "\n")
  print(stats::aggregate(cbind(Error, RT) ~ Condition, df, mean), digits = 3)
  print(round(stats::quantile(df$RT, c(0, 0.01, 0.1, 0.5, 0.9, 0.99, 1)), 3))
  t0 <- Sys.time()
  e <- tryCatch({ stanli_model(code = code, data = sdat); "built" },
                error = function(e) conditionMessage(e))
  cat("stanli:", e, "in", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")
}

# ---- fit ---------------------------------------------------------------------
key <- paste0(opt$model, "_", opt$seed, if (nzchar(opt$tag)) paste0("_", opt$tag))
res_file <- function(stem) file.path(res_dir, paste0("ml_", stem, "_", key, ".csv"))
handoff <- file.path(opt$out, paste0(key, "_cmdstanr.rds"))  # fit -> grad
if (mode != "emit") {
  cat(opt$model, "seed", opt$seed, "- stanli", as.character(utils::packageVersion("stanli")),
      "- CmdStan", as.character(cmdstanr::cmdstan_version()), "\n")
}
if (mode == "fit") {
  inp <- readRDS(path("_inputs.rds"))
  f <- inp$f; df <- inp$df; prior <- inp$prior; sv <- inp$sv
  n <- opt$iter
  timed <- function(expr) { t0 <- Sys.time(); force(expr); as.numeric(Sys.time() - t0, units = "secs") }

  set.seed(opt$seed)
  inits <- cogmod_inits(f, df)
  init_list <- lapply(1:4, function(i) inits(i))
  # Off the seed again before sampling (benchmarks/lnr_start/README.md, "Two
  # earlier arrays lost chains"): with processx < 3.9 two tasks of one seed
  # on one node, as two tagged arrays running at once can put there, would
  # draw the same chain-process tags and kill each other's chains. Both
  # engines take the seed explicitly, so the draws do not change.
  set.seed(NULL)

  # backend = "cmdstanr" so that stancode() is the code cmdstanr compiles.
  dummy <- brms::brm(f, data = df, prior = prior, stanvars = sv, backend = "cmdstanr",
                     empty = TRUE)
  code <- as.character(brms::stancode(dummy))
  sdat <- lapply(unclass(brms::standata(dummy)), identity)
  as_brmsfit <- function(stanfit) { d <- dummy; d$fit <- stanfit; brms:::rename_pars(d) }

  # What a user waits for before sampling starts: the C++ compile, against
  # stanli's parse and build. force_recompile, so that a cached executable
  # cannot make the first look free.
  build <- list()
  build$cmdstanr <- timed(cmod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(code, dir = opt$out), force_recompile = TRUE, quiet = TRUE))
  build$stanli <- timed(sm <- stanli_model(code = code, data = sdat))
  init_u <- do.call(rbind, lapply(init_list, unconstrain, model = sm))

  fc <- NULL
  arms <- list(
    cmdstanr = function() {
      fc <<- cmod$sample(data = sdat, chains = 4, parallel_chains = 4, iter_warmup = n,
                         iter_sampling = n, init = init_list, seed = opt$seed, refresh = 0,
                         show_messages = FALSE, output_dir = opt$out,
                         output_basename = paste0(key, "_cmdstanr"))
      as_brmsfit(brms:::read_csv_as_stanfit(fc$output_files(), model = cmod))
    },
    stanli = function() as_brmsfit(as_stanfit(sample_model(
      sm, chains = 4, seed = opt$seed, warmup = n, samples = n, init = init_u,
      parallel_chains = 4, refresh = 0))))
  arm_order <- if (opt$seed %% 2 == 1) names(arms) else rev(names(arms))

  fits <- list()
  chains <- list()
  for (a in arm_order) {
    cat("fit:", a, "\n")
    t0 <- Sys.time()
    fit <- arms[[a]]()
    wall <- as.numeric(Sys.time() - t0, units = "secs")
    el <- rstan::get_elapsed_time(fit$fit)
    np <- brms::nuts_params(fit)
    per <- function(p, fun = sum) as.numeric(tapply(np$Value[np$Parameter == p],
                                                    np$Chain[np$Parameter == p], fun))
    chains[[a]] <- data.frame(
      model = opt$model, seed = opt$seed, arm = a, position = match(a, arm_order),
      chain = seq_len(nrow(el)), warmup_s = el[, "warmup"], sample_s = el[, "sample"],
      leapfrog = per("n_leapfrog__"), divergent = per("divergent__"),
      treedepth_max = per("treedepth__", max), stepsize = per("stepsize__", mean),
      build_s = build[[a]], wall_s = wall)
    fits[[a]] <- fit
  }
  chains <- do.call(rbind, chains)
  saveRDS(fits, path("_fits_", opt$seed, if (nzchar(opt$tag)) paste0("_", opt$tag), ".rds"))

  # Every sampled quantity brms reports, fixed dpars (sd 0) left out. "pop"
  # is everything but the participants' own deviations (r_).
  pars <- do.call(rbind, lapply(names(fits), function(a) {
    dr <- posterior::subset_draws(brms::as_draws_df(fits[[a]]), regex = TRUE,
      variable = "^(b|sd|cor|r|poutlier|sigma|bias|boundary|ndt|drift|nu)")
    d <- as.data.frame(posterior::summarise_draws(dr, "mean", "sd", "mcse_mean", "rhat",
                                                  "ess_bulk", "ess_tail"))
    d <- d[is.finite(d$sd) & d$sd > 0, ]
    data.frame(model = opt$model, seed = opt$seed, arm = a,
               level = ifelse(grepl("^r_", d$variable), "r", "pop"), d)
  }))

  # The same posterior both ways: each mean against the other arm's, in
  # units of their joint Monte Carlo error.
  w <- merge(pars[pars$arm == "cmdstanr", c("variable", "mean", "mcse_mean")],
             pars[pars$arm == "stanli", c("variable", "mean", "mcse_mean")], by = "variable")
  z <- abs(w$mean.x - w$mean.y) / sqrt(w$mcse_mean.x^2 + w$mcse_mean.y^2)

  summ <- do.call(rbind, lapply(split(chains, chains$arm), function(x) {
    p <- pars[pars$arm == x$arm[1], ]
    pop <- p[p$level == "pop", ]
    chain_s <- x$warmup_s + x$sample_s
    # nuts_params() sees the sampling iterations only, so the cost of a
    # leapfrog step is the sampling time over them; dividing warmup + sampling
    # time by sampling leapfrogs (as this did until 2026-10-12) mixed in each
    # arm's luck in warmup.
    data.frame(model = opt$model, seed = opt$seed, arm = x$arm[1], position = x$position[1],
               N = nrow(df), n_par = ncol(init_u), build_s = x$build_s[1],
               chain_s_max = max(chain_s), wall_s = x$wall_s[1],
               leapfrog = sum(x$leapfrog), ms_per_leapfrog = 1000 * sum(x$sample_s) / sum(x$leapfrog),
               ess_bulk_min = min(pop$ess_bulk), ess_tail_min = min(pop$ess_tail),
               ess_bulk_min_r = min(p$ess_bulk[p$level == "r"], Inf),
               ess_bulk_per_wall_s = min(pop$ess_bulk) / x$wall_s[1],
               ess_bulk_per_1k_leapfrog = 1000 * min(pop$ess_bulk) / sum(x$leapfrog),
               rhat_max = max(p$rhat), divergent = sum(x$divergent),
               treedepth_max = max(x$treedepth_max), mean_z_max = max(z))
  }))
  summ <- summ[match(names(arms), summ$arm), ]
  summ$wall_vs_cmdstanr <- summ$wall_s / summ$wall_s[1]
  summ$ess_per_wall_s_vs_cmdstanr <- summ$ess_bulk_per_wall_s / summ$ess_bulk_per_wall_s[1]
  summ$ms_per_leapfrog_vs_cmdstanr <- summ$ms_per_leapfrog / summ$ms_per_leapfrog[1]
  utils::write.csv(chains, res_file("chains"), row.names = FALSE)
  utils::write.csv(pars, res_file("pars"), row.names = FALSE)
  utils::write.csv(summ, res_file("fit"), row.names = FALSE)
  print(summ, digits = 3, row.names = FALSE)
  saveRDS(list(csv = fc$output_files(), stan_file = cmod$stan_file(), sdat = sdat,
               init = init_list[[1]], n_par = ncol(init_u)), handoff)
}

# ---- grad --------------------------------------------------------------------
# DRAWS draws of the cmdstanr fit spread over its chains, unconstrained by
# CmdStan's model methods; both engines take the same vectors (the same
# program, so the same parameter order). Plus the first chain's init. The
# methods hang off a one-iteration fixed_param run of the same program,
# compiled again with them, as gradient_programs.R's gp_methods() does.
if (mode == "grad") {
  h <- readRDS(handoff)
  mod <- cmdstanr::cmdstan_model(h$stan_file, compile_model_methods = TRUE,
                                 force_recompile = TRUE, quiet = TRUE)
  fc <- mod$sample(data = h$sdat, init = list(h$init), chains = 1, iter_warmup = 1,
                   iter_sampling = 1, fixed_param = TRUE, refresh = 0, show_messages = FALSE,
                   seed = opt$seed)
  fc$init_model_methods(verbose = FALSE)
  U <- posterior::as_draws_matrix(fc$unconstrain_draws(files = h$csv, format = "draws_matrix"))
  rows <- round(seq(1, nrow(U), length.out = opt$draws))
  pts <- c(list(init = as.numeric(fc$unconstrain_variables(h$init))),
           setNames(lapply(rows, function(i) as.numeric(U[i, ])), paste0("draw", rows)))
  stopifnot(length(pts$init) == h$n_par, ncol(U) == h$n_par)
  # stanli single-threaded, as CmdStan's model methods are.
  sm1 <- stanli_model(code = paste(readLines(h$stan_file), collapse = "\n"), data = h$sdat,
                      threads_per_chain = 1)
  lpg <- list(
    cmdstan = function(q) { g <- fc$grad_log_prob(q); list(lp = attr(g, "log_prob"), grad = as.numeric(g)) },
    stanli = function(q) log_prob_grad(sm1, q))
  ulp <- function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_real_)
    if (!is.finite(a) || !is.finite(b)) return(if (identical(a, b)) 0 else Inf)
    if (a == b) return(0)
    abs(a - b) / 2^(floor(log2(max(abs(a), abs(b)))) - 52)
  }
  check <- do.call(rbind, lapply(names(pts), function(p) {
    r <- lapply(lpg, function(f) tryCatch(c(f(pts[[p]]), err = ""), error = function(e)
      list(lp = NA_real_, grad = rep(NA_real_, length(pts[[p]])), err = conditionMessage(e))))
    a <- r$cmdstan; b <- r$stanli
    both <- all(is.finite(a$grad)) && all(is.finite(b$grad))
    data.frame(model = opt$model, seed = opt$seed, point = p, lp_cmdstan = a$lp, lp_stanli = b$lp,
               lp_ulp = ulp(a$lp, b$lp), lp_absdiff = abs(a$lp - b$lp),
               grad_finite_cmdstan = all(is.finite(a$grad)), grad_finite_stanli = all(is.finite(b$grad)),
               grad_relerr = if (both) max(abs(b$grad - a$grad) / pmax(1, abs(a$grad))) else NA_real_,
               err = substr(paste0(a$err, b$err), 1, 120))
  }))
  utils::write.csv(check, res_file("check"), row.names = FALSE)
  cat(sprintf("\nexactness at %d points: lp identical %d, max %.3g ULP (abs %.3g); grad relerr max %.3g\n",
              nrow(check), sum(check$lp_ulp == 0, na.rm = TRUE), max(check$lp_ulp, na.rm = TRUE),
              max(check$lp_absdiff, na.rm = TRUE), max(check$grad_relerr, na.rm = TRUE)))

  # Cost: one block = one gradient at each of the draws (not the init, which
  # no sampler sits at), repeated `passes` times so that a block lasts about
  # half a second whatever the model costs (Windows' clock ticks in ~10 ms);
  # engines in random order within each rep.
  D <- pts[-1]
  timers <- list(cmdstan = function(k) for (i in seq_len(k)) for (q in D) fc$grad_log_prob(q),
                 stanli = function(k) for (i in seq_len(k)) for (q in D) log_prob_grad(sm1, q))
  t1 <- system.time(timers$cmdstan(1))[["elapsed"]]  # also the warm-up
  timers$stanli(1)
  passes <- max(1L, as.integer(ceiling(0.5 / max(t1, 0.005))))
  tm <- do.call(rbind, lapply(seq_len(opt$reps), function(r) {
    do.call(rbind, lapply(sample(names(timers)), function(nm) {
      st <- system.time(timers[[nm]](passes))
      data.frame(model = opt$model, seed = opt$seed, rep = r, engine = nm, passes = passes,
                 us_per_grad = st[["elapsed"]] / (passes * length(D)) * 1e6,
                 cpu_over_elapsed = (st[["user.self"]] + st[["sys.self"]]) / st[["elapsed"]])
    }))
  }))
  utils::write.csv(tm, res_file("time"), row.names = FALSE)
  w <- stats::reshape(tm[c("rep", "engine", "us_per_grad")], idvar = "rep", timevar = "engine",
                      direction = "wide")
  rel <- w$us_per_grad.stanli / w$us_per_grad.cmdstan
  cat(sprintf("per gradient: CmdStan %.0f us, stanli %.0f us, ratio %.2f [%.2f, %.2f]\n",
              stats::median(w$us_per_grad.cmdstan), stats::median(w$us_per_grad.stanli),
              stats::median(rel), stats::quantile(rel, 0.1), stats::quantile(rel, 0.9)))
}
