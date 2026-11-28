# cmdstanr, stanli and stanli + nutpie on multilevel decision models
# ==================================================================
#
# Question (2026-11-15): stanli 0.26.0 adds sample_nutpie()
# (seantalts/stanli#535): nutpie's sampler (nuts-rs, through nutpieR) driving
# stanli's own log density and gradient over its BridgeStan interface, with
# no C++ toolchain. stanli's NUTS already takes 0.53-0.75x cmdstanr's wall
# time on the multilevel race models, 0.97-1.11x on the DDMs
# (benchmarks/stanli/README.md, "0.21.0").
# Does nutpie's adaptation buy more ESS per second on top of stanli's
# gradient, and how far is each from cmdstanr?
#
# One task = one seed: simulate the data, then fit the same brms program from
# the same cogmod_inits() and seed in each arm, one arm after the other on the
# same node:
#
#   cmdstanr        brms's program compiled by CmdStan, sampled by cmdstanr
#   stanli          stanli::sample_model(), stanli's NUTS (Stan's algorithm)
#   nutpie          stanli::sample_nutpie(), diagonal mass matrix
#   nutpie_lowrank  stanli::sample_nutpie(adaptation = "low_rank"); not in
#                   the default arms
#
# The data: each model simulated from its own RNG without outliers, so it is
# correctly specified, with every participant's own value of the parameter
# that carries the choice, the second accumulator or the boundary, and the
# non-decision time - multilevel.R's models of the same name (`lnr`, `lba`,
# `rdm`, `ddm4`, `ddm5`): the vignette's formulas with participants added.
# 10 participants x 100 trials (50 per condition, Accuracy and Speed). The
# data seed is fixed (--data_seed), so every seed of a model fits the same
# data and only the samplers' randomness varies.
#
# Run from the package root:
#
#   Rscript benchmarks/nutpieR/bench.R --seed 1 [--model lnr|lba|rdm|ddm4|ddm5] [--participants 10]
#           [--trials 100] [--iter 500] [--chains 4] [--arms cmdstanr:stanli:nutpie]
#           [--nutpie_warmup same|default|N] [--data_seed 2026] [--tag NAME]
#           [--out benchmarks/results/nutpie]
#
# --arms is split at ":" or ",", since sbatch --export splits at commas.
# --iter is warmup and sampling per chain. nutpie takes the same warmup
# unless --nutpie_warmup says otherwise: "default" is nuts-rs's own (400 for
# the diagonal, 800 for low rank), an integer that many.
#
# Writes benchmarks/nutpieR/results/{fit,chains,pars}_<model>_<seed>[_<tag>].csv.
#
# Every arm ends the same way, so each pays the same overhead: the draws as
# a stanfit (brms's CSV reader for CmdStan, as brm(backend = "cmdstanr")
# does; stanli::as_stanfit() for the other two) through brms:::rename_pars()
# on an empty brm(), as multilevel.R. The wall clock runs from the call that
# starts sampling to the brmsfit; `build_s` is the compile (cmdstanr) or the
# build (stanli, shared by the nutpie arms) before it. sample_nutpie() hands
# the program to the runtime again inside its own call, so that cost is in
# its wall time.
#
# Warmup draws are kept in every arm (save_warmup), so that the leapfrog
# steps of warmup are counted too: nutpie and stanli's NUTS evaluate the same
# gradient, so leapfrog steps per effective draw compare their adaptations
# without a clock. Summaries and ESS are on the sampling draws only.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()

opt <- gp_args(list(out = "benchmarks/results/nutpie", model = "lnr", seed = 1L, iter = 500L,
                    chains = 4L, participants = 10L, trials = 100L, data_seed = 2026L,
                    arms = "cmdstanr:stanli:nutpie", nutpie_warmup = "same", tag = ""))
res_dir <- "benchmarks/nutpieR/results"
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
all_arms <- c("cmdstanr", "stanli", "nutpie", "nutpie_lowrank")
arms_wanted <- strsplit(opt$arms, "[,:]")[[1]]
if (!length(arms_wanted) || !all(arms_wanted %in% all_arms)) {
  stop("--arms must be some of ", paste(all_arms, collapse = ", "), call. = FALSE)
}
nutpie_warmup <- switch(opt$nutpie_warmup, same = opt$iter, default = NULL,
                        suppressWarnings(as.integer(opt$nutpie_warmup)))
if (!is.null(nutpie_warmup) && (is.na(nutpie_warmup) || nutpie_warmup < 1)) {
  stop("--nutpie_warmup must be same, default or a positive integer", call. = FALSE)
}
if (opt$trials %% 2 != 0) stop("--trials (per participant) must be even: two conditions", call. = FALSE)

gp_load(".")
suppressPackageStartupMessages(library(stanli))
if (any(startsWith(arms_wanted, "nutpie"))) {
  # sample_nutpie() checks this too, but only once the earlier arms have run.
  if (!exists("sample_nutpie", asNamespace("stanli")) || !requireNamespace("nutpieR", quietly = TRUE)) {
    stop("the nutpie arms need stanli >= 0.26.0 and nutpieR >= 1.9.0", call. = FALSE)
  }
}
pkg_version <- function(p) tryCatch(as.character(utils::packageVersion(p)), error = function(e) NA_character_)
versions <- c(brms = pkg_version("brms"), stanli = pkg_version("stanli"),
              nutpieR = pkg_version("nutpieR"), cmdstanr = pkg_version("cmdstanr"),
              cmdstan = as.character(cmdstanr::cmdstan_version()))
cat(opt$model, "seed", opt$seed, "-", paste(names(versions), versions, collapse = ", "), "\n")

# ---- Models ------------------------------------------------------------------

# Per model: the formula, and the truth the data are drawn from. `truth` is
# per dpar on its link scale: `b` = (intercept, Speed effect), `sd` = the
# participants' SD of (intercept, Speed effect), drawn independently; `fixed`
# a value on the natural scale, for a dpar fixed in the formula. `rng` maps
# the dpars, named as brms names them (`mu` is the first one), to the
# family's r<family>(). Copied from multilevel.R, which has two more models
# in the same form (LNR-6, DDM-7): another entry here is another model.
.MODELS <- list(
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
                                     bias = p$sigmabias, ndt = p$ndt)),
  # The DDM with its three variabilities at 0 (Stan's 4-argument
  # wiener_lpdf()), and with the drift's estimated (its 5-argument one).
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
                                     sigmabias = p$sigmabias, sigmandt = p$sigmandt))
)
if (!opt$model %in% names(.MODELS)) {
  stop("--model must be one of ", paste(names(.MODELS), collapse = ", "), call. = FALSE)
}
spec <- .MODELS[[opt$model]]

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

simulate <- function(spec, n_ppt, n_trials, seed) {
  set.seed(seed)
  d <- expand.grid(trial = seq_len(n_trials / 2), Condition = c("Accuracy", "Speed"),
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

# ---- Data, program, starting values -----------------------------------------

df <- simulate(spec, opt$participants, opt$trials, opt$data_seed)
f <- spec$formula()
prior <- suppressMessages(cogmod_priors(f, df))
sv <- cogmod_stanvars(f)
cat("N =", nrow(df), "\n")
print(stats::aggregate(cbind(Error, RT) ~ Condition, df, mean), digits = 3)

set.seed(opt$seed)
inits <- cogmod_inits(f, df)
init_list <- lapply(seq_len(opt$chains), function(i) inits(i))
# Off the seed again before sampling (benchmarks/lnr_start/README.md, "Two
# earlier arrays lost chains"): with processx < 3.9 two tasks of one seed on
# one node would draw the same chain-process tags and kill each other's
# chains. Every engine takes the seed explicitly, so the draws do not change.
set.seed(NULL)

# backend = "cmdstanr" so that stancode() is the code cmdstanr compiles; the
# stanli arms run the same text.
dummy <- brms::brm(f, data = df, prior = prior, stanvars = sv, backend = "cmdstanr",
                   empty = TRUE)
code <- as.character(brms::stancode(dummy))
sdat <- lapply(unclass(brms::standata(dummy)), identity)
as_brmsfit <- function(stanfit) { d <- dummy; d$fit <- stanfit; brms:::rename_pars(d) }
timed <- function(expr) { t0 <- Sys.time(); force(expr); as.numeric(Sys.time() - t0, units = "secs") }
key <- paste0(opt$model, "_", opt$seed, if (nzchar(opt$tag)) paste0("_", opt$tag))
res_file <- function(stem) file.path(res_dir, paste0(stem, "_", key, ".csv"))

# What a user waits for before sampling starts: the C++ compile, against
# stanli's parse and build. force_recompile, so that a cached executable
# cannot make the first look free.
build <- list()
if ("cmdstanr" %in% arms_wanted) {
  build$cmdstanr <- timed(cmod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(code, dir = opt$out), force_recompile = TRUE, quiet = TRUE))
}
build$stanli <- timed(sm <- stanli_model(code = code, data = sdat))
build$nutpie <- build$nutpie_lowrank <- build$stanli
init_u <- do.call(rbind, lapply(init_list, unconstrain, model = sm))

# ---- Arms --------------------------------------------------------------------

n <- opt$iter
nutpie_arm <- function(adaptation) function() {
  # nutpieR takes per-chain constrained inits and `cores` through `...`;
  # it does not promise which chain gets which start, only one each.
  as_brmsfit(as_stanfit(sample_nutpie(
    sm, chains = opt$chains, seed = opt$seed, warmup = nutpie_warmup, samples = n,
    adaptation = adaptation, save_warmup = TRUE, init = init_list, cores = opt$chains)))
}
arms <- list(
  cmdstanr = function() {
    fc <- cmod$sample(data = sdat, chains = opt$chains, parallel_chains = opt$chains,
                      iter_warmup = n, iter_sampling = n, save_warmup = TRUE, init = init_list,
                      seed = opt$seed, refresh = 0, show_messages = FALSE, output_dir = opt$out,
                      output_basename = paste0(key, "_cmdstanr"))
    as_brmsfit(brms:::read_csv_as_stanfit(fc$output_files(), model = cmod))
  },
  stanli = function() as_brmsfit(as_stanfit(sample_model(
    sm, chains = opt$chains, seed = opt$seed, warmup = n, samples = n, init = init_u,
    save_warmup = TRUE, parallel_chains = opt$chains, refresh = 0))),
  nutpie = nutpie_arm("diag"),
  nutpie_lowrank = nutpie_arm("low_rank"))[arms_wanted]
# The order rotates with the seed, so over any k consecutive seeds every arm
# runs once in each position.
k <- length(arms)
arm_order <- names(arms)[(seq_len(k) + opt$seed - 2) %% k + 1]
cat("arm order:", paste(arm_order, collapse = ", "), "\n")

fits <- list()
chains <- list()
for (a in arm_order) {
  cat("fit:", a, "\n")
  t0 <- Sys.time()
  fit <- arms[[a]]()
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  # nutpie times warmup and sampling together: its warmup column is NA and
  # its sampling column the whole run, the same for every chain.
  el <- rstan::get_elapsed_time(fit$fit)
  np <- brms::nuts_params(fit, inc_warmup = TRUE)
  n_warm <- fit$fit@sim$warmup2[1]
  per <- function(p, phase, fun = sum) {
    x <- np[np$Parameter == p & (np$Iteration > n_warm) == (phase == "sampling"), ]
    as.numeric(tapply(x$Value, x$Chain, fun))
  }
  chains[[a]] <- data.frame(
    model = opt$model, seed = opt$seed, arm = a, position = match(a, arm_order),
    chain = seq_len(nrow(el)), warmup_iter = n_warm, warmup_s = el[, "warmup"],
    sample_s = el[, "sample"], leapfrog_warmup = per("n_leapfrog__", "warmup"),
    leapfrog = per("n_leapfrog__", "sampling"), divergent = per("divergent__", "sampling"),
    treedepth_max = per("treedepth__", "sampling", max),
    stepsize = per("stepsize__", "sampling", mean), build_s = build[[a]], wall_s = wall)
  fits[[a]] <- fit
}
chains <- do.call(rbind, chains)

# ---- Summaries ---------------------------------------------------------------

# Every sampled quantity brms reports, fixed dpars (sd 0) left out. "pop" is
# everything but the participants' own deviations (r_). as_draws_df() leaves
# the warmup out. A correlated group term also comes through as the Stan
# matrix r_1[j, k], a copy of r_Participant[j, ...]: dropped, so that no
# deviation counts twice (the 2026-11-15 runs still have them).
pars <- do.call(rbind, lapply(names(fits), function(a) {
  dr <- posterior::subset_draws(brms::as_draws_df(fits[[a]]), regex = TRUE,
    variable = "^(b|sd|cor|r|poutlier|sigma|bias|boundary|ndt|drift|nu)")
  d <- as.data.frame(posterior::summarise_draws(dr, "mean", "sd", "mcse_mean", "rhat",
                                                "ess_bulk", "ess_tail"))
  d <- d[is.finite(d$sd) & d$sd > 0 & !grepl("^r_[0-9]+\\[", d$variable), ]
  data.frame(model = opt$model, seed = opt$seed, arm = a,
             level = ifelse(grepl("^r_", d$variable), "r", "pop"), d)
}))

# The same posterior every way: each mean against cmdstanr's, in units of
# their joint Monte Carlo error (against the first arm if cmdstanr is not run).
ref <- if ("cmdstanr" %in% names(fits)) "cmdstanr" else names(arms)[1]
mean_z <- function(a) {
  w <- merge(pars[pars$arm == ref, c("variable", "mean", "mcse_mean")],
             pars[pars$arm == a, c("variable", "mean", "mcse_mean")], by = "variable")
  max(abs(w$mean.x - w$mean.y) / sqrt(w$mcse_mean.x^2 + w$mcse_mean.y^2))
}

summ <- do.call(rbind, lapply(split(chains, chains$arm), function(x) {
  a <- x$arm[1]
  p <- pars[pars$arm == a, ]
  pop <- p[p$level == "pop", ]
  ess <- min(pop$ess_bulk)
  total_s <- x$build_s[1] + x$wall_s[1]
  data.frame(
    model = opt$model, seed = opt$seed, arm = a, position = x$position[1], N = nrow(df),
    n_par = ncol(init_u), chains = opt$chains, warmup_iter = x$warmup_iter[1], sampling_iter = n,
    build_s = x$build_s[1], wall_s = x$wall_s[1], total_s = total_s,
    # Sampling-phase cost of a leapfrog step; NA for nutpie, which reports no
    # phase times.
    ms_per_leapfrog = 1000 * sum(x$sample_s) / sum(x$leapfrog) * ifelse(anyNA(x$warmup_s), NA, 1),
    leapfrog_warmup = sum(x$leapfrog_warmup), leapfrog = sum(x$leapfrog),
    ess_bulk_min = ess, ess_tail_min = min(pop$ess_tail), ess_bulk_median = stats::median(pop$ess_bulk),
    ess_bulk_min_r = min(p$ess_bulk[p$level == "r"], Inf),
    ess_bulk_per_wall_s = ess / x$wall_s[1], ess_tail_per_wall_s = min(pop$ess_tail) / x$wall_s[1],
    ess_bulk_per_total_s = ess / total_s,
    ess_bulk_per_1k_leapfrog = 1000 * ess / sum(x$leapfrog),
    leapfrog_per_ess = (sum(x$leapfrog_warmup) + sum(x$leapfrog)) / ess,
    rhat_max = max(p$rhat), divergent = sum(x$divergent), treedepth_max = max(x$treedepth_max),
    stepsize = mean(x$stepsize), mean_z_max = if (a == ref) 0 else mean_z(a))
}))
summ <- summ[match(names(arms), summ$arm), ]
vs <- function(col) summ[[col]] / summ[[col]][summ$arm == ref]
summ$wall_vs_cmdstanr <- vs("wall_s")
summ$total_vs_cmdstanr <- vs("total_s")
summ$ess_per_wall_s_vs_cmdstanr <- vs("ess_bulk_per_wall_s")
summ$ess_tail_per_wall_s_vs_cmdstanr <- vs("ess_tail_per_wall_s")
summ$ess_per_total_s_vs_cmdstanr <- vs("ess_bulk_per_total_s")
summ$reference <- ref
for (v in names(versions)) summ[[paste0("v_", v)]] <- versions[[v]]
utils::write.csv(chains, res_file("chains"), row.names = FALSE)
utils::write.csv(pars, res_file("pars"), row.names = FALSE)
utils::write.csv(summ, res_file("fit"), row.names = FALSE)
print(summ[c("arm", "build_s", "wall_s", "ess_bulk_min", "ess_bulk_per_wall_s",
             "ess_per_wall_s_vs_cmdstanr", "ess_per_total_s_vs_cmdstanr", "leapfrog_per_ess",
             "rhat_max", "divergent", "mean_z_max")], digits = 3, row.names = FALSE)
