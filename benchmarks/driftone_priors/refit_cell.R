# One cell of the driftone prior check on Artemis (README.md): a vignette
# model refitted under arm A (the pushed tree's own priors), B (the
# proposal: R1-R6 of the 2026-10-12 prior review) or C (B with the tighter
# driftone fallback), 4 chains, the same seed - so the same starting values -
# in every arm.
#
# Usage: Rscript refit_cell.R <row of cells.csv>
# Env:   PC_DIR (the pushed tree), PC_WORK (this task's scratch directory)
idx <- as.integer(commandArgs(trailingOnly = TRUE)[1])
DIR <- Sys.getenv("PC_DIR")
WORK <- Sys.getenv("PC_WORK")
resdir <- file.path(DIR, "results")
dir.create(resdir, showWarnings = FALSE)
cells <- utils::read.csv(file.path(DIR, "cells.csv"), stringsAsFactors = FALSE)
name <- cells$model[idx]
arm <- cells$arm[idx]
seed <- cells$seed[idx]

suppressMessages(pkgload::load_all(file.path(DIR, "tree"), quiet = TRUE))
suppressMessages(library(brms))
log_msg <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), name, arm, seed, ..., "\n")
}
log_msg("brms", as.character(utils::packageVersion("brms")),
        "cmdstanr", as.character(utils::packageVersion("cmdstanr")),
        "cmdstan", cmdstanr::cmdstan_version())

ns <- asNamespace("cogmod")
patch <- function(nm, value) {
  if (bindingIsLocked(nm, ns)) unlockBinding(nm, ns)
  assign(nm, value, envir = ns)
}
if (arm %in% c("B", "C")) {
  ch <- get(".CHOICE", envir = ns)
  sh <- get(".SHIFTED", envir = ns)
  base <- get(".SHIFTED_BASE_PRIORS", envir = ns)
  dr <- switch(arm,
    B = c(link = "normal(-2, 3)", slope = "normal(0, 4)"),
    C = c(link = "normal(-2, 2)", slope = "normal(0, 3)")
  )
  ch$cogmod_rdm$prior$driftone <- c(link = dr[["link"]], nat = "lognormal(-2, 1.5)",
                                    slope = dr[["slope"]])
  ch$cogmod_lba2$prior$driftone <- c(link = dr[["link"]], nat = dr[["link"]],
                                     slope = dr[["slope"]])
  ch$cogmod_ddm$prior$boundary <- c(link = "normal(1, 1)", nat = "lognormal(0.27, 0.5)",
                                    slope = "normal(0, 0.5)")
  ch$cogmod_ddm$prior$bias <- c(link = "logistic(0, 1)", slope = "normal(0, 0.5)")
  for (fam in c("cogmod_invgaussian", "cogmod_exwald", "cogmod_bisa")) {
    if (is.null(sh[[fam]]$prior)) sh[[fam]]$prior <- list()
    sh[[fam]]$prior$boundary <- c(link = "normal(0.5, 1)", nat = "lognormal(0, 0.6)",
                                  slope = "normal(0, 0.5)")
  }
  base$poutlier <- c(link = "normal(-6, 1)", nat = "exponential(280)")
  patch(".CHOICE", ch)
  patch(".SHIFTED", sh)
  patch(".SHIFTED_BASE_PRIORS", base)
}

dm <- readRDS(file.path(DIR, "data", "dm.rds"))
sub <- readRDS(file.path(DIR, "data", "sub.rds"))
rt <- readRDS(file.path(DIR, "data", "rt.rds"))

f_rdm <- bf(RT | dec(Error) ~ Condition, driftone ~ Condition, sigmabias ~ 1,
            boundary ~ 1, ndt ~ 1, family = cogmod_rdm())
f_ddm <- function(bias, pout = FALSE) {
  fs <- list(RT | dec(Error) ~ Condition, boundary ~ Condition, bias, ndt ~ 1)
  if (pout) fs <- c(fs, list(poutlier ~ 1))
  do.call(bf, c(fs, list(sigmadrift = 0, sigmabias = 0, sigmandt = 0,
                         family = cogmod_ddm())))
}
spec <- list(
  rdm = list(data = dm, metric = "diag_e", f = f_rdm),
  rdm_sub = list(data = sub, metric = "diag_e", f = f_rdm),
  lba2 = list(data = dm, metric = "dense_e",
              f = bf(RT | dec(Error) ~ Condition, driftone ~ Condition, sigmazero = 1,
                     sigmaone ~ 1, sigmabias ~ 1, boundary ~ 1, ndt ~ 1,
                     family = cogmod_lba2())),
  ddm = list(data = dm, metric = "diag_e", f = f_ddm(bias ~ 1)),
  ddm_bias = list(data = dm, metric = "diag_e", f = f_ddm(bias ~ Condition)),
  ddm_pout = list(data = dm, metric = "diag_e", f = f_ddm(bias ~ 1, pout = TRUE)),
  wald4 = list(data = rt, metric = "diag_e",
               f = bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition,
                      sigmadrift ~ Condition, sigmandt = 0, family = cogmod_invgaussian())),
  exwald = list(data = rt, metric = "diag_e",
                f = bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition,
                       tau ~ Condition, family = cogmod_exwald())),
  bisa = list(data = rt, metric = "diag_e",
              f = bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition,
                     family = cogmod_bisa()))
)

softplus <- function(x) log1p(exp(x))
pp_race <- function(name, dd) {
  idx <- round(seq(1, nrow(dd), length.out = 100))
  out <- lapply(c("Accuracy", "Speed"), function(cond) {
    s <- as.numeric(cond == "Speed")
    r <- t(vapply(idx, function(i) {
      g <- function(v) dd[[v]][i]
      if (name %in% c("rdm", "rdm_sub")) {
        sim <- rcogmod_rdm(2000,
          vzero = softplus(g("b_Intercept") + s * g("b[1]")),
          vone = softplus(g("b_driftone_Intercept") + s * g("b_driftone[1]")),
          boundary = softplus(g("b_boundary_Intercept")),
          bias = softplus(g("b_sigmabias_Intercept")),
          ndt = exp(g("b_ndt_Intercept")), poutlier = 0)
      } else {
        sim <- rcogmod_lba2(2000,
          driftzero = g("b_Intercept") + s * g("b[1]"),
          driftone = g("b_driftone_Intercept") + s * g("b_driftone[1]"),
          sigmazero = 1, sigmaone = softplus(g("b_sigmaone_Intercept")),
          sigmabias = softplus(g("b_sigmabias_Intercept")),
          boundary = softplus(g("b_boundary_Intercept")),
          ndt = exp(g("b_ndt_Intercept")), poutlier = 0)
      }
      c(err = mean(sim$response == 1),
        rt_correct = stats::median(sim$rt[sim$response == 0]),
        rt_error = stats::median(sim$rt[sim$response == 1]))
    }, numeric(3)))
    data.frame(cond = cond, stat = colnames(r),
               median = apply(r, 2, stats::median, na.rm = TRUE),
               q5 = apply(r, 2, stats::quantile, 0.05, na.rm = TRUE),
               q95 = apply(r, 2, stats::quantile, 0.95, na.rm = TRUE))
  })
  do.call(rbind, out)
}

sp <- spec[[name]]
f <- sp$f
d <- sp$data
stem <- file.path(resdir, sprintf("%s__%s__s%d", name, arm, seed))
pri <- suppressMessages(suppressWarnings(cogmod_priors(f, d)))
pr <- as.data.frame(pri)
pr <- pr[nzchar(pr$prior) & (pr$dpar %in% c("driftone", "boundary", "bias", "poutlier") |
           pr$class %in% c("driftone", "boundary", "bias", "poutlier")),
         c("prior", "class", "coef", "group", "dpar")]
utils::write.csv(pr, paste0(stem, "__priors.csv"), row.names = FALSE)
log_msg("priors:", paste(sprintf("%s/%s/%s=%s", pr$dpar, pr$class, pr$coef, pr$prior),
                         collapse = "; "))
sv <- cogmod_stanvars(f)
code <- suppressWarnings(brms::make_stancode(f, data = d, prior = pri, stanvars = sv))
sdata <- suppressWarnings(brms::make_standata(f, data = d, prior = pri, stanvars = sv))
file <- cmdstanr::write_stan_file(code, dir = WORK)
log_msg("compiling")
mod <- cmdstanr::cmdstan_model(file, quiet = TRUE)
log_msg("sampling")
set.seed(seed)
init <- cogmod_inits(f, d)
fit <- mod$sample(
  data = unclass(sdata), seed = seed, chains = 4, parallel_chains = 4,
  iter_warmup = 250, iter_sampling = 250, init = init, refresh = 0,
  show_messages = FALSE, show_exceptions = FALSE, save_warmup = TRUE,
  metric = sp$metric, output_dir = WORK
)
sdg <- as.data.frame(fit$sampler_diagnostics(inc_warmup = TRUE, format = "draws_df"))
chains <- do.call(rbind, lapply(split(sdg, sdg$.chain), function(x) {
  w <- x$.iteration <= 250
  data.frame(chain = x$.chain[1],
             leap_warm = sum(x$n_leapfrog__[w]), leap_samp = sum(x$n_leapfrog__[!w]),
             div = sum(x$divergent__[!w]), div_warm = sum(x$divergent__[w]),
             td10 = sum(x$treedepth__[!w] >= 10),
             stepsize = x$stepsize__[!w][1])
}))
lp <- as.data.frame(fit$draws("lp__", format = "draws_df"))
chains$lp_mean <- as.numeric(tapply(lp$lp__, lp$.chain, mean))
chains$ebfmi <- fit$diagnostic_summary(quiet = TRUE)$ebfmi
tm <- fit$time()$chains
chains$time <- tm$total[match(chains$chain, tm$chain_id)]
chains$node <- Sys.info()[["nodename"]]
vars <- fit$metadata()$model_params
vars <- vars[grepl("^(b_|b\\[|Intercept|poutlier|sigmaone|tau|sd_)", vars) &
               !grepl("^(lprior|lp__|z_|r_)", vars)]
dd <- as.data.frame(fit$draws(vars, format = "draws_df"))
pars <- as.data.frame(posterior::summarise_draws(
  posterior::as_draws_df(dd[, vars, drop = FALSE]), median,
  ~stats::quantile(.x, probs = c(0.05, 0.95)), "rhat", "ess_bulk", "ess_tail"))
spread <- vapply(vars, function(v) {
  m <- tapply(dd[[v]], dd$.chain, stats::median)
  diff(range(m)) / stats::sd(dd[[v]])
}, numeric(1))
pars$chain_spread <- spread[pars$variable]
utils::write.csv(chains, paste0(stem, "__chains.csv"), row.names = FALSE)
if (name %in% c("rdm", "rdm_sub", "lba2")) {
  utils::write.csv(pp_race(name, dd), paste0(stem, "__pp.csv"), row.names = FALSE)
}
utils::write.csv(pars, paste0(stem, "__pars.csv"), row.names = FALSE)
log_msg(sprintf("done: div %d, max rhat %.3f, min ess_bulk %.0f, leapfrog %d",
                sum(chains$div), max(pars$rhat, na.rm = TRUE),
                min(pars$ess_bulk, na.rm = TRUE), sum(chains$leap_warm + chains$leap_samp)))
