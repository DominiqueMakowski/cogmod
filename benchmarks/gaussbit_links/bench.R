# Probit or logit for the choice half of cogmod_gaussbit()?
# ============================================================
#
# Question (2026-09-24): is the probit version actually cheaper to SAMPLE from
# than a "Gaussistic" one - Gaussian RT, logistic choice - especially in mixed
# models? "Gaussistic" is read two ways, and both are fitted:
#
#   gl_rho     the logit twin of the family: the same Gaussian copula, with
#              the choice marginal logistic instead of probit, so that
#              P(dec = 1 | t) = Phi((q + r u) / sqrt(1 - r^2)),
#              q = Phi^-1(logistic(mudec)). Defined in this script only.
#   mv_logit   what people fit: brms' native gaussian() + bernoulli("logit"),
#              independent.
#
# against
#
#   gp_rho     cogmod_gaussbit(), rho free
#   gp_0       cogmod_gaussbit(), rho = 0 (same model as mv_probit)
#   mv_probit  brms' native gaussian() + bernoulli("probit")
#
# Efficiency is split into its two factors, because they have different causes
# and different fixes:
#
#   leapfrogs per effective sample - the GEOMETRY; noise-free (a count)
#   microseconds per gradient      - the IMPLEMENTATION; timed in alternating
#                                    blocks at a posterior draw, as
#                                    gradient_cost.R does
#
# Their product is CPU time per effective sample. Wall-clock ESS/s from the
# fits themselves is reported too but is the noisiest of the three.
#
# Usage:
#   Rscript benchmarks/gaussbit_links/bench.R [--datasets 3] [--warmup 1000]
#     [--sampling 1000] [--out benchmarks/gaussbit_links/results]
#     [--settings fixed,mixed] [--models gp_rho,gl_rho,gp_0,mv_probit,mv_logit]

args <- commandArgs(trailingOnly = TRUE)
opt <- list(datasets = 3L, warmup = 1000L, sampling = 1000L,
            out = "benchmarks/gaussbit_links/results",
            settings = "fixed,mixed",
            models = "gp_rho,gl_rho,gp_0,mv_probit,mv_logit",
            reps = 21L, block = 20L)
for (i in 2 * seq_len(length(args) %/% 2) - 1) {
  key <- sub("^--", "", args[i])
  if (!key %in% names(opt)) stop("unknown option ", args[i])
  opt[[key]] <- methods::as(args[i + 1], class(opt[[key]]))
}
settings <- strsplit(opt$settings, ",")[[1]]
model_names <- strsplit(opt$models, ",")[[1]]
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

suppressPackageStartupMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(brms)
})
log_line <- function(...) {
  cat(format(Sys.time(), "%H:%M:%S"), ..., "\n")
  flush.console()
}


# The logit twin ----------------------------------------------------------

# q = Phi^-1(logistic(mudec)), formed in log space on whichever side keeps its
# digits: std_normal_log_qf() of log_inv_logit() is exact in the lower tail,
# and the upper tail is the same thing by symmetry. At rho = 0 the choice term
# is then log Phi(+/- q) = log_inv_logit(+/- mudec), the Bernoulli-logit.
gaussistic <- function() {
  brms::custom_family(
    "cogmod_gaussistic", dpars = c("mu", "sigma", "mudec", "rho"),
    links = c("identity", "log", "identity", "identity"),
    lb = c(NA, 0, NA, NA), ub = c(NA, NA, NA, NA), type = "real", vars = "dec[n]"
  )
}
gaussistic_stanvars <- function() {
  brms::stanvar(block = "functions", scode = paste0(cogmod:::.LOG_PHI_STAN_PRELUDE, "
real cogmod_gaussistic_lpdf(real Y, real mu, real sigma, real mudec, real rho,
                            int dec) {
    if (sigma <= 0) return negative_infinity();
    real u = (Y - mu) / sigma;
    real q = mudec < 0 ? std_normal_log_qf(log_inv_logit(mudec))
                       : -std_normal_log_qf(log_inv_logit(-mudec));
    real a = q * cosh(rho) + u * sinh(rho);
    return -0.5 * square(u) - log(sigma) - 0.91893853320467274
           + cogmod_log_Phi(dec == 1 ? a : -a);
}
"))
}
# R reference, for the check below
ldens_gaussistic <- function(t, k, mu, sigma, mudec, rho) {
  r <- tanh(rho)
  u <- (t - mu) / sigma
  z <- (stats::qnorm(stats::plogis(mudec)) + r * u) / sqrt(1 - r^2)
  # pnorm()'s lower.tail is not vectorised, so the side is chosen by sign
  stats::dnorm(t, mu, sigma, log = TRUE) + stats::pnorm((2 * k - 1) * z, log.p = TRUE)
}


# Data --------------------------------------------------------------------

# Simulated from cogmod_gaussbit(), so the probit models are the correctly
# specified ones; the logit models differ only in the shape of the link, which
# at these error rates (about 5-20%) is a rescaling by ~1.7 and not a misfit.
#
# `mixed`: 30 participants x 120 trials, two conditions. Random intercepts and
# condition slopes on the RT mean and on the choice, random intercepts on
# sigma. With ~10% errors that is about six errors per participant per cell,
# so the choice slopes' SD is weakly identified - the funnel a mixed model of
# accuracy usually has, and the place a difference in link could show.
# `fixed`: the same 3600 trials, no participant structure.
simulate <- function(setting, seed) {
  set.seed(seed)
  n_id <- 30; n_trial <- 120
  d <- expand.grid(trial = seq_len(n_trial), id = seq_len(n_id))
  d$Condition <- factor(ifelse(d$trial %% 2 == 0, "hard", "easy"))
  x <- as.numeric(d$Condition == "hard")
  if (setting == "mixed") {
    re <- function(sd) stats::rnorm(n_id, 0, sd)[d$id]
    mu <- 0.6 + re(0.08) + (0.05 + re(0.03)) * x
    sigma <- exp(log(0.15) + re(0.2))
    mudec <- -1.3 + re(0.4) + (0.4 + re(0.2)) * x
  } else {
    mu <- 0.6 + 0.05 * x
    sigma <- 0.15
    mudec <- -1.3 + 0.4 * x
  }
  sim <- rcogmod_gaussbit(nrow(d), mu = mu, sigma = sigma, mudec = mudec, rho = 0.3)
  d$RT <- sim$rt
  d$Error <- sim$response
  d$id <- factor(d$id)
  d
}


# Models ------------------------------------------------------------------

spec <- function(name, setting) {
  mixed <- setting == "mixed"
  re_mu <- if (mixed) " + (Condition | id)" else ""
  re_sigma <- if (mixed) " + (1 | id)" else ""
  f_rt <- as.formula(paste("RT ~ Condition", re_mu))
  f_sigma <- as.formula(paste("sigma ~ 1", re_sigma))
  f_dec <- as.formula(paste("mudec ~ Condition", re_mu))
  custom <- function(fam, rho) {
    lhs <- as.formula(paste("RT | dec(Error) ~ Condition", re_mu))
    if (identical(rho, "free")) {
      bf(lhs, f_sigma, f_dec, rho ~ 1, family = fam)
    } else {
      bf(lhs, f_sigma, f_dec, rho = 0, family = fam)
    }
  }
  native <- function(link) {
    bf(f_rt, f_sigma) +
      bf(as.formula(paste("Error ~ Condition", re_mu)), family = bernoulli(link)) +
      set_rescor(FALSE)
  }
  switch(name,
    gp_rho = list(f = custom(cogmod_gaussbit(), "free"),
                  sv = cogmod_gaussbit_stanvars(), prior = "cogmod"),
    gp_0 = list(f = custom(cogmod_gaussbit(), 0),
                sv = cogmod_gaussbit_stanvars(), prior = "cogmod"),
    # the family's own priors, restated: brms' defaults, the bernoulli()
    # intercept on mudec, normal(0, 0.5) on rho
    gl_rho = list(f = custom(gaussistic(), "free"), sv = gaussistic_stanvars(),
                  prior = c(prior(student_t(3, 0, 2.5), class = "Intercept", dpar = "mudec"),
                            prior(normal(0, 0.5), class = "Intercept", dpar = "rho"))),
    mv_probit = list(f = native("probit"), sv = NULL, prior = NULL),
    mv_logit = list(f = native("logit"), sv = NULL, prior = NULL)
  )
}

program <- function(name, setting, d) {
  s <- spec(name, setting)
  # The simulated RTs are Gaussian, so a participant drawn with a fast mean and
  # a wide sigma puts a handful of trials (2-3 in 3600) below zero, and
  # cogmod_priors()'s data check warns about them, as it should on real data.
  # Here they are the model's own tail, and the likelihood has density there.
  prior <- if (identical(s$prior, "cogmod")) {
    withCallingHandlers(
      cogmod_priors(s$f, d),
      warning = function(w) {
        if (grepl("zero or negative", conditionMessage(w))) invokeRestart("muffleWarning")
      }
    )
  } else {
    s$prior
  }
  if (is.null(prior)) prior <- brms::empty_prior()
  list(
    code = as.character(make_stancode(s$f, data = d, prior = prior, stanvars = s$sv,
                                      backend = "cmdstanr")),
    data = lapply(unclass(make_standata(s$f, data = d, prior = prior, stanvars = s$sv)),
                  identity)
  )
}


# Check the logit twin against its R reference ----------------------------

check_gaussistic <- function() {
  code <- paste0("functions {\n", gaussistic_stanvars()[[1]]$scode, "\n}")
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code), quiet = TRUE)
  mod$expose_functions()
  f <- mod$functions$cogmod_gaussistic_lpdf
  g <- expand.grid(t = c(0.2, 0.6, 1.5, 5), mudec = c(-6, -1.3, 0, 0.8, 4),
                   rho = c(-1, 0, 0.4, 2), k = 0:1)
  stan <- mapply(function(t, m, r, k) f(t, 0.6, 0.15, m, r, as.integer(k)),
                 g$t, g$mudec, g$rho, g$k)
  ref <- ldens_gaussistic(g$t, g$k, 0.6, 0.15, g$mudec, g$rho)
  err <- max(abs(stan - ref) / pmax(1, abs(ref)))
  # at rho = 0 it has to be Gaussian x Bernoulli-logit
  z <- g$rho == 0
  lg <- stats::dnorm(g$t[z], 0.6, 0.15, log = TRUE) +
    stats::plogis(ifelse(g$k[z] == 1, g$mudec[z], -g$mudec[z]), log.p = TRUE)
  err0 <- max(abs(stan[z] - lg))
  log_line(sprintf("gaussistic Stan vs R: max rel err %.1e; vs Bernoulli-logit at rho = 0: %.1e",
                   err, err0))
  stopifnot(err < 1e-10, err0 < 1e-10)
}


# Fitting -----------------------------------------------------------------

# Population-level quantities in the raw CmdStan names: `b` (the response's
# slopes; `b_<dpar>` for the others and for the mv model's responses), the
# `b_*_Intercept` generated quantities, group-level SDs and correlations. The
# diagonal of a correlation matrix is a constant and has no ESS.
pop_vars <- function(vars) {
  keep <- grepl("^(b$|b\\[|b_|sd_|Cor_)", vars)
  diag <- grepl("^Cor_\\d+\\[(\\d+),\\1\\]$", vars)
  vars[keep & !diag]
}

# Which of them belong to the choice half. brms numbers the group-level terms
# in formula order - mu, sigma, then the choice - in both the custom families
# and the mv models, so the choice's SDs and correlations are sd_3 and Cor_3.
is_choice <- function(vars) grepl("mudec|Error|rho|^sd_3|^Cor_3", vars)

fit_one <- function(mod, prog, setting, name, ds, seed) {
  fit <- mod$sample(data = prog$data, chains = 4, parallel_chains = 4,
                    iter_warmup = opt$warmup, iter_sampling = opt$sampling,
                    seed = seed, refresh = 0, show_messages = FALSE,
                    show_exceptions = FALSE, save_warmup = TRUE)
  sd <- fit$sampler_diagnostics(inc_warmup = TRUE, format = "draws_df")
  warm <- sd$.iteration <= opt$warmup
  vars <- pop_vars(fit$metadata()$stan_variables)
  vars <- pop_vars(posterior::variables(fit$draws(vars)))
  s <- posterior::summarise_draws(fit$draws(vars), "ess_bulk", "ess_tail", "rhat")
  tm <- fit$time()$chains
  lf_samp <- sum(sd$n_leapfrog__[!warm])
  choice <- is_choice(s$variable)
  rho_est <- if ("b_rho_Intercept" %in% s$variable) {
    mean(fit$draws("b_rho_Intercept"))
  } else NA_real_
  row <- data.frame(
    setting = setting, model = name, dataset = ds,
    divergent = sum(sd$divergent__[!warm]),
    max_treedepth_hits = sum(sd$treedepth__[!warm] >= 10),
    stepsize = mean(sd$stepsize__[!warm]),
    leapfrog_per_iter = mean(sd$n_leapfrog__[!warm]),
    leapfrog_warmup = sum(sd$n_leapfrog__[warm]),
    min_ess_bulk = min(s$ess_bulk), min_ess_tail = min(s$ess_tail),
    median_ess_bulk = stats::median(s$ess_bulk),
    min_ess_bulk_choice = min(s$ess_bulk[choice]),
    min_ess_bulk_rt = min(s$ess_bulk[!choice]),
    worst_var = s$variable[which.min(s$ess_bulk)],
    max_rhat = max(s$rhat),
    sampling_s = sum(tm$sampling), warmup_s = sum(tm$warmup),
    rho_z = rho_est
  )
  row$leapfrog_per_ess <- lf_samp / row$min_ess_bulk
  row$ess_per_s <- row$min_ess_bulk / row$sampling_s
  list(row = row, fit = fit)
}


# Gradient cost at a posterior draw, alternating blocks -------------------

time_gradients <- function(fits) {
  up <- lapply(fits, function(f) {
    f$init_model_methods(verbose = FALSE)
    um <- f$unconstrain_draws(format = "draws_matrix")
    as.numeric(um[nrow(um), ])
  })
  tb <- function(f, u, n) {
    t0 <- Sys.time()
    for (i in seq_len(n)) f$grad_log_prob(u, jacobian = TRUE)
    as.numeric(Sys.time() - t0, units = "secs")
  }
  for (k in names(fits)) tb(fits[[k]], up[[k]], 5)
  secs <- setNames(vector("list", length(fits)), names(fits))
  for (r in seq_len(opt$reps)) {
    for (k in if (r %% 2) names(fits) else rev(names(fits))) {
      secs[[k]] <- c(secs[[k]], tb(fits[[k]], up[[k]], opt$block))
    }
  }
  vapply(secs, function(s) 1e6 * stats::median(s) / opt$block, numeric(1))
}


# Main --------------------------------------------------------------------

if ("gl_rho" %in% model_names) check_gaussistic()

rows <- list()
grad <- list()
for (setting in settings) {
  d1 <- simulate(setting, seed = 1)
  mods <- list()
  for (nm in model_names) {
    p <- program(nm, setting, d1)
    stan <- file.path(opt$out, sprintf("%s_%s.stan", setting, nm))
    writeLines(p$code, stan)
    log_line("compiling", setting, nm)
    mods[[nm]] <- cmdstanr::cmdstan_model(stan, compile_model_methods = TRUE,
                                          dir = opt$out, quiet = TRUE)
  }
  first_fits <- list()
  for (ds in seq_len(opt$datasets)) {
    d <- if (ds == 1) d1 else simulate(setting, seed = ds)
    for (nm in model_names) {
      p <- program(nm, setting, d)
      res <- fit_one(mods[[nm]], p, setting, nm, ds, seed = 100 + ds)
      r <- res$row
      rows[[length(rows) + 1]] <- r
      if (ds == 1) first_fits[[nm]] <- res$fit
      log_line(sprintf("%-6s %-9s ds%d  div %3d  lf/iter %5.1f  minESS %5.0f  lf/ESS %6.1f  samp %6.1fs  rho %s",
                       setting, nm, ds, r$divergent, r$leapfrog_per_iter, r$min_ess_bulk,
                       r$leapfrog_per_ess, r$sampling_s,
                       if (is.na(r$rho_z)) "-" else sprintf("%.2f", r$rho_z)))
      utils::write.csv(do.call(rbind, rows), file.path(opt$out, "fits.csv"), row.names = FALSE)
    }
  }
  log_line("timing gradients,", setting)
  g <- time_gradients(first_fits)
  grad[[setting]] <- data.frame(setting = setting, model = names(g), us_per_gradient = g)
  utils::write.csv(do.call(rbind, grad), file.path(opt$out, "gradients.csv"), row.names = FALSE)
  print(grad[[setting]])
}

res <- do.call(rbind, rows)
gr <- do.call(rbind, grad)
agg <- stats::aggregate(
  cbind(divergent, leapfrog_per_iter, min_ess_bulk, min_ess_tail, leapfrog_per_ess,
        ess_per_s, max_rhat) ~ setting + model,
  data = res, FUN = stats::median
)
agg <- merge(agg, gr, by = c("setting", "model"))
agg$cpu_ms_per_ess <- agg$leapfrog_per_ess * agg$us_per_gradient / 1000
utils::write.csv(agg, file.path(opt$out, "summary.csv"), row.names = FALSE)
log_line("done")
print(agg, digits = 3)
