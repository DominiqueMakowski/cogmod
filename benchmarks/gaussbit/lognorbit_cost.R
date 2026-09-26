# Would a descriptive "lognorbit" baseline - shifted lognormal RT + probit
# choice, correlated tanh(rho) on log(RT - ndt), with the choice families'
# outlier component - be cheaper to fit than the LNR? Same data (the
# decision_making vignette's), the vignette's LNR formula, the same ndt and
# poutlier priors; three fits each, alternating order, then grad_log_prob() in
# 21 alternating blocks. See README.md. About 15 minutes.
#
#   Rscript benchmarks/gaussbit/lognorbit_cost.R          # or `show` to print the program
#
# The lognorbit density is a prototype written here, never a family: it had 3
# warmup rejections from a NaN in log_mix() (an unguarded u * sinh(rho)
# overflow), which a real version would have to handle.
suppressMessages({
  pkgload::load_all(".", quiet = TRUE)
  library(brms); library(cmdstanr)
})
out <- file.path("benchmarks", "gaussbit", "results")
dir.create(out, showWarnings = FALSE)
args <- commandArgs(trailingOnly = TRUE)
stage <- if (length(args)) args[1] else "all"

data(speed_acc, package = "rtdists")
df <- data.frame(
  Participant = as.integer(as.character(speed_acc$id)),
  Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[as.character(speed_acc$condition)]),
  RT = speed_acc$rt,
  Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat))
)
df <- df[df$Participant %in% c(1, 2, 3) & df$RT <= 2, ]

# ---- LNR, the vignette's formula ---------------------------------------------
f_lnr <- bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1,
            sigmaone ~ 1, sigmabias = 0, ndt ~ Condition, family = cogmod_lnr())
pr_lnr <- suppressMessages(cogmod_priors(f_lnr, df))
sv_lnr <- cogmod_stanvars(f_lnr)
init_lnr <- cogmod_inits(f_lnr, df)

# ---- Lognorbit prototype -----------------------------------------------------
# Shifted lognormal RT + probit choice (Gaussian copula on log(RT - ndt)),
# mixed with the choice families' outlier component, written the way
# .choice_lpdf() writes the LNR's.
lc <- log(2) - 0.5 * log(2 * pi * .POUTLIER_SCALE^2) - log(2)
lnb_code <- paste0(.LOG_PHI_STAN_PRELUDE, sprintf("
real cogmod_lognorbit_lpdf(real Y, real mu, real sigma, real mudec, real rho,
                           real ndt, real poutlier, int dec) {
    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) return negative_infinity();
    if (dec < 0 || dec > 1) return negative_infinity();
    if (Y <= 0) return negative_infinity();
    real lp_out = %s - %s * square(Y);
    real t_adj = Y - ndt;
    if (t_adj <= 0) return log(poutlier) + lp_out;
    real lt = log(t_adj);
    real u = (lt - mu) / sigma;
    real a = mudec * cosh(rho) + u * sinh(rho);
    real lp_dec = -0.5 * square(u) - log(sigma) - 0.91893853320467274 - lt
                  + cogmod_log_Phi(dec == 1 ? a : -a);
    return log_mix(poutlier, lp_out, lp_dec);
}
", formatC(lc, format = "g", digits = 17), formatC(1 / (2 * .POUTLIER_SCALE^2), format = "g", digits = 15)))
lognorbit <- custom_family(
  "cogmod_lognorbit", dpars = c("mu", "sigma", "mudec", "rho", "ndt", "poutlier"),
  links = c("identity", "log", "identity", "identity", "log", "logit"),
  lb = c(NA, 0, NA, NA, 0, 0), ub = c(NA, NA, NA, NA, NA, 1),
  type = "real", vars = "dec[n]"
)
f_lnb <- bf(RT | dec(Error) ~ Condition, sigma ~ Condition, mudec ~ Condition,
            rho ~ Condition, ndt ~ Condition, family = lognorbit)
sv_lnb <- stanvar(scode = lnb_code, block = "functions")
# ndt and poutlier priors copied from the LNR's; the rest weakly informative
# on the scales the lognormal half lives on (log seconds of decision time).
keep <- pr_lnr[(pr_lnr$dpar == "ndt" | pr_lnr$class == "poutlier") & nzchar(pr_lnr$prior), ]
pr_lnb <- c(
  prior(normal(-1, 1), class = Intercept),
  prior(normal(0, 1), class = b),
  prior(normal(-0.7, 0.7), class = Intercept, dpar = sigma),
  prior(normal(0, 0.5), class = b, dpar = sigma),
  prior(student_t(3, 0, 2.5), class = Intercept, dpar = mudec),
  prior(normal(0, 1), class = b, dpar = mudec),
  prior(normal(0, 0.5), class = Intercept, dpar = rho),
  prior(normal(0, 0.5), class = b, dpar = rho),
  keep
)
init_lnb <- function(chain_id = 1) {
  j <- function() stats::rnorm(1, 0, 0.1)
  list(Intercept = -1 + j(), b = array(0, 1), Intercept_sigma = log(0.4) + j(),
       b_sigma = array(0, 1), Intercept_mudec = qnorm(mean(df$Error)) + j(),
       b_mudec = array(0, 1), Intercept_rho = j(), b_rho = array(0, 1),
       Intercept_ndt = log(0.25) + j(), b_ndt = array(0, 1), poutlier = 0.01)
}

models <- list(
  LNR = list(f = f_lnr, pr = pr_lnr, sv = sv_lnr, init = init_lnr),
  Lognorbit = list(f = f_lnb, pr = pr_lnb, sv = sv_lnb, init = init_lnb)
)
if (stage == "show") {
  print(keep)
  cat(make_stancode(f_lnb, data = df, prior = pr_lnb, stanvars = sv_lnb))
  quit(save = "no")
}

mods <- list(); sdat <- list()
for (k in names(models)) {
  m <- models[[k]]
  code <- make_stancode(m$f, data = df, prior = m$pr, stanvars = m$sv)
  sdat[[k]] <- unclass(make_standata(m$f, data = df, prior = m$pr, stanvars = m$sv))
  mods[[k]] <- cmdstan_model(write_stan_file(code, dir = out, basename = k),
                             compile_model_methods = TRUE, quiet = TRUE,
                             # model methods cannot be stood up on an
                             # executable left over from an earlier run
                             force_recompile = TRUE)
}

# ---- Fits: three seeds each, alternating order ---------------------------------
runs <- list(c("LNR", 1), c("Lognorbit", 1), c("Lognorbit", 2), c("LNR", 2), c("LNR", 3), c("Lognorbit", 3))
rows <- list(); fits <- list()
for (r in runs) {
  k <- r[1]; s <- as.integer(r[2])
  set.seed(s)
  fit <- mods[[k]]$sample(data = sdat[[k]], init = models[[k]]$init, seed = s,
                          chains = 4, parallel_chains = 4, iter_warmup = 500,
                          iter_sampling = 500, refresh = 0, show_messages = FALSE)
  if (is.null(fits[[k]])) fits[[k]] <- fit
  sm <- fit$summary(NULL, c("ess_bulk", "ess_tail", "rhat"))
  sm <- sm[!grepl("^lp__|^lprior|^Intercept$|^Intercept_", sm$variable) | grepl("^b_|^Intercept", sm$variable), ]
  sm <- sm[!sm$variable %in% c("lp__", "lprior"), ]
  sm <- sm[is.finite(sm$ess_bulk) & is.finite(sm$ess_tail), ]
  ess <- pmin(sm$ess_bulk, sm$ess_tail)
  sd <- fit$sampler_diagnostics(format = "df")
  tm <- fit$time()$chains
  rows[[length(rows) + 1]] <- data.frame(
    model = k, seed = s,
    wall_s = round(fit$time()$total, 1),
    cpu_s = round(sum(tm$total), 1),
    warmup_s = round(sum(tm$warmup), 1), sampling_s = round(sum(tm$sampling), 1),
    leapfrog_per_iter = round(mean(sd$n_leapfrog__), 1),
    min_ess = round(min(ess)), worst = sm$variable[which.min(ess)],
    max_rhat = round(max(sm$rhat, na.rm = TRUE), 3),
    divergent = sum(sd$divergent__),
    cpu_ms_per_ess = round(1000 * sum(tm$sampling) / min(ess), 1)
  )
  print(rows[[length(rows)]], row.names = FALSE)
}
res <- do.call(rbind, rows)
write.csv(res, file.path(out, "lognorbit_fits.csv"), row.names = FALSE)

# ---- Gradient cost, alternating blocks ---------------------------------------
up <- list()
for (k in names(fits)) {
  fits[[k]]$init_model_methods(verbose = FALSE)
  up[[k]] <- as.numeric(fits[[k]]$unconstrain_draws(format = "draws_matrix")[1, ])
}
tb <- function(k, n) { t0 <- Sys.time(); for (i in seq_len(n)) fits[[k]]$grad_log_prob(up[[k]], jacobian = TRUE); as.numeric(Sys.time() - t0, units = "secs") }
for (k in names(fits)) tb(k, 20)
secs <- list(LNR = numeric(0), Lognorbit = numeric(0))
for (r in 1:21) {
  ord <- if (r %% 2) c("LNR", "Lognorbit") else c("Lognorbit", "LNR")
  for (k in ord) secs[[k]] <- c(secs[[k]], tb(k, 100))
}
us <- vapply(secs, function(s) 1e6 * median(s) / 100, numeric(1))
cat("\nus per gradient (N =", nrow(df), "):", round(us), " ratio LNR/Lognorbit:", round(us[1] / us[2], 2), "\n")
cat("per observation (us):", round(us / nrow(df), 3), "\n")
write.csv(data.frame(model = names(us), us_per_gradient = round(us), N = nrow(df)),
          file.path(out, "lognorbit_gradients.csv"), row.names = FALSE)
print(res, row.names = FALSE)

# Sanity: posterior means of the lognorbit ndt against the ML fit (0.341, 0.299)
dr <- fits$Lognorbit$draws(c("Intercept_ndt", "b_ndt[1]", "Intercept_rho", "b_rho[1]", "poutlier"), format = "df")
cat("\nLognorbit posterior means: ndt (Acc, Speed; brms Intercept is centred) ",
    round(colMeans(as.data.frame(dr)[, 1:5]), 3), "\n")
