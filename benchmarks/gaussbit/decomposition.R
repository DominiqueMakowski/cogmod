# Where does the elpd gap between cogmod_gaussbit() and the LNR come from?
# Maximum likelihood on the decision_making vignette's data (speed_acc,
# participants 1-3, RT <= 2 s), every parameter ~ Condition, stepping from the
# Gaussian + probit model to the LNR one assumption at a time. See README.md.
#
#   Rscript benchmarks/gaussbit/decomposition.R
#
# Needs only the package as it stands: the gaussbit density is written out
# below, since the family itself is no longer in the tree (gaussbit.patch).

suppressMessages(pkgload::load_all(".", quiet = TRUE))
out <- file.path("benchmarks", "gaussbit", "results")
dir.create(out, showWarnings = FALSE)

data(speed_acc, package = "rtdists")
df <- data.frame(
  Participant = as.integer(as.character(speed_acc$id)),
  Condition = factor(as.character(speed_acc$condition)),
  RT = speed_acc$rt,
  Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat))
)
df <- df[df$Participant %in% c(1, 2, 3) & df$RT <= 2, ]
df$Condition <- droplevels(df$Condition)
X <- model.matrix(~ Condition, df)
spd <- df$Condition == "speed"
minrt <- tapply(df$RT, df$Condition, min)

# The gaussbit log-density, as it was in R/model_gaussbit.R: Gaussian t, probit
# choice, correlation tanh(rho) between the two.
gb_ldens <- function(t, k, mu, sigma, mudec, rho) {
  u <- (t - mu) / sigma
  a <- mudec * cosh(rho) + u * sinh(rho)
  ld <- -0.5 * u^2 - log(sigma) - 0.5 * log(2 * pi) + pnorm((2 * k - 1) * a, log.p = TRUE)
  ld[!is.finite(ld)] <- -Inf
  ld
}
gb_rng <- function(n, mu, sigma, mudec, rho) {
  u <- rnorm(n)
  w <- mudec + tanh(rho) * u + rnorm(n) / cosh(rho)
  data.frame(rt = mu + sigma * u, response = as.numeric(w > 0))
}
log_mix2 <- function(a, b) pmax(a, b) + log1p(exp(-abs(a - b)))
# The choice families' outlier component: half-Normal(0.2 s) RT, 1/2 per choice.
lp_out <- log(2) + dnorm(df$RT, 0, .POUTLIER_SCALE, log = TRUE) + log(0.5)

fit <- function(nll, st) {
  f <- function(p) { v <- nll(p); if (is.finite(v)) v else 1e10 }
  o <- optim(st, f, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  o <- optim(o$par, f, method = "Nelder-Mead", control = list(maxit = 20000, reltol = 1e-12))
  o <- optim(o$par, f, method = "BFGS", control = list(maxit = 10000, reltol = 1e-12))
  list(ll = -o$value, par = o$par, k = length(st))
}

# Gaussian or log-RT + probit; `shift` adds ndt ~ Condition, `outl` the
# outlier component. Parameters: mu (2), log sigma (2), mudec (2), rho (2),
# log ndt (2), logit poutlier (1).
lp_gb <- function(p, log_rt, rho = TRUE, shift = FALSE, outl = FALSE) {
  ndt <- if (shift) exp(X %*% p[9:10]) else 0
  y <- df$RT - ndt
  ld <- rep(-Inf, nrow(df)); ok <- y > 0
  r <- if (rho) (X %*% p[7:8])[ok] else 0
  yy <- if (log_rt) log(y[ok]) else y[ok]
  ld[ok] <- gb_ldens(yy, df$Error[ok], (X %*% p[1:2])[ok], exp(X %*% p[3:4])[ok],
                     (X %*% p[5:6])[ok], r) - (if (log_rt) log(y[ok]) else 0)
  if (!outl) return(ld)
  po <- plogis(p[11])
  log_mix2(log1p(-po) + ld, log(po) + lp_out)
}
st <- function(y) c(mean(y), 0, log(sd(y)), 0, qnorm(mean(df$Error)), 0, 0, 0)

g_free <- fit(function(p) -sum(lp_gb(p, FALSE)), st(df$RT))
g_0 <- fit(function(p) -sum(lp_gb(c(p, 0, 0), FALSE, rho = FALSE)), st(df$RT)[1:6])
l_free <- fit(function(p) -sum(lp_gb(p, TRUE)), st(log(df$RT)))
l_0 <- fit(function(p) -sum(lp_gb(c(p, 0, 0), TRUE, rho = FALSE)), st(log(df$RT))[1:6])
l_sh <- fit(function(p) -sum(lp_gb(p, TRUE, shift = TRUE, outl = TRUE)),
            c(l_free$par, log(0.15), 0, -5))

# The native brms equivalent of the free-rho Gaussian fit: heteroscedastic
# Gaussian on the RT, probit on the choice with the RT (x Condition) entered
# as a predictor.
o_rt <- optim(c(mean(df$RT), 0, log(sd(df$RT)), 0),
              function(p) -sum(dnorm(df$RT, X %*% p[1:2], exp(X %*% p[3:4]), log = TRUE)),
              method = "BFGS", control = list(reltol = 1e-12))
native_ll <- -o_rt$value +
  as.numeric(logLik(glm(Error ~ Condition * RT, family = binomial("probit"), data = df)))

# The vignette's LNR: nuzero ~ C, nuone ~ C, sigmazero, sigmaone, ndt ~ C,
# poutlier. Three starts; the best is kept.
lp_lnr <- function(p) dcogmod_lnr(df$RT, nuzero = X %*% p[1:2], nuone = X %*% p[3:4],
  sigmazero = exp(p[5]), sigmaone = exp(p[6]), ndt = exp(X %*% p[7:8]),
  response = df$Error, poutlier = plogis(p[9]), log = TRUE)
lnr <- NULL
for (s in list(c(1, 0, -0.5, 0, log(0.3), log(0.5), log(0.15), 0, -5),
               c(2, 0.5, 0, 0, log(0.3), log(0.4), log(0.2), 0, -6),
               c(0.5, 0.3, -1, 0.3, log(0.5), log(0.6), log(0.1), 0, -4))) {
  r <- suppressWarnings(fit(function(p) -sum(lp_lnr(p)), s))
  if (is.null(lnr) || r$ll > lnr$ll) lnr <- r
}

res <- data.frame(
  model = c("gaussbit, rho = 0 (= gaussian() + bernoulli('probit'))",
            "gaussbit, rho free",
            "native: gaussian() + probit(Error ~ Condition * RT)",
            "log-RT + probit, rho = 0",
            "log-RT + probit, rho free",
            "log-RT + probit + ndt ~ Condition + poutlier ('lognorbit')",
            "LNR (vignette formula)"),
  logLik = c(g_0$ll, g_free$ll, native_ll, l_0$ll, l_free$ll, l_sh$ll, lnr$ll),
  k = c(6, 8, 8, 6, 8, 11, 9)
)
res$AIC <- -2 * res$logLik + 2 * res$k
res$logLik <- round(res$logLik, 2); res$AIC <- round(res$AIC, 1)
print(res, row.names = FALSE)
write.csv(res, file.path(out, "decomposition.csv"), row.names = FALSE)

cat("\nrho (Fisher z), accuracy / speed: Gaussian", round(c(g_free$par[7], sum(g_free$par[7:8])), 3),
    "| lognorbit", round(c(l_sh$par[7], sum(l_sh$par[7:8])), 3), "\n")
cat("ndt, accuracy / speed: lognorbit", round(exp(c(l_sh$par[9], sum(l_sh$par[9:10]))), 3),
    "| LNR", round(exp(c(lnr$par[7], sum(lnr$par[7:8]))), 3), "\n")

cat("\nLNR minus lognorbit, summed log-likelihood by condition x response:\n")
print(round(tapply(lp_lnr(lnr$par) - lp_gb(l_sh$par, TRUE, shift = TRUE, outl = TRUE),
                   list(df$Condition, Error = df$Error), sum), 1))

# Error-minus-correct mean RT, observed and implied (simulated at the ML fit).
set.seed(1); n <- 2e5
cat("\nmean RT correct, error, error - correct (s):\n")
for (cnd in 0:1) {
  lab <- levels(df$Condition)[cnd + 1]; o <- df[df$Condition == lab, ]
  i <- c(1, cnd)
  g <- function(par, k) sum(par[k] * i)
  sl <- gb_rng(n, g(l_sh$par, 1:2), exp(g(l_sh$par, 3:4)), g(l_sh$par, 5:6), g(l_sh$par, 7:8))
  sl$rt <- exp(sl$rt) + exp(g(l_sh$par, 9:10))
  sr <- rcogmod_lnr(n, g(lnr$par, 1:2), g(lnr$par, 3:4), exp(lnr$par[5]), exp(lnr$par[6]),
                    exp(g(lnr$par, 7:8)))
  m <- function(rt, r) { a <- tapply(rt, r, mean); round(c(a, a[2] - a[1]), 3) }
  cat(sprintf("%-8s observed %s | lognorbit %s | LNR %s\n", lab,
              paste(m(o$RT, o$Error), collapse = " "),
              paste(m(sl$rt, sl$response), collapse = " "),
              paste(m(sr$rt, sr$response), collapse = " ")))
}
