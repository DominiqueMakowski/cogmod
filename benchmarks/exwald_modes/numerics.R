# Is the valley between the modes a property of the ex-Wald, or of how it is
# computed? The valley sits on the seam mu^2 = 2 / tau where the density
# switches from the Faddeeva branch to the closed form, so both branches and
# the switch are checked:
#
#  1. .dexwald_raw() against brute-force numerical convolution (adaptive
#     quadrature split at the Wald mode, rel.tol 1e-12), trial by trial, at
#     every cell of every mode and along the straight line WW -> EW.
#  2. The log-likelihood as the Accuracy drift steps through the seam.
#  3. The Stan program (model methods): log_prob against the R twin, and its
#     gradient against central differences, on both sides of the seam.
#
#   Rscript benchmarks/exwald_modes/numerics.R      (needs modes.R's 0.3.5 .rds)
#
# Writes results/numerics_*.csv.
source("benchmarks/exwald_modes/common.R")
md <- readRDS(file.path(OUT, "modes_035.rds"))
modes <- stats::setNames(md$par, md$label)

conv_ldens <- function(t, mu, a, tau) {
  g <- 1 / tau
  f <- function(u) exp(log(a) - 0.5 * (log(2 * pi) + 3 * log(u)) - (a - mu * u)^2 / (2 * u)) *
    g * exp(-g * (t - u))
  m <- a / mu * (sqrt(1 + (3 / (2 * a * mu))^2) - 3 / (2 * a * mu))
  br <- sort(unique(c(0, pmin(c(m / 4, m, 2 * m, t / 2), t), t)))
  log(sum(vapply(seq_len(length(br) - 1), function(j) {
    stats::integrate(f, br[j], br[j + 1], rel.tol = 1e-12, abs.tol = 0, subdivisions = 2000L)$value
  }, numeric(1))))
}
brute <- function(where, u, n_trials = 150) {
  n <- natpar(u)
  do.call(rbind, lapply(1:2, function(k) {
    t <- Y[D$Condition == n$cond[k]] - n$ndt[k]
    t <- t[t > 0]
    set.seed(1)
    t <- t[sample(length(t), min(n_trials, length(t)))]
    d1 <- .dexwald_raw(t, n$mu[k], n$boundary[k], n$tau[k])
    d2 <- vapply(t, conv_ldens, numeric(1), mu = n$mu[k], a = n$boundary[k], tau = n$tau[k])
    data.frame(where = where, cell = n$cond[k], k2 = n$mu[k]^2 - 2 / n$tau[k],
               max_abs_dlogdens = max(abs(d1 - d2)))
  }))
}
b1 <- do.call(rbind, lapply(names(modes), function(m) brute(paste("mode", m), modes[[m]])))
b2 <- do.call(rbind, lapply(seq(0.1, 0.9, by = 0.1), function(s) {
  brute(sprintf("line WW->EW s=%.1f", s), (1 - s) * modes$WW + s * modes$EW, 60)[1, ]
}))
bf <- rbind(b1, b2)
print(bf, row.names = FALSE)
utils::write.csv(bf, file.path(OUT, "numerics_bruteforce.csv"), row.names = FALSE)

# Through the seam: the Accuracy cell at tau 0.17, boundary 0.55, ndt 0.345
# (the neighbourhood of the valley), Speed at its WW values; the Accuracy
# drift stepped through sqrt(2 / tau) by a relative offset.
nS <- natpar(modes$WW)[2, ]
mk_u <- function(muA) {
  pair <- function(eA, eS) { b <- eS - eA; c(b, eA + MX * b) }
  c(pair(softplus_inv(muA), softplus_inv(nS$mu)), pair(softplus_inv(0.55), softplus_inv(nS$boundary)),
    pair(softplus_inv(0.17), softplus_inv(nS$tau)), pair(log(0.345), log(nS$ndt)),
    stats::qlogis(nS$poutlier))
}
prog <- make_program("0.3.5")
mmdir <- file.path(tempdir(), "exw_mm")
dir.create(mmdir, showWarnings = FALSE)
mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(prog$code, dir = mmdir),
                               compile_model_methods = TRUE, force_recompile = TRUE, quiet = TRUE)
fit <- mod$sample(data = prog$sdata, init = list(as_init(modes$WW)), chains = 1, iter_warmup = 1,
                  iter_sampling = 1, fixed_param = TRUE, refresh = 0, show_messages = FALSE)
invisible(fit$metadata())
fit$init_model_methods(verbose = FALSE)
mu0 <- sqrt(2 / 0.17)
seam <- do.call(rbind, lapply(c(-1e-2, -1e-4, -1e-6, -1e-8, -1e-10, 1e-10, 1e-8, 1e-6, 1e-4, 1e-2), function(eps) {
  u <- mk_u(mu0 * (1 + eps))
  g <- fit$grad_log_prob(u, jacobian = TRUE)
  fd <- vapply(1:9, function(j) {
    e <- replace(numeric(9), j, 1e-6)
    (fit$log_prob(u + e, jacobian = TRUE) - fit$log_prob(u - e, jacobian = TRUE)) / 2e-6
  }, numeric(1))
  data.frame(rel_offset = eps, k2A = (mu0 * (1 + eps))^2 - 2 / 0.17,
             loglik_R = loglik_fun(u), lp_stan_minus_R = fit$log_prob(u, jacobian = TRUE) - lp_fun(u, "0.3.5"),
             dlp_dIntercept = g[2], max_rel_grad_err = max(abs(g - fd) / pmax(1, abs(fd))),
             grad_finite = all(is.finite(g)))
}))
print(format(seam, digits = 12), row.names = FALSE)
utils::write.csv(seam, file.path(OUT, "numerics_seam.csv"), row.names = FALSE)

at_modes <- do.call(rbind, lapply(names(modes), function(m) {
  u <- modes[[m]]
  data.frame(mode = m, lp_R = lp_fun(u, "0.3.5"), lp_stan = fit$log_prob(u, jacobian = TRUE),
             grad_norm = sqrt(sum(fit$grad_log_prob(u, jacobian = TRUE)^2)))
}))
print(format(at_modes, digits = 10), row.names = FALSE)
utils::write.csv(at_modes, file.path(OUT, "numerics_modes.csv"), row.names = FALSE)
