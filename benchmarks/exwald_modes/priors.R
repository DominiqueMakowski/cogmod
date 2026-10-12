# Which prior moves the mass back onto WW? Each of 0.3.5's four modes
# re-optimised under a variant of its priors, then Laplace log-mass. The
# variants change only the drift `mu` (softplus link) - its intercept, whose
# prior brms otherwise sets to student_t(3, 0, 2.5), and its slope, which is
# otherwise flat - or the centre of `tau`'s intercept prior.
#
#   Rscript benchmarks/exwald_modes/priors.R       (needs modes.R's 0.3.5 .rds)
#
# Writes results/prior_variants.csv. On 0.3.4's modes Laplace agrees with
# bridge sampling to 0.03 in log-mass for WW and EW but only 0.4 for EE, the
# widest ridge (README.md), so treat a share as a ranking rather than a
# measurement.
source("benchmarks/exwald_modes/common.R")
md <- readRDS(file.path(OUT, "modes_035.rds"))
V <- list(
  "0.3.5"                         = list(),
  "tau centre -2.5 (0.079 s)"     = list(tau = -2.5),
  "tau centre -3 (0.049 s)"       = list(tau = -3),
  "mu slope normal(0, 2)"         = list(slope = 2),
  "mu slope normal(0, 1)"         = list(slope = 1),
  "mu intercept normal(3, 1)"     = list(icpt = 1),
  "intercept (3, 1.5) + slope (0, 2)" = list(icpt = 1.5, slope = 2),
  "intercept (3, 1) + slope (0, 2)"   = list(icpt = 1, slope = 2),
  "intercept (3, 1) + slope (0, 1) = P1" = list(icpt = 1, slope = 1)
)
lp_v <- function(u, v) {
  lp <- lp_fun(u, "0.3.5")
  if (!is.null(v$icpt)) {
    lp <- lp - (stats::dt(u[2] / 2.5, 3, log = TRUE) - log(2.5)) + stats::dnorm(u[2], 3, v$icpt, log = TRUE)
  }
  if (!is.null(v$slope)) lp <- lp + stats::dnorm(u[1], 0, v$slope, log = TRUE)
  if (!is.null(v$tau)) {
    lp <- lp - stats::dnorm(u[6], -1.5, 0.7, log = TRUE) + stats::dnorm(u[6], v$tau, 0.7, log = TRUE)
  }
  lp
}
out <- NULL
for (vn in names(V)) {
  nlp <- function(u) { r <- -lp_v(u, V[[vn]]); if (!is.finite(r)) 1e10 else r }
  for (k in seq_along(md$par)) {
    o <- stats::optim(md$par[[k]], nlp, method = "BFGS",
                      control = list(maxit = 3000, reltol = 1e-15, ndeps = rep(1e-5, 9)))
    ev <- eigen(stats::optimHess(o$par, nlp, control = list(ndeps = rep(1e-4, 9))),
                symmetric = TRUE, only.values = TRUE)$values
    n <- natpar(o$par)
    out <- rbind(out, data.frame(variant = vn, mode = md$label[k], still = mode_label(o$par) == md$label[k],
                                 lp = -o$value, loglik = loglik_fun(o$par), min_eig = min(ev),
                                 logmass = -o$value - 0.5 * sum(log(pmax(ev, 1e-12))) + 4.5 * log(2 * pi),
                                 muA = n$mu[1], muS = n$mu[2], aA = n$boundary[1], aS = n$boundary[2],
                                 tauA = n$tau[1], tauS = n$tau[2], ndtA = n$ndt[1], ndtS = n$ndt[2]))
  }
}
out$share <- stats::ave(out$logmass, out$variant, FUN = function(l) exp(l - max(l)) / sum(exp(l - max(l))))
num <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, 6)
print(out[, c("variant", "mode", "still", "lp", "loglik", "min_eig", "share")], row.names = FALSE)
utils::write.csv(out, file.path(OUT, "prior_variants.csv"), row.names = FALSE)
