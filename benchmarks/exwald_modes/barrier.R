# How high the wall between two modes is. lp__ profiled along one cell's tau
# (link scale): at each value the other 8 coordinates are optimised, sweeping
# from each mode and keeping the better. Every path between the modes has to
# cross each level set of tau, so the profile's minimum bounds the saddle from
# above and `barrier_*` from below.
#
#   Rscript benchmarks/exwald_modes/barrier.R <prior>      (needs modes.R's .rds)
#
# Writes results/barrier_<prior>.csv: WW -> EW along the Accuracy tau (the
# split the driftone_priors refits showed) and WW -> WE along the Speed tau.
source("benchmarks/exwald_modes/common.R")
prior <- commandArgs(trailingOnly = TRUE)[1]
tag <- gsub("[^0-9A-Za-z]", "", prior)
md <- readRDS(file.path(OUT, sprintf("modes_%s.rds", tag)))
modes <- stats::setNames(md$par, md$label)
# tau's two coordinates re-expressed as the cells' link values, one profiled
to_u <- function(v8, eta, cell) {
  etaA <- if (cell == 1) eta else v8[5]
  etaS <- if (cell == 1) v8[5] else eta
  b <- etaS - etaA
  c(v8[1:4], b, etaA + MX * b, v8[6:8])
}
from_u <- function(u, cell) {
  eta <- c(u[6] - MX * u[5], u[6] + (1 - MX) * u[5])
  c(u[1:4], eta[3 - cell], u[7:9])
}
profile <- function(from, to, cell, n = 25) {
  grid <- softplus_inv(exp(seq(log(natpar(modes[[from]])$tau[cell]),
                               log(natpar(modes[[to]])$tau[cell]), length.out = n)))
  nlp <- function(v8, eta) { r <- -lp_fun(to_u(v8, eta, cell), prior); if (!is.finite(r)) 1e10 else r }
  sweep <- function(u0, gr) {
    v <- from_u(u0, cell)
    t(vapply(gr, function(e) {
      o <- stats::optim(v, nlp, eta = e, method = "BFGS",
                        control = list(maxit = 2000, reltol = 1e-13, ndeps = rep(1e-5, 8)))
      v <<- o$par
      n <- natpar(to_u(o$par, e, cell))
      c(-o$value, n$mu[cell]^2 - 2 / n$tau[cell])
    }, numeric(2)))
  }
  f <- sweep(modes[[from]], grid)
  b <- sweep(modes[[to]], rev(grid))[n:1, ]
  use_f <- f[, 1] >= b[, 1]
  data.frame(prior = prior, path = paste(from, "->", to), cell = c("Accuracy", "Speed")[cell],
             tau = softplus(grid), lp = ifelse(use_f, f[, 1], b[, 1]),
             k2 = ifelse(use_f, f[, 2], b[, 2]),
             lp_from = lp_fun(modes[[from]], prior), lp_to = lp_fun(modes[[to]], prior))
}
res <- rbind(profile("WW", "EW", 1), profile("WW", "WE", 2))
for (p in unique(res$path)) {
  z <- res[res$path == p, ]
  j <- which.min(z$lp)
  cat(sprintf("%s (%s tau): lp %.2f -> %.2f, profile minimum %.2f at tau %.3f (k2 %+.1f): barrier >= %.1f / %.1f\n",
              p, z$cell[1], z$lp_from[1], z$lp_to[1], z$lp[j], z$tau[j], z$k2[j],
              z$lp_from[1] - z$lp[j], z$lp_to[1] - z$lp[j]))
}
num <- vapply(res, is.numeric, logical(1))
res[num] <- lapply(res[num], signif, 7)
utils::write.csv(res, file.path(OUT, sprintf("barrier_%s.csv", tag)), row.names = FALSE)
