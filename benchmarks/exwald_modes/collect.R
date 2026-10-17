# Per-chain summaries of fit.R runs -> results/chains.csv, and the per-run
# R-hat / bulk ESS -> results/runs.csv.
#
#   Rscript benchmarks/exwald_modes/collect.R <run dir> [<run dir> ...]
#
# A run dir holds the draws.rds fit.R writes. The mode of a chain is the label
# of its post-warmup draws (common.R's mode_label(), W = Wald tail,
# E = exponential tail, Accuracy first); `capture_iter` is the first warmup
# iteration from which the chain's label stayed what it ended as.
source("benchmarks/exwald_modes/common.R")
dirs <- commandArgs(trailingOnly = TRUE)
chains <- NULL
runs <- NULL
for (r in dirs) {
  x <- readRDS(file.path(r, "draws.rds"))
  W <- x$warmup
  dr <- x$draws
  sdg <- x$sdg
  u <- draws_to_u(dr)
  lab <- apply(u, 1, mode_label)
  for (k in sort(unique(dr$.chain))) {
    sel <- dr$.chain == k
    post <- sel & dr$.iteration > W
    g <- sdg[sdg$.chain == k, ]
    lk <- lab[sel]
    final <- names(which.max(table(lab[post])))
    changed <- which(lk != final)
    n <- natpar(apply(u[post, , drop = FALSE], 2, stats::median))
    ini <- x$inits[[k]]
    iu <- c(ini$b, ini$Intercept, ini$b_boundary, ini$Intercept_boundary, ini$b_tau,
            ini$Intercept_tau, ini$b_ndt, ini$Intercept_ndt, stats::qlogis(ini$poutlier))
    ni <- natpar(iu)
    chains <- rbind(chains, data.frame(
      run = basename(r), prior = x$prior, chain = k, origin = x$origin[k],
      mode = final, mode_share = mean(lab[post] == final),
      capture_iter = if (length(changed)) max(changed) + 1 else 1,
      lp_mean = mean(dr$lp__[post]),
      leap_warm = sum(g$n_leapfrog__[g$.iteration <= W]),
      leap_samp = sum(g$n_leapfrog__[g$.iteration > W]),
      div = sum(g$divergent__[g$.iteration > W]),
      stepsize = g$stepsize__[g$.iteration > W][1],
      muA = n$mu[1], aA = n$boundary[1], tauA = n$tau[1], ndtA = n$ndt[1],
      muS = n$mu[2], aS = n$boundary[2], tauS = n$tau[2], ndtS = n$ndt[2],
      init_tauA = ni$tau[1], init_k2A = ni$mu[1]^2 - 2 / ni$tau[1],
      init_ndtA = ni$ndt[1], minutes = x$time$chains$total[k] / 60))
  }
  post <- dr[dr$.iteration > W, ]
  pd <- posterior::as_draws_df(post[, c("b[1]", "Intercept", "b_boundary[1]", "Intercept_boundary",
                                        "b_tau[1]", "Intercept_tau", "b_ndt[1]", "Intercept_ndt",
                                        "poutlier", "lp__", ".chain", ".iteration", ".draw")])
  s <- posterior::summarise_draws(pd, "rhat", "ess_bulk")
  runs <- rbind(runs, data.frame(run = basename(r), prior = x$prior, chains = length(x$origin),
                                 max_rhat = max(s$rhat), min_ess_bulk = min(s$ess_bulk)))
}
num <- vapply(chains, is.numeric, logical(1))
chains[num] <- lapply(chains[num], signif, 6)
utils::write.csv(chains, file.path(OUT, "chains.csv"), row.names = FALSE)
utils::write.csv(runs, file.path(OUT, "runs.csv"), row.names = FALSE)
print(chains[, c("run", "chain", "origin", "mode", "capture_iter", "lp_mean", "leap_warm", "leap_samp")],
      row.names = FALSE)
print(runs, row.names = FALSE)
