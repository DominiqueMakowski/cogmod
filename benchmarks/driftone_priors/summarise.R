# The tables of README.md, from the stacked results collect.R writes. Run from
# the package root:
#
#   Rscript benchmarks/driftone_priors/summarise.R [model ...]
#
# Per model and arm: sampler health summed over seeds (divergences after
# warmup, treedepth hits, chains with a collapsed step size, R-hat, the chain
# spread - the range of the chains' medians in pooled SDs, which is where a
# split shows first), cost (leapfrog steps per sampling draw) and efficiency
# (bulk ESS per 1,000 leapfrog steps, per parameter); the posterior medians
# with the prior quantile of each where the prior is a normal; and, for the
# race families, the error rate and median RTs the decision process predicts.
args <- commandArgs(trailingOnly = TRUE)
res <- file.path("benchmarks/driftone_priors", "results")
chains <- utils::read.csv(file.path(res, "chains.csv"))
pars <- utils::read.csv(file.path(res, "pars.csv"), check.names = FALSE)
pp <- utils::read.csv(file.path(res, "pp.csv"))
priors <- utils::read.csv(file.path(res, "priors.csv"))
models <- if (length(args)) args else unique(pars$model)
options(width = 200)

leap <- stats::aggregate(cbind(leap = leap_warm + leap_samp, samp = leap_samp) ~ model + arm + seed,
                         chains, sum)
leap$per_draw <- leap$samp / (4 * 250)

prior_q <- function(model, arm, variable, x) {
  dp <- cls <- NULL
  if (grepl("^Intercept_", variable)) { dp <- sub("^Intercept_", "", variable); cls <- "Intercept" }
  if (grepl("^b_[a-z]+\\[", variable)) { dp <- sub("^b_([a-z]+)\\[.*", "\\1", variable); cls <- "b" }
  if (is.null(dp)) return(NA_real_)
  r <- priors[priors$model == model & priors$arm == arm & priors$dpar == dp &
                priors$class == cls & (is.na(priors$coef) | !nzchar(priors$coef)), ]
  m <- regmatches(r$prior[1], regexec("^normal\\(([-0-9.]+), ([0-9.]+)\\)$", r$prior[1]))[[1]]
  if (length(m) != 3) return(NA_real_)
  stats::pnorm(x, as.numeric(m[2]), as.numeric(m[3]))
}

for (mdl in models) {
  cat("\n================", mdl, "================\n")
  ch <- chains[chains$model == mdl, ]
  pa <- pars[pars$model == mdl, ]
  health <- do.call(rbind, lapply(split(ch, ch$arm), function(x) {
    a <- x$arm[1]
    p <- pa[pa$arm == a, ]
    l <- leap[leap$model == mdl & leap$arm == a, ]
    data.frame(arm = a, seeds = length(unique(x$seed)), chains = nrow(x),
               div = sum(x$div), td10 = sum(x$td10), stuck = sum(x$stepsize < 1e-3),
               max_rhat = round(max(p$rhat, na.rm = TRUE), 3),
               max_chain_spread = round(max(p$chain_spread, na.rm = TRUE), 2),
               leapfrog_per_draw = round(stats::median(l$per_draw), 1))
  }))
  print(health, row.names = FALSE)

  keep <- !grepl("^(Intercept$|b\\[[2-9]|b_Intercept$)", pa$variable)
  pa <- pa[keep, ]
  pa <- merge(pa, leap[leap$model == mdl, c("arm", "seed", "leap")], by = c("arm", "seed"))
  tab <- do.call(rbind, lapply(split(pa, list(pa$variable, pa$arm), drop = TRUE), function(y) {
    data.frame(variable = y$variable[1], arm = y$arm[1],
               median = signif(stats::median(y$median), 4),
               q5 = signif(stats::median(y[["5%"]]), 4), q95 = signif(stats::median(y[["95%"]]), 4),
               ess_bulk = round(stats::median(y$ess_bulk)),
               ess_per_kleap = round(stats::median(y$ess_bulk / (y$leap / 1000)), 2))
  }))
  tab$prior_q <- signif(mapply(prior_q, mdl, tab$arm, tab$variable, tab$median), 2)
  print(tab[order(tab$variable, tab$arm), ], row.names = FALSE)

  p <- pp[pp$model == mdl, ]
  if (nrow(p)) {
    agg <- stats::aggregate(cbind(median, q5, q95) ~ cond + stat + arm, p, stats::median)
    agg[, c("median", "q5", "q95")] <- signif(agg[, c("median", "q5", "q95")], 3)
    cat("\npredicted by the decision process (median over seeds):\n")
    print(agg[order(agg$cond, agg$stat, agg$arm), ], row.names = FALSE)
  }
}
