# Is the W/E bimodality a shortage of data, or a property of the model?
#
# One condition at a time (the four-parameter ex-Wald plus poutlier, 0.3.5's
# intercept priors), on two kinds of data:
#
#  - simulated from the ex-Wald itself, at each of the four cell solutions the
#    real fit found (Accuracy W, Accuracy E, Speed W, Speed E), with 500 to
#    32,000 trials. The truth is then known, and an ex-Wald.
#  - the real trials of each condition, subsampled to 25% and 50%.
#
# For each data set: optimise lp from a W start and from an E start, label
# what each reached, and if they reached different optima report the gap, the
# Laplace share of the mass and, for the first replicate, the barrier along
# `tau` between them. If the second optimum were a shortage of data it would
# disappear as the trials grow; if it is a property of the model it stays,
# and its gap and barrier grow with the trials.
#
#   Rscript benchmarks/exwald_modes/data_size.R <job>
#
# <job> is one of Accuracy:W, Accuracy:E, Speed:W, Speed:E (simulated, truth
# = that cell solution), real (the subsamples), or combine, which stacks the
# jobs' CSVs into results/data_size.csv. The first five can run in parallel.
source("benchmarks/exwald_modes/common.R")
job <- commandArgs(trailingOnly = TRUE)[1]
if (identical(job, "combine")) {
  fs <- list.files(OUT, pattern = "^data_size_.+\\.csv$", full.names = TRUE)
  x <- do.call(rbind, lapply(fs, utils::read.csv))
  utils::write.csv(x, file.path(OUT, "data_size.csv"), row.names = FALSE)
  unlink(fs)
  quit(save = "no")
}

# lp of one condition on the unconstrained scale: (eta_mu, eta_boundary,
# eta_tau, log ndt, logit poutlier), 0.3.5's priors on the intercepts
lp1 <- function(u, y) {
  p <- stats::plogis(u[5])
  ll <- sum(dcogmod_exwald(y, mu = softplus(u[1]), boundary = softplus(u[2]), tau = softplus(u[3]),
                           ndt = exp(u[4]), poutlier = p, log = TRUE))
  ll + stats::dt(u[1] / 2.5, 3, log = TRUE) - log(2.5) + stats::dnorm(u[2], 0.5, 1, log = TRUE) +
    stats::dnorm(u[3], -1.5, 0.7, log = TRUE) + stats::dnorm(u[4], -1.2, 0.5, log = TRUE) +
    stats::dexp(p, 280, log = TRUE) - stats::pexp(1, 280, log.p = TRUE) + log(p) + log1p(-p)
}
to_u1 <- function(x, p = 0.002) {
  c(softplus_inv(x[1]), softplus_inv(x[2]), softplus_inv(x[3]), log(x[4]), stats::qlogis(p))
}
label1 <- function(u) if (softplus(u[1])^2 - 2 / softplus(u[3]) > 0) "E" else "W"
opt1 <- function(u0, y) {
  nlp <- function(u) { v <- -lp1(u, y); if (!is.finite(v)) 1e10 else v }
  o <- stats::optim(u0, nlp, method = "BFGS",
                    control = list(maxit = 3000, reltol = 1e-14, ndeps = rep(1e-6, 5)))
  ev <- eigen(stats::optimHess(o$par, nlp, control = list(ndeps = rep(1e-5, 5))),
              symmetric = TRUE, only.values = TRUE)$values
  list(par = o$par, lp = -o$value, type = label1(o$par), min_eig = min(ev),
       logmass = -o$value - 0.5 * sum(log(pmax(ev, 1e-12))) + 2.5 * log(2 * pi))
}
# lp profiled along tau between two optima: the floor bounds the saddle
barrier1 <- function(uW, uE, y, n_grid = 13) {
  grid <- seq(uW[3], uE[3], length.out = n_grid)
  sweep <- function(u0, gr) {
    v <- u0[-3]
    vapply(gr, function(e) {
      nlp <- function(w) { r <- -lp1(c(w[1:2], e, w[3:4]), y); if (!is.finite(r)) 1e10 else r }
      o <- stats::optim(v, nlp, method = "BFGS",
                        control = list(maxit = 2000, reltol = 1e-13, ndeps = rep(1e-6, 4)))
      v <<- o$par
      -o$value
    }, numeric(1))
  }
  pmax(sweep(uW, grid), rev(sweep(uE, rev(grid))))
}

# The W and E solutions of the real fit, cell by cell (0.3.5's modes): truths
# for the simulations, and the two starts for every data set
md <- readRDS(file.path(OUT, "modes_035.rds"))
cells <- list()
for (k in seq_along(md$label)) {
  n <- natpar(md$par[[k]])
  for (j in 1:2) {
    key <- paste(n$cond[j], substr(md$label[k], j, j))
    if (is.null(cells[[key]])) cells[[key]] <- c(n$mu[j], n$boundary[j], n$tau[j], n$ndt[j])
  }
}

one <- function(y, cond, what, truth, rep) {
  oW <- opt1(to_u1(cells[[paste(cond, "W")]]), y)
  oE <- opt1(to_u1(cells[[paste(cond, "E")]]), y)
  two <- oW$type == "W" && oE$type == "E"
  bar <- if (two && rep == 1) min(barrier1(oW$par, oE$par, y)) else NA
  r <- data.frame(data = what, cell = cond, truth = truth, trials = length(y), rep = rep,
                  from_W_start = oW$type, from_E_start = oE$type, two_modes = two,
                  lp_W = oW$lp, lp_E = if (two) oE$lp else NA,
                  lp_W_minus_E = if (two) oW$lp - oE$lp else NA,
                  share_W = if (two) 1 / (1 + exp(oE$logmass - oW$logmass)) else NA,
                  barrier_from_W = oW$lp - bar, barrier_from_E = if (two) oE$lp - bar else NA,
                  min_eig_W = oW$min_eig, min_eig_E = if (two) oE$min_eig else NA)
  cat(sprintf("%s %s truth %s n=%d rep %d: %s/%s  W-E %.2f  barrier %.1f / %.1f\n", what, cond, truth,
              length(y), rep, oW$type, oE$type, r$lp_W_minus_E, r$barrier_from_W, r$barrier_from_E))
  r
}

out <- NULL
if (job == "real") {
  for (cond in c("Accuracy", "Speed")) {
    yall <- Y[D$Condition == cond]
    for (frac in c(0.25, 0.5, 1)) for (rep in if (frac < 1) 1:3 else 1) {
      set.seed(rep * 100 + round(frac * 100))
      y <- if (frac < 1) sample(yall, round(frac * length(yall))) else yall
      out <- rbind(out, one(y, cond, sprintf("real %d%%", round(100 * frac)), "?", rep))
    }
  }
} else {
  cond <- sub(":.*", "", job)
  truth <- sub(".*:", "", job)
  th <- cells[[paste(cond, truth)]]
  for (n in c(500, 2000, 8000, 32000)) for (rep in 1:3) {
    set.seed(1e5 * (cond == "Speed") + 1e4 * (truth == "E") + 1000 * rep + n %% 997)
    y <- rcogmod_exwald(n, mu = th[1], boundary = th[2], tau = th[3], ndt = th[4], poutlier = 0.0015)
    out <- rbind(out, one(y, cond, "simulated", truth, rep))
  }
}
num <- vapply(out, is.numeric, logical(1))
out[num] <- lapply(out[num], signif, 6)
utils::write.csv(out, file.path(OUT, sprintf("data_size_%s.csv", gsub("[^A-Za-z]", "", job))),
                 row.names = FALSE)
