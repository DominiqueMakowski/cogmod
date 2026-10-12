# What the W and E solutions of a cell disagree about: per-trial log-likelihood
# differences between the two modes, summed within RT bins, and the RT
# quantiles each predicts against the observed ones. Speed is compared WW vs
# WE and Accuracy WW vs EW, so the other cell is held at its W solution.
#
#   Rscript benchmarks/exwald_modes/shape.R <prior>     (needs modes.R's .rds)
#
# Writes results/shape_bins.csv and results/shape_quantiles.csv.
source("benchmarks/exwald_modes/common.R")
prior <- commandArgs(trailingOnly = TRUE)[1]
md <- readRDS(file.path(OUT, sprintf("modes_%s.rds", gsub("[^0-9A-Za-z]", "", prior))))
modes <- stats::setNames(md$par, md$label)
bins <- NULL
qs <- NULL
probs <- c(0.001, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 0.999)
for (k in 1:2) {
  cond <- c("Accuracy", "Speed")[k]
  y <- Y[D$Condition == cond]
  pw <- natpar(modes$WW)[k, ]
  pe <- natpar(modes[[if (k == 1) "EW" else "WE"]])[k, ]
  lw <- dcogmod_exwald(y, pw$mu, pw$boundary, pw$tau, pw$ndt, pw$poutlier, log = TRUE)
  le <- dcogmod_exwald(y, pe$mu, pe$boundary, pe$tau, pe$ndt, pe$poutlier, log = TRUE)
  br <- stats::quantile(y, c(0, 0.01, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.99, 1))
  bin <- cut(y, br, include.lowest = TRUE)
  bins <- rbind(bins, data.frame(cell = cond, rt_bin = levels(bin), trials = as.vector(table(bin)),
                                 E_minus_W = as.vector(tapply(le - lw, bin, sum))))
  set.seed(1)
  sw <- rcogmod_exwald(1e6, pw$mu, pw$boundary, pw$tau, pw$ndt, 0)
  se <- rcogmod_exwald(1e6, pe$mu, pe$boundary, pe$tau, pe$ndt, 0)
  qs <- rbind(qs, data.frame(cell = cond, prob = probs, observed = stats::quantile(y, probs, names = FALSE),
                             W = stats::quantile(sw, probs, names = FALSE),
                             E = stats::quantile(se, probs, names = FALSE)))
  cat(sprintf("%s: E - W log-lik %.2f; tail decay rates (1/s) W: Wald %.2f, exp %.2f | E: Wald %.2f, exp %.2f\n",
              cond, sum(le - lw), pw$mu^2 / 2, 1 / pw$tau, pe$mu^2 / 2, 1 / pe$tau))
}
bins$E_minus_W <- round(bins$E_minus_W, 3)
qs[, 3:5] <- round(qs[, 3:5], 4)
print(bins, row.names = FALSE)
print(qs, row.names = FALSE)
utils::write.csv(bins, file.path(OUT, "shape_bins.csv"), row.names = FALSE)
utils::write.csv(qs, file.path(OUT, "shape_quantiles.csv"), row.names = FALSE)
