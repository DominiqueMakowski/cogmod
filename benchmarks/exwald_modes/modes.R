# The posterior's local modes under one prior set: optimisation of lp__
# (common.R's R twin) from every cogmod_inits() start a 4-chain fit with seeds
# 1-4 would use, plus one hand-placed start in each of the four W/E
# combinations (from those 16 alone, 0.3.5 and P1 never reach EW), distinct
# optima kept, each with its Hessian and Laplace log-mass. Writes
# results/modes_<prior>.csv and .rds (the vectors, for fit.R's modes: inits).
#
#   Rscript benchmarks/exwald_modes/modes.R <prior>
source("benchmarks/exwald_modes/common.R")
prior <- commandArgs(trailingOnly = TRUE)[1]
nlp <- function(u) { v <- -lp_fun(u, prior); if (!is.finite(v)) 1e10 else v }
opt <- function(u0) {
  o <- stats::optim(u0, nlp, method = "Nelder-Mead", control = list(maxit = 3000, reltol = 1e-10))
  stats::optim(o$par, nlp, method = "BFGS",
               control = list(maxit = 3000, reltol = 1e-15, ndeps = rep(1e-5, 9)))
}
starts <- list()
for (s in 1:4) {
  set.seed(s)
  ini <- suppressWarnings(cogmod_inits(FORMULA, D))
  invisible(ini(1))
  for (ch in 1:4) {
    i <- ini(ch)
    starts[[length(starts) + 1]] <- c(i$b, i$Intercept, i$b_boundary, i$Intercept_boundary,
                                      i$b_tau, i$Intercept_tau, i$b_ndt, i$Intercept_ndt,
                                      stats::qlogis(i$poutlier))
  }
}
# Hand-placed starts, natural values per cell: W is a skewed Wald with a short
# exponential stage, E a near-Gaussian Wald with a long one and a lower ndt
cell <- list(W = list(Accuracy = c(2.5, 0.73, 0.06, 0.33), Speed = c(3.6, 0.74, 0.037, 0.306)),
             E = list(Accuracy = c(7, 1.5, 0.21, 0.27), Speed = c(10, 2, 0.12, 0.22)))
for (combo in c("WW", "WE", "EW", "EE")) {
  a <- cell[[substr(combo, 1, 1)]]$Accuracy
  s <- cell[[substr(combo, 2, 2)]]$Speed
  pair <- function(eA, eS) { b <- eS - eA; c(b, eA + MX * b) }
  starts[[length(starts) + 1]] <- c(pair(softplus_inv(a[1]), softplus_inv(s[1])),
                                    pair(softplus_inv(a[2]), softplus_inv(s[2])),
                                    pair(softplus_inv(a[3]), softplus_inv(s[3])),
                                    pair(log(a[4]), log(s[4])), stats::qlogis(0.0018))
}
res <- lapply(starts, function(u0) { o <- opt(u0); list(par = o$par, lp = -o$value) })
lab <- vapply(res, function(r) mode_label(r$par), character(1))
cat("modes reached from the 16 cogmod_inits() starts:\n")
print(table(lab[1:16]))
cat("and from the hand-placed WW, WE, EW, EE starts:", lab[17:20], "\n")
keep <- !duplicated(lab)
modes <- res[keep]
labs <- lab[keep]
tab <- do.call(rbind, lapply(seq_along(modes), function(k) {
  u <- modes[[k]]$par
  H <- stats::optimHess(u, nlp, control = list(ndeps = rep(1e-4, 9)))
  ev <- eigen(H, symmetric = TRUE, only.values = TRUE)$values
  n <- natpar(u)
  data.frame(prior = prior, mode = labs[k], from_inits = sum(lab[1:16] == labs[k]),
             lp = modes[[k]]$lp, loglik = loglik_fun(u), lprior = lprior_fun(u, prior),
             min_eig = min(ev), laplace_logmass = modes[[k]]$lp - 0.5 * sum(log(ev)) + 4.5 * log(2 * pi),
             muA = n$mu[1], aA = n$boundary[1], tauA = n$tau[1], ndtA = n$ndt[1],
             muS = n$mu[2], aS = n$boundary[2], tauS = n$tau[2], ndtS = n$ndt[2],
             poutlier = n$poutlier[1], k2A = n$mu[1]^2 - 2 / n$tau[1], k2S = n$mu[2]^2 - 2 / n$tau[2])
}))
tab$laplace_share <- exp(tab$laplace_logmass - max(tab$laplace_logmass)) /
  sum(exp(tab$laplace_logmass - max(tab$laplace_logmass)))
tab <- tab[order(-tab$lp), ]
num <- vapply(tab, is.numeric, logical(1))
tab[num] <- lapply(tab[num], signif, 6)
print(tab, row.names = FALSE)
tag <- gsub("[^0-9A-Za-z]", "", prior)
utils::write.csv(tab, file.path(OUT, sprintf("modes_%s.csv", tag)), row.names = FALSE)
saveRDS(list(label = labs, par = lapply(modes, `[[`, "par")),
        file.path(OUT, sprintf("modes_%s.rds", tag)))
