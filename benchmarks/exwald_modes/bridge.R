# The posterior mass in each mode, from HMC draws: bridge sampling (normal
# and warp3, 3 repetitions each) on the post-warmup draws of every chain that
# sat in the mode, against common.R's lp__ twin. The modes are far enough
# apart that each mode's proposal has no mass in the others, so each estimate
# is that mode's share of the normalising constant.
#
#   Rscript benchmarks/exwald_modes/bridge.R <prior> <run dir> [<run dir> ...]
#
# Writes results/bridge_<prior>.csv.
source("benchmarks/exwald_modes/common.R")
args <- commandArgs(trailingOnly = TRUE)
prior <- args[1]
U <- NULL
lab <- NULL
for (r in args[-1]) {
  x <- readRDS(file.path(r, "draws.rds"))
  stopifnot(identical(x$prior, prior))
  dr <- x$draws[x$draws$.iteration > x$warmup, ]
  u <- draws_to_u(dr)
  U <- rbind(U, u)
  lab <- c(lab, apply(u, 1, mode_label))
}
set.seed(2026)
cat("draws per mode:\n")
print(table(lab))
lpost <- function(pars, data) lp_fun(as.numeric(pars), prior)
res <- NULL
for (m in sort(unique(lab))) {
  S <- U[lab == m, , drop = FALSE]
  for (meth in c("normal", "warp3")) {
    bs <- bridgesampling::bridge_sampler(
      samples = S, log_posterior = lpost, data = NULL, method = meth, silent = TRUE,
      lb = stats::setNames(rep(-Inf, 9), PN), ub = stats::setNames(rep(Inf, 9), PN),
      repetitions = 3)
    res <- rbind(res, data.frame(prior = prior, mode = m, method = meth, draws = nrow(S),
                                 logml = mean(bs$logml), logml_sd = stats::sd(bs$logml)))
  }
}
res$share <- NA
for (meth in unique(res$method)) {
  k <- res$method == meth
  res$share[k] <- exp(res$logml[k] - max(res$logml[k])) / sum(exp(res$logml[k] - max(res$logml[k])))
}
num <- vapply(res, is.numeric, logical(1))
res[num] <- lapply(res[num], signif, 7)
print(res, row.names = FALSE)
utils::write.csv(res, file.path(OUT, sprintf("bridge_%s.csv", gsub("[^0-9A-Za-z]", "", prior))),
                 row.names = FALSE)
