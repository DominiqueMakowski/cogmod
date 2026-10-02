# Self-contained reproducer for the stanli issue (README.md, "Upstream"):
# nothing from cogmod, only stanli. Run with Rscript from anywhere.
#
# 1. log_mix() inside control flow that depends on a parameter fails to
#    compile; log_sum_exp() written out in the same place compiles.
# 2. A comparison used as a value (`log(y > mu)`) fails to compile, branch or
#    not; the same comparison in a ternary or an `if` compiles.
# 3. A per-observation function with an `if` on a parameter costs far more
#    per gradient than the same arithmetic without it, even when the branch
#    is never taken. This is what a package of custom brms families pays:
#    numerically careful densities branch (tail switches, support checks).
# 4. On the full LNR program (bench.R), log_sum_exp() in place of log_mix()
#    cost 5x and an early return on data alone 4x. Here, alone, they cost
#    1.2x and nothing (results/repro_issue.txt): the cliff needs the rest of
#    the program, so only the emitted LNR programs reproduce it.
suppressPackageStartupMessages(library(stanli))
cat("stanli", as.character(packageVersion("stanli")), "\n\n")

program <- function(body) sprintf('
functions {
  real f(real y, real mu) {
%s
  }
}
data { int N; vector[N] y; }
parameters { real mu; real<lower=0> sigma; }
model {
  mu ~ normal(0, 1);
  sigma ~ normal(0, 1);
  for (n in 1:N) target += f(y[n], mu) + normal_lpdf(y[n] | mu, sigma);
}', body)

try_build <- function(label, body, data) {
  r <- tryCatch({ stanli_model(code = program(body), data = data); "compiles" },
                error = function(e) paste("FAILS:", conditionMessage(e)))
  cat(sprintf("  %-34s %s\n", label, r))
}

set.seed(1)
small <- list(N = 3, y = c(0.5, 1, 2))

cat("1. log_mix() in a parameter-dependent branch\n")
try_build("log_mix(), no branch",
  "    return log_mix(0.1, -1.0, normal_lpdf(y | mu, 1));", small)
try_build("log_mix(), inside `if (y - mu <= 0)`",
  "    if (y - mu <= 0) return log(0.1);\n    return log_mix(0.1, -1.0, normal_lpdf(y | mu, 1));", small)
try_build("log_sum_exp() written out, same `if`",
  "    if (y - mu <= 0) return log(0.1);\n    return log_sum_exp(log(0.1) - 1.0, log1m(0.1) + normal_lpdf(y | mu, 1));", small)

cat("\n2. A comparison as a value\n")
try_build("log(y > mu)", "    return log(y > mu);", small)
try_build("y > mu ? 0 : -1e3", "    return y > mu ? 0 : -1e3;", small)
try_build("step(y - mu)", "    return step(y - mu);", small)

cat("\n3. Cost of a parameter-dependent branch, N = 5000, per gradient\n")
big <- list(N = 5000, y = rnorm(5000, 1, 1))
# log Phi(y - mu) on the erfc route; the branched copy adds the kind of
# far-tail switch a careful implementation needs. Every y here is within a
# few units of mu, so the branch is never taken: same value, same gradient.
bodies <- c(
  no_branch = "    return log(0.5 * erfc(-(y - mu) * 0.7071067811865476));",
  branch = paste0("    if (y - mu < -25) return -0.5 * square(y - mu) - log(mu - y) - 0.9189385332046727;\n",
                  "    return log(0.5 * erfc(-(y - mu) * 0.7071067811865476));")
)
models <- lapply(bodies, function(b) stanli_model(code = program(b), data = big))
q <- c(0.9, log(1.1))
g <- lapply(models, log_prob_grad, q = q)
cat(sprintf("  same lp: %s (%.10f vs %.10f); same gradient: %s\n",
            isTRUE(all.equal(g[[1]]$lp, g[[2]]$lp)), g[[1]]$lp, g[[2]]$lp,
            isTRUE(all.equal(g[[1]]$grad, g[[2]]$grad))))
K <- 200
Q <- lapply(seq_len(K), function(k) q + rnorm(2, 0, 0.05))
for (m in models) for (k in 1:20) log_prob_grad(m, Q[[k]])  # warm up
tm <- replicate(11, sapply(models, function(m)
  system.time(for (k in seq_len(K)) log_prob_grad(m, Q[[k]]))[["elapsed"]] / K * 1e6))
med <- apply(tm, 1, stats::median)
cat(sprintf("  no branch %.0f us, branch %.0f us: %.1fx\n",
            med[["no_branch"]], med[["branch"]], med[["branch"]] / med[["no_branch"]]))

# 4. The two costs seen on the full LNR program, tried alone.
cat("\n4. Same value, no parameter-dependent branch, N = 5000, per gradient\n")
pos <- list(N = 5000, y = abs(rnorm(5000, 1, 1)) + 0.01)
mix_bodies <- c(
  log_mix = "    return log_mix(0.1, -12.5 * square(y), normal_lpdf(y | mu, 1));",
  log_sum_exp = "    return log_sum_exp(log(0.1) - 12.5 * square(y), log1m(0.1) + normal_lpdf(y | mu, 1));",
  log_mix_data_check = paste0("    if (y <= 0) return negative_infinity();\n",
                              "    return log_mix(0.1, -12.5 * square(y), normal_lpdf(y | mu, 1));")
)
mm <- lapply(mix_bodies, function(b) stanli_model(code = program(b), data = pos))
lp <- sapply(mm, function(m) log_prob_grad(m, q)$lp)
cat(sprintf("  same lp: %s\n", isTRUE(all.equal(unname(lp), rep(lp[[1]], 3)))))
for (m in mm) for (k in 1:20) log_prob_grad(m, Q[[k]])
tm <- replicate(11, sapply(mm, function(m)
  system.time(for (k in seq_len(K)) log_prob_grad(m, Q[[k]]))[["elapsed"]] / K * 1e6))
med <- apply(tm, 1, stats::median)
for (nm in names(med)) cat(sprintf("  %-20s %5.0f us  (%.1fx)\n", nm, med[[nm]], med[[nm]] / med[["log_mix"]]))
