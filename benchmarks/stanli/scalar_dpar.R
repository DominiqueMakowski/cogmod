# What an intercept-only dpar costs as brms writes it. For `sigma ~ 1` brms
# builds vector[N] sigma = rep_vector(0.0, N), adds the intercept and applies
# the inverse link to all N elements; the same program with sigma declared
# `real` computes it once, and a vectorised lpdf then gets a scalar. Measured
# on a built-in family here (bench.R's `orig_s` does the same for the LNR).
#
#   Rscript benchmarks/stanli/scalar_dpar.R [--out DIR] [--reps 21]
#
# 2026-10-04, Windows laptop (i7-1265U), CmdStan 2.38, brms 2.23.1, N = 5000:
#   log density: same to 1.2e-14 relative (identical at 14 of 200 points; not
#   bitwise, as normal_lpdf() takes N * log(sigma) once instead of summing N
#   logs); gradient to 5.9e-14
#   us per gradient: vector 620, real 240; real / vector 0.41 [0.32, 0.47]
source("benchmarks/gradient_programs.R")  # gp_args(), gp_compile(), gp_methods()
opt <- gp_args(list(out = "benchmarks/results/stanli/scalar_dpar", reps = 21L))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
suppressPackageStartupMessages(library(brms))

set.seed(1)
N <- 5000
d <- data.frame(x = rnorm(N))
d$y <- 1 + 0.5 * d$x + rnorm(N, 0, 0.8)
f <- bf(y ~ x, sigma ~ 1)
code <- as.character(make_stancode(f, data = d))
old <- "vector[N] sigma = rep_vector(0.0, N);"
stopifnot(lengths(regmatches(code, gregexpr(old, code, fixed = TRUE))) == 1)
writeLines(code, file.path(opt$out, "vector.stan"))
writeLines(sub(old, "real sigma = 0;", code, fixed = TRUE), file.path(opt$out, "real.stan"))
cmdstanr::write_stan_json(lapply(unclass(make_standata(f, data = d)), identity),
                          file.path(opt$out, "data.json"))
fit <- lapply(c(vector = "vector", real = "real"), function(v)
  gp_methods(gp_compile(file.path(opt$out, paste0(v, ".stan"))), file.path(opt$out, "data.json"))$fit)

set.seed(2)
P <- lapply(1:500, function(k) rnorm(3, 0, 0.5))  # b, Intercept, Intercept_sigma
g <- lapply(fit, function(x) lapply(P[1:200], x$grad_log_prob))
lp <- sapply(g, function(x) vapply(x, attr, 0, "log_prob"))
cat(sprintf("log density: max relative difference %.2g (identical at %d of %d points)\n",
            max(abs(lp[, 1] - lp[, 2]) / abs(lp[, 1])), sum(lp[, 1] == lp[, 2]), nrow(lp)))
cat(sprintf("gradient: max relative difference %.2g\n", max(mapply(function(a, b)
  max(abs(a - b) / pmax(1, abs(a))), g$vector, g$real))))

block <- function(x) system.time(for (q in P) x$grad_log_prob(q))[["elapsed"]] / length(P) * 1e6
for (x in fit) block(x)  # warm up
r <- t(replicate(opt$reps, { us <- c(vector = NA, real = NA)
  for (v in sample(names(fit))) us[[v]] <- block(fit[[v]]); us }))
q <- stats::quantile(r[, "real"] / r[, "vector"], c(0.1, 0.5, 0.9))
cat(sprintf("us per gradient: vector %.0f, real %.0f; real / vector %.2f [%.2f, %.2f]\n",
            median(r[, "vector"]), median(r[, "real"]), q[2], q[1], q[3]))
