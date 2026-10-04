# stanli 0.19.0 runs chains one at a time whenever sample_model() is given
# `init`, whatever parallel_chains says. Self-contained (stanli only); run
# with Rscript from anywhere. Found 2026-10-03 when bench.R's fit started its
# stanli arms from cogmod_inits(): wall time came out as the sum of the chain
# times, on Windows and on Linux. sample_model() hands init_vec and
# parallel_chains to the native sampler unchanged, so the switch is there.
#
# Output on the Windows laptop (i7-1265U):
#   init none   wall  2.49 s, chains max  2.48, sum  7.59: wall/sum 0.33
#   init vector wall  8.11 s, chains max  2.59, sum  8.09: wall/sum 1.00
#   init matrix wall  8.73 s, chains max  3.42, sum  8.70: wall/sum 1.00
suppressPackageStartupMessages(library(stanli))
cat("stanli", as.character(packageVersion("stanli")), "\n")
set.seed(1)
N <- 5000
d <- list(N = N, y = abs(rnorm(N, 1, 1)) + 0.05)
code <- "data { int N; vector[N] y; }
parameters { real mu; real<lower=0> sigma; real<lower=0, upper=1> p; }
model { mu ~ normal(0, 1); sigma ~ normal(0, 1); p ~ beta(1, 20);
  for (n in 1:N) target += log_mix(p, -12.5 * square(y[n]), normal_lpdf(y[n] | mu, sigma)); }"
m <- stanli_model(code = code, data = d)
u <- unconstrain(m, list(mu = 1, sigma = 1, p = 0.05))
inits <- list(none = NULL, vector = u, matrix = rbind(u, u + 0.1, u - 0.1, u + 0.2))
for (nm in names(inits)) {
  t0 <- Sys.time()
  f <- sample_model(m, chains = 4, seed = 1, warmup = 150, samples = 150, parallel_chains = 4,
                    refresh = 0, init = inits[[nm]])
  w <- as.numeric(Sys.time() - t0, units = "secs")
  ch <- f$report$warmup_seconds + f$report$sampling_seconds
  cat(sprintf("init %-6s wall %5.2f s, chains max %5.2f, sum %5.2f: wall/sum %.2f\n",
              nm, w, max(ch), sum(ch), w / sum(ch)))
}
