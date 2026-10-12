# Does stanli change the R process's CPU affinity, which processes R starts
# afterwards (CmdStan's chains, under cmdstanr) inherit? Prints the allowed
# CPU list of R itself and of a fresh child before stanli is loaded, after
# library(stanli), after stanli_model(), and after sample_model() with 4
# parallel chains, with and without inits; then times four busy children
# against one, before stanli and after (about 1 if they get a core each,
# about 4 if they share one). Run on a compute node with 4 CPUs.
# Written 2026-10-12, when CmdStan ran 5x slower on one node in the 0.20.0
# fit while stanli in the same R process did not.
allowed <- function(pid) sub("^Cpus_allowed_list:\\s*", "",
  grep("^Cpus_allowed_list", readLines(sprintf("/proc/%s/status", pid)), value = TRUE))
child <- function() system("grep Cpus_allowed_list /proc/self/status | cut -f2", intern = TRUE)
show <- function(stage) cat(sprintf("%-22s R %-12s child %-12s R threads %s\n", stage, allowed(Sys.getpid()),
  child(), length(list.files(sprintf("/proc/%d/task", Sys.getpid())))))
loop <- "i=0; while [ $i -lt 2000000 ]; do i=$((i+1)); done"
busy <- function(stage) {
  t1 <- system.time(system(sprintf("sh -c '%s'", loop)))[["elapsed"]]
  t4 <- system.time(system(sprintf("for k in 1 2 3 4; do sh -c '%s' & done; wait", loop)))[["elapsed"]]
  cat(sprintf("%-22s one busy child %.1f s, four at once %.1f s: ratio %.2f\n", stage, t1, t4, t4 / t1))
}
cat("node", Sys.info()[["nodename"]], "\n")
show("start")
busy("before stanli")
suppressPackageStartupMessages(library(stanli))
show("library(stanli)")
code <- "data { int N; vector[N] y; } parameters { real mu; real<lower=0> sigma; }
model { y ~ normal(mu, sigma); }"
m <- stanli_model(code = code, data = list(N = 1000, y = rnorm(1000)))
show("stanli_model()")
busy("after stanli_model()")
u <- unconstrain(m, list(mu = 0, sigma = 1))
f <- sample_model(m, chains = 4, warmup = 200, samples = 200, init = rbind(u, u, u, u),
                  parallel_chains = 4, refresh = 0)
show("sample_model(init)")
f <- sample_model(m, chains = 4, warmup = 200, samples = 200, parallel_chains = 4, refresh = 0)
show("sample_model()")
busy("after sample_model()")
