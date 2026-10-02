# What keeps the select LNR fast in stanli: the table in README.md, "The
# select pattern". Every variant has the same log density; only the way the
# branches, the mixture and the checks are written changes.
#
#   Rscript benchmarks/stanli/bisect.R [--out DIR]
#
# Run after `bench.R emit` (reads lnr_orig.stan, the data and the init from
# DIR, default benchmarks/results/stanli). stanli only, no compiling; about a
# minute. Writes results/bisect.csv.

source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(out = "benchmarks/results/stanli"))
suppressPackageStartupMessages(library(stanli))
inp <- readRDS(file.path(opt$out, "inputs.rds"))
x <- readLines(file.path(opt$out, "lnr_orig.stan"))
a <- grep("^functions [{]", x)
b <- grep("^data [{]", x)

# sel_log_Phi() as lnr_select_functions.stan has it.
sf <- readLines("benchmarks/stanli/lnr_select_functions.stan")
s0 <- grep("^real sel_log_Phi", sf)
selphi <- sf[s0:(s0 + grep("^}", sf[(s0 + 1):length(sf)])[1])]

lpdf <- function(ret, phi = "naive", pre = character(), tclamp = "1e-300") c(
  "functions {",
  "real naive_log_Phi(real x) { return log(0.5 * erfc(-x * 0.7071067811865476)); }",
  selphi,
  "real cogmod_lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real ndt, real poutlier, int dec) {",
  pre,
  "  real lp_out = 0.6904993792294275 - 12.5 * square(Y);",
  "  real nu_w = dec == 0 ? mu : nuone;", "  real s_w = dec == 0 ? sigmazero : sigmaone;",
  "  real nu_l = dec == 0 ? nuone : mu;", "  real s_l = dec == 0 ? sigmaone : sigmazero;",
  sprintf("  real t = fmax(Y - ndt, %s);", tclamp),
  sprintf("  real lp_dec = lognormal_lpdf(t | -nu_w, s_w) + %s_log_Phi((-nu_l - log(t)) / s_l);", phi),
  "  real w = step(Y - ndt);",
  paste0("  return ", ret, ";"), "}", "}")

MIX <- "log_mix(poutlier, lp_out, lp_dec)"
MASK <- "log_mix(poutlier, lp_out, lp_dec + (w - 1) * 1e300)"
BLEND <- paste0("w * ", MIX, " + (1 - w) * (log(poutlier) + lp_out)")
CHECKS <- c("  if (dec < 0 || dec > 1) return negative_infinity();",
            "  if (Y <= 0) return negative_infinity();")
variants <- list(
  bf = lpdf(MIX),
  sel_phi = lpdf(MIX, phi = "sel"),
  sel = lpdf(MASK, phi = "sel"),
  sel_checks = lpdf(MASK, phi = "sel", pre = CHECKS),
  bf_lse = lpdf("log_sum_exp(log(poutlier) + lp_out, log1m(poutlier) + lp_dec)"),
  bf_blend = lpdf(BLEND),
  sel_phi_blend = lpdf(BLEND, phi = "sel")
)
build <- function(funs) {
  y <- c(x[seq_len(a - 1)], funs, x[b:length(x)])
  y <- sub("sigmaone[n], sigmabias, ndt[n]", "sigmaone[n], ndt[n]", y, fixed = TRUE)
  stanli_model(code = paste(y, collapse = "\n"), data = inp$sdat)
}

ms <- lapply(variants, build)
q0 <- unconstrain(ms$bf, inp$init)
K <- 100
set.seed(3)
P <- lapply(seq_len(K), function(k) q0 + stats::rnorm(length(q0), 0, 0.1))
for (m in ms) for (k in 1:10) log_prob_grad(m, P[[k]])  # warm up
rows <- list()
for (r in 1:11) for (nm in sample(names(ms))) {
  st <- system.time(for (k in seq_len(K)) log_prob_grad(ms[[nm]], P[[k]]))
  rows[[length(rows) + 1]] <- data.frame(rep = r, variant = nm, us_per_grad = st[["elapsed"]] / K * 1e6)
}
tm <- do.call(rbind, rows)
res <- data.frame(variant = names(ms),
                  lp_init = vapply(ms, function(m) log_prob_grad(m, q0)$lp, numeric(1)),
                  us_per_grad = as.numeric(tapply(tm$us_per_grad, tm$variant, stats::median)[names(ms)]))
print(res, digits = 10, row.names = FALSE)
utils::write.csv(res, "benchmarks/stanli/results/bisect.csv", row.names = FALSE)
