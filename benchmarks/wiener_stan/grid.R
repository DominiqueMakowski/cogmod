# Function-level accuracy of the Stan-language Wiener densities against an
# independent reference, beside Stan's own implementations.
#
#   Rscript benchmarks/wiener_stan/grid.R --base <tree> [--batches 8] [--n 5000] [--out DIR]
#
# One Stan program holds the base tree's prelude and both variants of this
# directory (renamed apart), and a data switch picks the density:
#   base     cogmod_ddm_decision_lpdf() of the base tree (certificate + router
#            -> classic, else wiener5), as cogmod calls it today
#   classic  Stan's 4-parameter wiener_lpdf(), raw (sv = 0 only)
#   w5       Stan's 5-parameter wiener_lpdf(), raw
#   direct   wfpt_direct.stan
#   rec      wfpt_rec.stan
#   unrolled wfpt_unrolled.stan
# Each trial's t, v, a, w and sv are parameters of their own, so one gradient
# gives every trial's partials and a generated quantity every trial's value.
# Two modes: sv = 0 (the DDM-4 call) and sv > 0 (the DDM-5 call).
#
# Reference (`ref_ll()` below, R, independent of the package's .ddm_lfpt()):
# the small-time series with 61 terms below u = 0.05, its k = 0 and k = -1
# terms paired through expm1() as in wfpt_rec.stan; the large-time series
# with 200 terms from u = 0.05, sin() taken at the start point nearer its
# boundary. Both converge to double precision on their side of 0.05; where
# they overlap (u in [0.05, 0.4]) they are compared and the gap reported.
# Partials of the reference by central differences on the natural scale; a
# partial whose differences at h and 2h disagree by more than 1e-7
# (relative) is too noisy to judge and is left out of the comparison.
#
# Points as benchmarks/ddm_router/grid_ab.R: three batches in four
# log-uniform over wide ranges, one like a DDM-4 posterior's trials; sv
# log-uniform over [0.01, 5], or softplus(N(0.25, 0.3)) in the latter. Batch
# 7 (of every 8) puts the boundary between 1e-160 and 1e-3 or 200 and 1e80,
# where no fit should go but a tail proposal can: there only finiteness is
# judged, since a log density of -1e191 has no useful reference.
source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile()
opt <- gp_args(list(base = "", batches = 8L, n = 5000L,
                    out = file.path(Sys.getenv("TEMP"), "wiener_stan_grid")))
stopifnot(nzchar(opt$base))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

gp_load(opt$base)
prelude <- get(".DDM_STAN_PRELUDE", envir = asNamespace("cogmod"))
lfpt_pkg <- get(".ddm_lfpt", envir = asNamespace("cogmod"))
variant_code <- function(v, suf) {
  x <- paste(readLines(file.path("benchmarks", "wiener_stan", paste0("wfpt_", v, ".stan"))), collapse = "\n")
  x <- gsub("cogmod_wfpt0(", paste0("cogmod_wfpt0", suf, "("), x, fixed = TRUE)
  gsub("cogmod_wfpt_lpdf(", paste0("cogmod_wfpt", suf, "_lpdf("), x, fixed = TRUE)
}
program <- paste0("functions {\n", prelude, "\n", variant_code("direct", "_direct"), "\n",
                  variant_code("rec", "_rec"), "\n", variant_code("unrolled", "_unrolled"), "
real dispatch(int variant, real t, real v, real a, real w, real sv) {
  if (variant == 1) return cogmod_ddm_decision_lpdf(t | v, a, w, sv, 0, 0);
  if (variant == 2) return wiener_lpdf(t + 1e-10 | a, 1e-10, w, v);
  if (variant == 3) return wiener_lpdf(t + 1e-10 | a, 1e-10, w, v, sv);
  if (variant == 4) return cogmod_wfpt_direct_lpdf(t | a, w, v, sv);
  if (variant == 5) return cogmod_wfpt_rec_lpdf(t | a, w, v, sv);
  return cogmod_wfpt_unrolled_lpdf(t | a, w, v, sv);
}
}
data { int<lower=1> N; int<lower=1, upper=6> variant; int<lower=0, upper=1> use_sv; }
parameters {
  vector<lower=0>[N] t;
  vector[N] v;
  vector<lower=0>[N] a;
  vector<lower=0, upper=1>[N] w;
  vector<lower=0>[N] sv;
}
model {
  for (n in 1:N) {
    if (use_sv) target += dispatch(variant, t[n], v[n], a[n], w[n], sv[n]);
    else target += dispatch(variant, t[n], v[n], a[n], w[n], 0);
  }
}
generated quantities {
  vector[N] ll;
  for (n in 1:N) {
    if (use_sv) ll[n] = dispatch(variant, t[n], v[n], a[n], w[n], sv[n]);
    else ll[n] = dispatch(variant, t[n], v[n], a[n], w[n], 0);
  }
}
")
sf <- file.path(opt$out, "grid.stan")
writeLines(program, sf)
cat("compiling\n")
mod <- gp_compile(sf)

N <- opt$n
init0 <- list(t = rep(0.6, N), v = rep(1, N), a = rep(1.5, N), w = rep(0.5, N), sv = rep(0.5, N))
fit_for <- function(variant, use_sv) {
  dj <- file.path(opt$out, sprintf("data_%d_%d.json", variant, use_sv))
  cmdstanr::write_stan_json(list(N = N, variant = variant, use_sv = use_sv), dj)
  fit <- mod$sample(data = dj, init = list(init0), chains = 1, iter_warmup = 1, iter_sampling = 1,
                    fixed_param = TRUE, refresh = 0, show_messages = FALSE)
  fit$init_model_methods(verbose = FALSE)
  fit
}
modes <- list(sv0 = c(base = 1L, classic = 2L, w5 = 3L, direct = 4L, rec = 5L, unrolled = 6L),
              sv = c(base = 1L, w5 = 3L, direct = 4L, rec = 5L, unrolled = 6L))
fits <- lapply(names(modes), function(m) lapply(modes[[m]], fit_for, use_sv = as.integer(m == "sv")))
names(fits) <- names(modes)

# ---- Reference ---------------------------------------------------------------
ref_small <- function(u, wl, w, m = 30) {
  r <- 2 / u
  s <- -expm1(-w * r) - w * (1 + exp(-w * r))
  for (k in 1:m) s <- s + (wl + 2 * k) * exp(-k * (k + wl) * r)
  for (j in 2:m) s <- s + (wl - 2 * j) * exp(-j * (j - wl) * r)
  log(s) - wl^2 / (2 * u) - 0.5 * log(2 * pi) - 1.5 * log(u)
}
ref_large <- function(u, wl, w, K = 200) {
  flip <- wl > 0.5
  th <- pi * ifelse(flip, w, wl)
  s <- 0
  for (k in 1:K) s <- s + ifelse(flip & k %% 2 == 0, -1, 1) * k * exp(-(k^2 - 1) * pi^2 * u / 2) * sin(k * th)
  log(pi) + log(s) - pi^2 * u / 2
}
ref_lfpt0 <- function(u, wl, w) {
  out <- numeric(length(u)); sm <- u < 0.05
  if (any(sm)) out[sm] <- ref_small(u[sm], wl[sm], w[sm])
  if (any(!sm)) out[!sm] <- ref_large(u[!sm], wl[!sm], w[!sm])
  out
}
# Stan's convention: upper boundary, w from the lower one.
ref_ll <- function(t, v, a, w, sv) {
  wl <- 1 - w
  ref_lfpt0(t / a^2, wl, w) - 2 * log(a) +
    (sv^2 * a^2 * wl^2 + 2 * v * a * wl - v^2 * t) / (2 * (1 + sv^2 * t)) - 0.5 * log1p(sv^2 * t)
}
ref_grad <- function(p) {
  steps <- list(t = 1e-5 * p$t, v = 1e-5 * pmax(1, abs(p$v)), a = 1e-5 * p$a,
                w = 1e-5 * pmin(p$w, 1 - p$w), sv = 1e-5 * pmax(p$sv, 1e-3))
  f <- function(q) ref_ll(q$t, q$v, q$a, q$w, q$sv)
  g <- matrix(NA_real_, nrow(p), 5, dimnames = list(NULL, names(steps)))
  noise <- g
  for (j in names(steps)) {
    cd <- function(h) { up <- p; dn <- p; up[[j]] <- p[[j]] + h; dn[[j]] <- p[[j]] - h; (f(up) - f(dn)) / (2 * h) }
    g1 <- cd(steps[[j]]); g2 <- cd(2 * steps[[j]])
    g[, j] <- g1
    noise[, j] <- abs(g1 - g2) / pmax(1, abs(g1))
  }
  list(g = g, noise = noise)
}

# ---- Points ------------------------------------------------------------------
points <- function(b) {
  set.seed(b)
  if (b %% 8 == 7) {  # boundaries no fit should reach: only finiteness is judged
    a <- exp(c(runif(N / 2, log(1e-160), log(1e-3)), runif(N / 2, log(200), log(1e80))))
    t <- exp(runif(N, log(2e-10), log(100)))
    w <- plogis(runif(N, -12, 12))
    v <- sign(runif(N, -1, 1)) * exp(runif(N, log(1e-3), log(60)))
    sv <- exp(runif(N, log(0.01), log(5)))
    return(data.frame(t = t, v = v, a = a, w = w, sv = sv, kind = "extreme"))
  }
  if (b %% 4 == 0) {
    a <- log1p(exp(rnorm(N, 1.1, 0.4)))
    t <- runif(N, 0.05, 2.5)
    v <- rnorm(N, 1.4, 0.6) * sample(c(-1, 1), N, TRUE)
    w <- plogis(rnorm(N, 0, 0.3))
    sv <- log1p(exp(rnorm(N, 0.25, 0.3)))
  } else {
    a <- exp(runif(N, log(1e-3), log(200)))
    t <- exp(runif(N, log(1e-4), log(50))) * a^2 * exp(runif(N, -1, 1))
    w <- plogis(runif(N, -12, 12))
    v <- sign(runif(N, -1, 1)) * exp(runif(N, log(1e-3), log(60)))
    sv <- exp(runif(N, log(0.01), log(5)))
  }
  t <- pmax(t, 2e-10)
  data.frame(t = t, v = v, a = a, w = w, sv = sv, kind = if (b %% 4 == 0) "posterior" else "wide")
}

# ---- Evaluate ----------------------------------------------------------------
eval_one <- function(fit, p) {
  up <- fit$unconstrain_variables(list(t = p$t, v = p$v, a = p$a, w = p$w, sv = p$sv))
  g <- tryCatch(as.numeric(fit$grad_log_prob(up, jacobian = FALSE)), error = function(e) rep(NA_real_, 5 * N))
  ll <- tryCatch(as.numeric(fit$constrain_variables(up, generated_quantities = TRUE)$ll),
                 error = function(e) rep(NA_real_, N))
  G <- matrix(g, N)
  # natural-scale partials: t, a, sv through exp(); w through inv_logit().
  G[, 1] <- G[, 1] / p$t; G[, 3] <- G[, 3] / p$a; G[, 4] <- G[, 4] / (p$w * (1 - p$w)); G[, 5] <- G[, 5] / p$sv
  list(ll = ll, G = G)
}
rows <- list(); worst <- list(); overlap <- c(); pkg_gap <- c(); ext <- list()
for (b in seq_len(opt$batches)) {
  p <- points(b)
  for (m in names(modes)) {
    q <- p
    if (m == "sv0") q$sv <- 0
    if (p$kind[1] == "extreme") {
      keep <- if (m == "sv0") 1:4 else 1:5
      for (v in names(modes[[m]])) {
        r <- eval_one(fits[[m]][[v]], q)
        gfin <- apply(is.finite(r$G[, keep, drop = FALSE]), 1, all)
        ext[[length(ext) + 1]] <- data.frame(batch = b, mode = m, variant = v, n = N,
          value_finite_grad_nonfinite = sum(is.finite(r$ll) & !gfin),
          value_neginf_grad_nonfinite = sum(r$ll == -Inf & !gfin, na.rm = TRUE),
          value_neginf_grad_finite = sum(r$ll == -Inf & gfin, na.rm = TRUE),
          value_nan = sum(is.nan(r$ll)))
      }
      next
    }
    ref <- ref_ll(q$t, q$v, q$a, q$w, q$sv)
    rg <- ref_grad(q)
    if (any(!is.finite(ref))) cat(sprintf("  %s: reference not finite at %d points\n", m, sum(!is.finite(ref))))
    keep <- if (m == "sv0") 1:4 else 1:5
    u <- q$t / q$a^2
    ov <- u >= 0.05 & u <= 0.4
    if (any(ov)) overlap <- c(overlap, max(abs(ref_small(u[ov], 1 - q$w[ov], q$w[ov]) -
                                                 ref_large(u[ov], 1 - q$w[ov], q$w[ov])) / pmax(1, abs(ref[ov]))))
    pk <- lfpt_pkg(q$t, -q$v, q$a, 1 - q$w, q$sv)
    fin <- is.finite(pk) & is.finite(ref)
    pkg_gap <- c(pkg_gap, max(abs(pk - ref)[fin] / pmax(1, abs(ref[fin]))))
    region <- cut(ref, c(-Inf, -1000, -30, Inf), labels = c("far", "tail", "body"))
    for (v in names(modes[[m]])) {
      r <- eval_one(fits[[m]][[v]], q)
      dll <- abs(r$ll - ref) / pmax(1, abs(ref))
      gerr <- abs(r$G[, keep, drop = FALSE] - rg$g[, keep, drop = FALSE]) / pmax(1, abs(rg$g[, keep, drop = FALSE]))
      gerr[rg$noise[, keep, drop = FALSE] > 1e-7] <- NA
      gmax <- apply(gerr, 1, function(x) if (all(is.na(x))) NA_real_ else max(x, na.rm = TRUE))
      nonfin_ll <- !is.finite(r$ll)
      nonfin_g <- !apply(is.finite(r$G[, keep, drop = FALSE]), 1, all)
      for (rg_name in levels(region)) {
        i <- which(region == rg_name)
        if (!length(i)) next
        rows[[length(rows) + 1]] <- data.frame(
          batch = b, kind = p$kind[1], mode = m, variant = v, region = rg_name, n = length(i),
          nonfinite_ll = sum(nonfin_ll[i]), nonfinite_grad = sum(nonfin_g[i] & !nonfin_ll[i]),
          ll_max = suppressWarnings(max(dll[i][is.finite(dll[i])], -Inf)),
          ll_p99 = if (any(is.finite(dll[i]))) stats::quantile(dll[i][is.finite(dll[i])], 0.99, names = FALSE) else NA,
          grad_max = suppressWarnings(max(gmax[i][is.finite(gmax[i])], -Inf)),
          grad_p99 = if (any(is.finite(gmax[i]))) stats::quantile(gmax[i][is.finite(gmax[i])], 0.99, names = FALSE) else NA,
          grad_judged = sum(is.finite(gmax[i])))
      }
      o <- order(-ifelse(is.finite(dll), dll, Inf))[1:3]
      worst[[length(worst) + 1]] <- cbind(mode = m, variant = v, q[o, 1:5], u = u[o], ref = ref[o],
                                          ll = r$ll[o], rel_dll = dll[o], rel_dgrad = gmax[o])
    }
  }
  cat(sprintf("batch %d (%s) done\n", b, p$kind[1]))
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(opt$out, "grid.csv"), row.names = FALSE)
utils::write.csv(do.call(rbind, worst), file.path(opt$out, "grid_worst.csv"), row.names = FALSE)
cat(sprintf("\nreference: small- vs large-time series on u in [0.05, 0.4], max rel gap %.3g\n", max(overlap)))
cat(sprintf("package .ddm_lfpt() vs reference, max rel gap %.3g\n", max(pkg_gap)))
agg <- do.call(rbind, lapply(split(res, list(res$mode, res$variant, res$region), drop = TRUE), function(x)
  data.frame(mode = x$mode[1], variant = x$variant[1], region = x$region[1], n = sum(x$n),
             nonfinite_ll = sum(x$nonfinite_ll), nonfinite_grad = sum(x$nonfinite_grad),
             ll_max = max(x$ll_max), ll_p99 = max(x$ll_p99, na.rm = TRUE),
             grad_max = max(x$grad_max), grad_p99 = max(x$grad_p99, na.rm = TRUE),
             grad_judged = sum(x$grad_judged))))
agg <- agg[order(agg$mode, agg$region, agg$variant), ]
utils::write.csv(agg, file.path(opt$out, "grid_summary.csv"), row.names = FALSE)
print(agg, row.names = FALSE, digits = 3)
if (length(ext)) {
  ext <- do.call(rbind, ext)
  utils::write.csv(ext, file.path(opt$out, "grid_extreme.csv"), row.names = FALSE)
  cat("\nboundaries from 1e-160 to 1e80 (finiteness only):\n")
  print(ext, row.names = FALSE)
}
