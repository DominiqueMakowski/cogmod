# Does a change to the DDM's classic-path routing change any trial's value
# or partials? Function-level A/B of cogmod_ddm_decision_lpdf() at
# sigmadrift = sigmabias = sigmandt = 0 (the DDM-4 call site), base tree
# against the current one, through CmdStan's model methods.
#
#   Rscript benchmarks/ddm_router/grid_ab.R --base <worktree of the base> [--pkg .]
#                                           [--batches 40] [--n 5000] [--out DIR]
#
# One Stan program per variant holds the prelude and nothing else: t, v, a
# and w are parameter vectors with one entry per trial, so one gradient gives
# every trial's own partials (no trial shares a parameter), and a generated
# quantity gives every trial's own value. The variants:
#   base     the base tree's prelude
#   new      the current tree's (left out if its text is the base's)
#   cert     the base with the certificate spliced in by text (README.md),
#            left out when `new` exists - it is then what `new` should be
#   w5       the base with the classic branch removed: every trial through
#            the 5-argument wiener_lpdf() at sv = 0
# `new` / `cert` are expected to equal `base` bit for bit; `w5` is reported.
#
# Points: log-uniform over wide ranges (as certificate_check.R), so that the
# certificate's edges, the router's -600 and the small-time cut are all hit,
# plus one batch in four like a DDM-4 posterior's trials. Each point is
# classified with an R transcription of the router and the certificate.
source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()
opt <- gp_args(list(base = "", pkg = ".", batches = 40L, n = 5000L,
                    out = file.path(Sys.getenv("TEMP"), "ddm_router_grid")))
stopifnot(nzchar(opt$base))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

prelude <- function(dir) {
  gp_load(dir)
  get(".DDM_STAN_PRELUDE", envir = asNamespace("cogmod"))
}
p_base <- prelude(opt$base)
p_new <- prelude(opt$pkg)

route_old <- paste0(
  "    if (sigmadrift == 0 && t > 0.01 * square(boundary)\n",
  "        && cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600) {\n",
  "      return wiener_lpdf(y | boundary, tau0, w, v);\n",
  "    }\n")
route_cert <- paste0(
  "    if (sigmadrift == 0) {\n",
  "      real a2 = square(boundary);\n",
  "      if (t > 0.01 * a2\n",
  "          && ((boundary < 100 && w > 1e-3 && w < 1 - 1e-3 && t < 50 * a2\n",
  "               && v * (v * t + 2 * boundary * w) < 600)\n",
  "              || cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600)) {\n",
  "        return wiener_lpdf(y | boundary, tau0, w, v);\n",
  "      }\n",
  "    }\n")
subs1 <- function(code, old, new) {
  n <- lengths(regmatches(code, gregexpr(old, code, fixed = TRUE)))
  if (n != 1) return(NULL)
  sub(old, new, code, fixed = TRUE)
}
variants <- list(base = p_base)
if (!identical(p_new, p_base)) variants$new <- p_new else variants$cert <- subs1(p_base, route_old, route_cert)
variants$w5 <- subs1(p_base, route_old, "")
variants <- Filter(Negate(is.null), variants)
cat("variants:", paste(names(variants), collapse = ", "), "\n")

program <- function(fns) paste0("functions {\n", fns, "\n}\n", "
data { int<lower=1> N; }
parameters {
  vector<lower=0>[N] t;
  vector[N] v;
  vector<lower=0>[N] a;
  vector<lower=0, upper=1>[N] w;
}
model {
  for (n in 1:N) target += cogmod_ddm_decision_lpdf(t[n] | v[n], a[n], w[n], 0, 0, 0);
}
generated quantities {
  vector[N] ll;
  for (n in 1:N) ll[n] = cogmod_ddm_decision_lpdf(t[n] | v[n], a[n], w[n], 0, 0, 0);
}
")

N <- opt$n
dj <- file.path(opt$out, "data.json")
cmdstanr::write_stan_json(list(N = N), dj)
init0 <- list(t = rep(0.6, N), v = rep(1, N), a = rep(1.5, N), w = rep(0.5, N))
fits <- lapply(names(variants), function(v) {
  sf <- file.path(opt$out, paste0(v, ".stan"))
  writeLines(program(variants[[v]]), sf)
  cat("compiling", v, "\n")
  mod <- gp_compile(sf)
  fit <- mod$sample(data = dj, init = list(init0), chains = 1, iter_warmup = 1, iter_sampling = 1,
                    fixed_param = TRUE, refresh = 0, show_messages = FALSE)
  fit$init_model_methods(verbose = FALSE)
  fit
})
names(fits) <- names(variants)

# ---- Points ------------------------------------------------------------------
router <- function(t, v, a, w) {  # cogmod_ddm_log_density_scale() at sv = sw = st0 = 0
  w_lo <- pmax(w, 1e-12); w_hi <- pmin(w, 1 - 1e-12)
  t_lo <- pmax(t, 1e-12); tau_hi <- t / a^2; tau_lo <- t_lo / a^2
  ws <- pmin(pmax(sqrt(tau_hi), w_lo), w_hi)
  small <- log(ws) - ws^2 / (2 * tau_hi) - 0.5 * log(2 * pi * tau_hi^3)
  wl <- pmin(pmax(0.5, w_lo), w_hi)
  large <- log(pi) - pi^2 * tau_lo / 2 + log(sin(pi * wl))
  num <- -v^2 * t_lo - 2 * v * a * ifelse(v > 0, w_lo, w_hi)
  -2 * log(a) + pmin(small, large) + num / 2
}
cert <- function(t, v, a, w) {
  a2 <- a^2
  t > 0.01 * a2 & a < 100 & w > 1e-3 & w < 1 - 1e-3 & t < 50 * a2 & v * (v * t + 2 * a * w) < 600
}
points <- function(b) {
  set.seed(b)
  if (b %% 4 == 0) {  # like multilevel.R's ddm4 truth
    a <- log1p(exp(rnorm(N, 1.1, 0.4)))
    t <- runif(N, 0.05, 2.5)
    v <- rnorm(N, 1.4, 0.6) * sample(c(-1, 1), N, TRUE)
    w <- plogis(rnorm(N, 0, 0.3))
  } else {
    a <- exp(runif(N, log(1e-3), log(200)))
    t <- exp(runif(N, log(1e-4), log(50))) * a^2 * exp(runif(N, -1, 1))
    w <- plogis(runif(N, -12, 12))
    v <- sign(runif(N, -1, 1)) * exp(runif(N, log(1e-3), log(60)))
  }
  t <- pmax(t, 2e-10)  # above the 1e-10 offset, where the density is defined
  data.frame(t = t, v = v, a = a, w = w)
}
classify <- function(p) {
  small <- !(p$t > 0.01 * p$a^2)
  ce <- !small & cert(p$t, p$v, p$a, p$w)
  r <- router(p$t, p$v, p$a, p$w)
  ifelse(small, "small_time", ifelse(ce, "certified", ifelse(r > -600, "router_classic", "router_w5")))
}

# ---- Evaluate ----------------------------------------------------------------
eval_one <- function(fit, up) {
  g <- tryCatch(fit$grad_log_prob(up, jacobian = TRUE), error = function(e) conditionMessage(e))
  if (is.character(g)) return(list(err = g))
  gq <- tryCatch(fit$constrain_variables(up, generated_quantities = TRUE)$ll,
                 error = function(e) conditionMessage(e))
  list(lp = attr(g, "log_prob"), grad = as.numeric(g), ll = gq)
}
rows <- list(); cls_all <- character(0)
for (b in seq_len(opt$batches)) {
  p <- points(b)
  cls <- classify(p); cls_all <- c(cls_all, cls)
  up <- fits$base$unconstrain_variables(list(t = p$t, v = p$v, a = p$a, w = p$w))
  r <- lapply(fits, eval_one, up = up)
  for (v in setdiff(names(r), "base")) {
    x <- r[[v]]; y <- r$base
    if (!is.null(x$err) || !is.null(y$err)) {
      rows[[length(rows) + 1]] <- data.frame(batch = b, variant = v, lp_identical = NA,
        grad_identical = NA, ll_identical = NA, ll_maxabs = NA, grad_maxrel = NA,
        err = paste(c(x$err, y$err), collapse = " | "))
      next
    }
    # Per trial: the gradient holds (t, v, a, w) blocks of N each, so trial
    # n's partials are entries n, N + n, 2N + n, 3N + n.
    G <- matrix(x$grad, N); H <- matrix(y$grad, N)
    fin <- is.finite(x$ll) & is.finite(y$ll)
    rows[[length(rows) + 1]] <- data.frame(batch = b, variant = v,
      lp_identical = identical(x$lp, y$lp), grad_identical = identical(x$grad, y$grad),
      ll_identical = identical(x$ll, y$ll),
      ll_maxabs = if (any(fin)) max(abs(x$ll - y$ll)[fin]) else NA,
      grad_maxrel = max((abs(G - H) / pmax(1, abs(H)))[fin, ], na.rm = TRUE), err = "")
    if (v %in% c("new", "cert") && !identical(x$ll, y$ll)) {
      bad <- which(!(x$ll == y$ll | (is.na(x$ll) & is.na(y$ll))))
      print(cbind(p[bad[seq_len(min(5, length(bad)))], ], class = cls[bad[seq_len(min(5, length(bad)))]]))
    }
  }
  last <- do.call(rbind, utils::tail(rows, length(fits) - 1))
  cat(sprintf("batch %2d: %s\n", b, paste(sprintf("%s lp/grad/values identical %s/%s/%s", last$variant,
    last$lp_identical, last$grad_identical, last$ll_identical), collapse = "; ")))
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(opt$out, "grid_ab.csv"), row.names = FALSE)
cat("\npoints per class:\n"); print(table(cls_all))
for (v in unique(res$variant)) {
  x <- res[res$variant == v, ]
  cat(sprintf("%-5s %d batches (%d errored): lp identical %d, grad identical %d, values identical %d; max |dll| %.3g, max rel dgrad %.3g\n",
              v, nrow(x), sum(nzchar(x$err)), sum(x$lp_identical, na.rm = TRUE),
              sum(x$grad_identical, na.rm = TRUE), sum(x$ll_identical, na.rm = TRUE),
              max(x$ll_maxabs, na.rm = TRUE), max(x$grad_maxrel, na.rm = TRUE)))
}
