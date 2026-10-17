# The DDM-4 routing change in context: multilevel.R's ddm4 program (10
# participants x 2 conditions x 50 trials), from the base tree and the
# current one, compared at posterior draws and tail points, then timed per
# gradient in CmdStan and stanli.
#
#   Rscript benchmarks/stanli/multilevel.R emit --model ddm4     # once, for the inputs
#   Rscript benchmarks/ddm_router/program_ab.R --base <worktree of the base> [--pkg .]
#                                              [--reps 21] [--draws 40] [--out DIR]
#
# Variants, as in grid_ab.R: `base`, `new` (or, when the trees agree, `cert`
# spliced into the base by text), `w5` (no classic branch: every trial
# through the 5-argument wiener_lpdf() at sv = 0) and `noroute` (the router
# dropped from the branch condition: unsafe in the tails, the upper bound on
# what the certificate can save). Points: DRAWS posterior draws of a short
# CmdStan fit of `base` (2 chains x 300 + 100, cached in OUT), then as many
# again moved by N(0, 1.5) and N(0, 3) on the unconstrained scale, which
# reach slow trials, extreme drifts and narrow boundaries. Timing: one block
# = one gradient at each posterior draw; blocks alternate over variants and
# engines in random order, REPS rounds, single-threaded.
source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile()
opt <- gp_args(list(base = "", pkg = ".", reps = 21L, draws = 40L,
                    out = file.path(Sys.getenv("TEMP"), "ddm_router_program")))
stopifnot(nzchar(opt$base))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
inputs <- normalizePath("benchmarks/results/stanli_ml/ddm4_inputs.rds")
d4 <- readRDS(inputs)

prog <- function(dir) {
  gp_load(dir)
  prior <- suppressMessages(cogmod_priors(d4$f, d4$df)); sv <- cogmod_stanvars(d4$f)
  list(code = as.character(brms::make_stancode(d4$f, data = d4$df, prior = prior, stanvars = sv)),
       sdat = lapply(unclass(brms::make_standata(d4$f, data = d4$df, prior = prior, stanvars = sv)), identity),
       inits = cogmod_inits(d4$f, d4$df))
}
p_base <- prog(opt$base)
p_new <- prog(opt$pkg)
suppressPackageStartupMessages(library(stanli))
cat("stanli", as.character(packageVersion("stanli")), "- CmdStan", as.character(cmdstanr::cmdstan_version()), "\n")

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
route_noroute <- paste0(
  "    if (sigmadrift == 0 && t > 0.01 * square(boundary)) {\n",
  "      return wiener_lpdf(y | boundary, tau0, w, v);\n",
  "    }\n")
subs1 <- function(code, old, new) {
  n <- lengths(regmatches(code, gregexpr(old, code, fixed = TRUE)))
  if (n != 1) return(NULL)
  sub(old, new, code, fixed = TRUE)
}
codes <- list(base = p_base$code)
if (!identical(p_new$code, p_base$code)) codes$new <- p_new$code else codes$cert <- subs1(p_base$code, route_old, route_cert)
codes$w5 <- subs1(p_base$code, route_old, "")
codes$noroute <- subs1(p_base$code, route_old, route_noroute)
codes <- Filter(Negate(is.null), codes)
cat("variants:", paste(names(codes), collapse = ", "), "\n")
stopifnot(identical(p_new$sdat, p_base$sdat))
sdat <- p_base$sdat
dj <- file.path(opt$out, "ddm4.json"); cmdstanr::write_stan_json(sdat, dj)
init1 <- p_base$inits(1)

fits <- list(); mods <- list(); sms <- list()
for (v in names(codes)) {
  sf <- file.path(opt$out, paste0("ddm4_", v, ".stan")); writeLines(codes[[v]], sf)
  cat("compiling", v, "\n")
  mods[[v]] <- gp_compile(sf)
  fc <- mods[[v]]$sample(data = dj, init = list(init1), chains = 1, iter_warmup = 1, iter_sampling = 1,
                         fixed_param = TRUE, refresh = 0, show_messages = FALSE)
  fc$init_model_methods(verbose = FALSE)
  fits[[v]] <- fc
  sms[[v]] <- stanli_model(code = codes[[v]], data = sdat, threads_per_chain = 1)
}

# ---- Points ------------------------------------------------------------------
draws_file <- file.path(opt$out, "draws_u.rds")
if (!file.exists(draws_file)) {
  cat("fitting base for posterior draws\n")
  fit <- mods$base$sample(data = dj, init = p_base$inits, chains = 2, parallel_chains = 2,
                          iter_warmup = 300, iter_sampling = 100, seed = 1, refresh = 0,
                          show_messages = FALSE)
  U <- posterior::as_draws_matrix(fits$base$unconstrain_draws(files = fit$output_files(),
                                                              format = "draws_matrix"))
  saveRDS(matrix(as.numeric(U), nrow(U)), draws_file)
}
U <- readRDS(draws_file)
rows <- round(seq(1, nrow(U), length.out = opt$draws))
D <- lapply(rows, function(i) as.numeric(U[i, ]))
set.seed(7)
tails <- c(lapply(D, function(q) q + rnorm(length(q), 0, 1.5)),
           lapply(D, function(q) q + rnorm(length(q), 0, 3)))
pts <- c(setNames(D, paste0("draw", seq_along(D))), setNames(tails, paste0("tail", seq_along(tails))))

# ---- Exactness ---------------------------------------------------------------
ev <- list(
  cmdstan = function(v, q) tryCatch({ g <- fits[[v]]$grad_log_prob(q); list(lp = attr(g, "log_prob"), grad = as.numeric(g)) },
                                    error = function(e) list(lp = NA_real_, grad = NA_real_, err = conditionMessage(e))),
  stanli = function(v, q) tryCatch(log_prob_grad(sms[[v]], q),
                                   error = function(e) list(lp = NA_real_, grad = NA_real_, err = conditionMessage(e))))
chk <- list()
for (e in names(ev)) for (nm in names(pts)) {
  b <- ev[[e]]("base", pts[[nm]])
  for (v in setdiff(names(codes), "base")) {
    x <- ev[[e]](v, pts[[nm]])
    fin <- is.finite(b$lp) && is.finite(x$lp) && all(is.finite(b$grad)) && all(is.finite(x$grad))
    chk[[length(chk) + 1]] <- data.frame(engine = e, point = nm, variant = v, lp_base = b$lp,
      lp_identical = identical(x$lp, b$lp), grad_identical = identical(x$grad, b$grad),
      lp_absdiff = abs(x$lp - b$lp),
      grad_relerr = if (fin) max(abs(x$grad - b$grad) / pmax(1, abs(b$grad))) else NA_real_,
      finite = fin)
  }
}
chk <- do.call(rbind, chk)
utils::write.csv(chk, file.path(opt$out, "program_check.csv"), row.names = FALSE)
cat("\nexactness against base (", length(pts), "points per engine):\n")
for (k in split(chk, list(chk$engine, chk$variant), drop = TRUE)) {
  cat(sprintf("  %-7s %-7s lp identical %3d/%d, grad identical %3d/%d, non-finite %d; max |dlp| %.3g, max rel dgrad %.3g\n",
              k$engine[1], k$variant[1], sum(k$lp_identical), nrow(k), sum(k$grad_identical), nrow(k),
              sum(!k$finite), max(k$lp_absdiff, na.rm = TRUE), max(k$grad_relerr, na.rm = TRUE)))
}

# ---- Cost --------------------------------------------------------------------
blk <- function(e, v) {
  if (e == "cmdstan") function(k) for (i in seq_len(k)) for (q in D) fits[[v]]$grad_log_prob(q)
  else function(k) for (i in seq_len(k)) for (q in D) log_prob_grad(sms[[v]], q)
}
cells <- expand.grid(engine = names(ev), variant = names(codes), stringsAsFactors = FALSE)
passes <- vapply(seq_len(nrow(cells)), function(i) {
  f <- blk(cells$engine[i], cells$variant[i]); f(1)
  t1 <- system.time(f(1))[["elapsed"]]
  max(1L, as.integer(ceiling(0.5 / max(t1, 0.005))))
}, integer(1))
tm <- list()
for (r in seq_len(opt$reps)) for (i in sample(nrow(cells))) {
  st <- system.time(blk(cells$engine[i], cells$variant[i])(passes[i]))[["elapsed"]]
  tm[[length(tm) + 1]] <- data.frame(rep = r, engine = cells$engine[i], variant = cells$variant[i],
                                     us = st / (passes[i] * length(D)) * 1e6)
}
tm <- do.call(rbind, tm)
utils::write.csv(tm, file.path(opt$out, "program_time.csv"), row.names = FALSE)
cat("\nus per gradient, median over", opt$reps, "rounds; ratio to base within round, median [min, max]:\n")
for (e in names(ev)) for (v in names(codes)) {
  x <- tm[tm$engine == e & tm$variant == v, ]; b <- tm[tm$engine == e & tm$variant == "base", ]
  rr <- x$us[order(x$rep)] / b$us[order(b$rep)]
  cat(sprintf("  %-7s %-7s %6.0f us  ratio %.3f [%.3f, %.3f]\n", e, v, stats::median(x$us),
              stats::median(rr), min(rr), max(rr)))
}
