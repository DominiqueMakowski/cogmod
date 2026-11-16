# The Stan-language Wiener densities in context: multilevel.R's ddm4 or ddm5
# program (10 participants x 2 conditions x 50 trials) from the base tree and
# from one patched tree per variant (apply_variant.R), compared at posterior
# draws and tail points, then timed per gradient in CmdStan and stanli.
#
#   Rscript benchmarks/stanli/multilevel.R emit --model ddm5   # once, for the inputs
#   Rscript benchmarks/wiener_stan/program_ab.R --base <tree> --protos direct=<tree>,rec=<tree>
#                                               --model ddm4 [--reps 21] [--draws 40] [--out DIR]
#
# As benchmarks/ddm_router/program_ab.R, except that the variants are not
# expected to be bit-identical to the base: the report is the largest log
# density and gradient gap, and the count of points where either side is not
# finite. Points: DRAWS posterior draws of a short CmdStan fit of the base
# (2 chains x 300 + 100, cached in OUT), then as many again moved by
# N(0, 1.5) and N(0, 3) on the unconstrained scale. Timing: one block = one
# gradient at each posterior draw; blocks alternate over variants and engines
# in random order, REPS rounds, single-threaded.
source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile()
opt <- gp_args(list(base = "", protos = "", model = "ddm4", reps = 21L, draws = 40L,
                    out = file.path(Sys.getenv("TEMP"), "wiener_stan_program")))
stopifnot(nzchar(opt$base), nzchar(opt$protos), opt$model %in% c("ddm4", "ddm5"))
out <- file.path(opt$out, opt$model)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
d <- readRDS(normalizePath(sprintf("benchmarks/results/stanli_ml/%s_inputs.rds", opt$model)))

prog <- function(dir) {
  gp_load(dir)
  prior <- suppressMessages(cogmod_priors(d$f, d$df)); sv <- cogmod_stanvars(d$f)
  list(code = as.character(brms::make_stancode(d$f, data = d$df, prior = prior, stanvars = sv)),
       sdat = lapply(unclass(brms::make_standata(d$f, data = d$df, prior = prior, stanvars = sv)), identity),
       inits = cogmod_inits(d$f, d$df))
}
pr <- strsplit(strsplit(opt$protos, ",")[[1]], "=")
trees <- c(base = opt$base, stats::setNames(vapply(pr, `[`, "", 2), vapply(pr, `[`, "", 1)))
progs <- lapply(trees, prog)
for (v in names(progs)[-1]) stopifnot(identical(progs[[v]]$sdat, progs$base$sdat))
suppressPackageStartupMessages(library(stanli))
cat("stanli", as.character(packageVersion("stanli")), "- CmdStan", as.character(cmdstanr::cmdstan_version()),
    "- model", opt$model, "- variants:", paste(names(trees), collapse = ", "), "\n")
sdat <- progs$base$sdat
dj <- file.path(out, "data.json"); cmdstanr::write_stan_json(sdat, dj)
init1 <- progs$base$inits(1)

fits <- list(); mods <- list(); sms <- list()
for (v in names(progs)) {
  sf <- file.path(out, paste0(opt$model, "_", v, ".stan")); writeLines(progs[[v]]$code, sf)
  cat("compiling", v, "\n")
  mods[[v]] <- gp_compile(sf)
  fc <- mods[[v]]$sample(data = dj, init = list(init1), chains = 1, iter_warmup = 1, iter_sampling = 1,
                         fixed_param = TRUE, refresh = 0, show_messages = FALSE)
  fc$init_model_methods(verbose = FALSE)
  fits[[v]] <- fc
  sms[[v]] <- stanli_model(code = progs[[v]]$code, data = sdat, threads_per_chain = 1)
}

# ---- Points ------------------------------------------------------------------
draws_file <- file.path(out, "draws_u.rds")
if (!file.exists(draws_file)) {
  cat("fitting base for posterior draws\n")
  fit <- mods$base$sample(data = dj, init = progs$base$inits, chains = 2, parallel_chains = 2,
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

# ---- Agreement ---------------------------------------------------------------
ev <- list(
  cmdstan = function(v, q) tryCatch({ g <- fits[[v]]$grad_log_prob(q); list(lp = attr(g, "log_prob"), grad = as.numeric(g)) },
                                    error = function(e) list(lp = NA_real_, grad = NA_real_)),
  stanli = function(v, q) tryCatch(log_prob_grad(sms[[v]], q),
                                   error = function(e) list(lp = NA_real_, grad = NA_real_)))
fin <- function(x) is.finite(x$lp) && all(is.finite(x$grad))
chk <- list()
for (e in names(ev)) for (nm in names(pts)) {
  b <- ev[[e]]("base", pts[[nm]])
  for (v in names(progs)) {
    x <- if (v == "base") b else ev[[e]](v, pts[[nm]])
    c0 <- ev$cmdstan(v, pts[[nm]])  # the same variant in CmdStan, for the engine gap
    both <- fin(b) && fin(x)
    chk[[length(chk) + 1]] <- data.frame(engine = e, point = nm, variant = v, lp_base = b$lp, lp = x$lp,
      finite_base = fin(b), finite = fin(x),
      lp_absdiff = if (both) abs(x$lp - b$lp) else NA_real_,
      grad_relerr = if (both) max(abs(x$grad - b$grad) / pmax(1, abs(b$grad))) else NA_real_,
      engine_lp_absdiff = if (fin(x) && fin(c0)) abs(x$lp - c0$lp) else NA_real_)
  }
}
chk <- do.call(rbind, chk)
utils::write.csv(chk, file.path(out, "program_check.csv"), row.names = FALSE)
cat("\nagainst base (", length(pts), "points per engine; draws / tails):\n")
for (k in split(chk, list(chk$engine, chk$variant), drop = TRUE)) {
  for (set in c("draw", "tail")) {
    s <- k[startsWith(k$point, set), ]
    cat(sprintf("  %-7s %-7s %-5s non-finite %2d (base %2d); max |dlp| %.3g, max rel dgrad %.3g; vs CmdStan max |dlp| %.3g\n",
                s$engine[1], s$variant[1], set, sum(!s$finite), sum(!s$finite_base),
                suppressWarnings(max(s$lp_absdiff, na.rm = TRUE)), suppressWarnings(max(s$grad_relerr, na.rm = TRUE)),
                suppressWarnings(max(s$engine_lp_absdiff, na.rm = TRUE))))
  }
}

# ---- Cost --------------------------------------------------------------------
blk <- function(e, v) {
  if (e == "cmdstan") function(k) for (i in seq_len(k)) for (q in D) fits[[v]]$grad_log_prob(q)
  else function(k) for (i in seq_len(k)) for (q in D) log_prob_grad(sms[[v]], q)
}
cells <- expand.grid(engine = names(ev), variant = names(progs), stringsAsFactors = FALSE)
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
utils::write.csv(tm, file.path(out, "program_time.csv"), row.names = FALSE)
cat("\nus per gradient, median over", opt$reps, "rounds; ratio to the same engine's base within round, median [min, max]:\n")
for (e in names(ev)) for (v in names(progs)) {
  x <- tm[tm$engine == e & tm$variant == v, ]; b <- tm[tm$engine == e & tm$variant == "base", ]
  rr <- x$us[order(x$rep)] / b$us[order(b$rep)]
  cat(sprintf("  %-7s %-7s %6.0f us  ratio %.3f [%.3f, %.3f]\n", e, v, stats::median(x$us),
              stats::median(rr), min(rr), max(rr)))
}
cb <- tm[tm$engine == "cmdstan" & tm$variant == "base", ]
for (v in names(progs)) {
  x <- tm[tm$engine == "stanli" & tm$variant == v, ]
  cat(sprintf("  stanli %-7s / CmdStan base: %.3f\n", v, stats::median(x$us[order(x$rep)] / cb$us[order(cb$rep)])))
}
