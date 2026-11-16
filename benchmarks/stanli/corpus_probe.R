# Why cm_invgaussian and cm_ddm are the two cogmod rows below 1x in stanli's benchmark
# =====================================================================================
#
# stanli 0.23.2's docs/benchmarks.md (run of 2026-10-21, Apple M3 Ultra) has
# every cogmod family faster per gradient than CmdStan except cm_ddm (0.975x)
# and cm_invgaussian (0.756x; CmdStan time / stanli time). Both are fallback
# rows: their reference fit timed out, so both engines were timed at stanli's
# ten synthetic points. stanli's corpus (tests/cogmod, cogmod 732be30)
# estimates every parameter. For the Wald that sends each trial through the
# 64-point drift quadrature (`cogmod_wald_sv_lquad()`), a configuration none
# of the articles fit: they all have `sigmandt = 0`.
#
# Run from the package root:
#
#   Rscript benchmarks/stanli/corpus_probe.R emit
#   Rscript benchmarks/stanli/corpus_probe.R grad   --family invgaussian [--variants a,b] [--reps 5] [--tag NAME]
#   Rscript benchmarks/stanli/corpus_probe.R replay [--variant sv_st0] [--set post|fixed] [--tag NAME]
#   Rscript benchmarks/stanli/corpus_probe.R fit    [--variant sv_st0] [--warmup 5] [--tag NAME]
#   Rscript benchmarks/stanli/corpus_probe.R ab     [--a sv_st0] [--b q3_while] [--reps 11]
#   Rscript benchmarks/stanli/corpus_probe.R diag
#
#   emit    fetch cm_invgaussian and cm_ddm (.stan, .json) from stanli's
#           tests/cogmod at the commit that added them (the files their
#           benchmark timed: same SHA-256 as in its points file), and write
#           cogmod's own programs on the same data, one per configuration,
#           plus three rewrites of the quadrature
#   grad    CmdStan (model methods) against stanli, per gradient, at two point
#           sets: "fixed", stanli's own evaluation_point()
#           (tools/benchmark_timer.hpp), and "post", 40 points near the
#           data-generating values. REPS alternating blocks; each block cycles
#           through the set, so every gradient is at a new point, as in a fit
#   replay  stanli only, at ten of the post points or at the fixed ones: one
#           point repeated, two points alternated, one point and a 1e-3 step
#           from it alternated, the ten cycled, and the ten each repeated 3
#           and 20 times in a row (the shape of stanli's batches)
#   fit     stanli only: a short NUTS run, wall time per leapfrog step
#   ab      CmdStan only: two variants in alternating blocks
#   diag    which form stanli's region map picks for each program, and why
#
# Every mode runs under the caller's STANLI_STRUCTURED_LOOPS and records it:
# set it to 0 for the lowering without retained loops.
#
# Writes benchmarks/stanli/results/corpus_{grad,check,replay,fit,ab}[_<tag>].csv;
# the programs go to benchmarks/results/stanli/corpus/. Numbers in the
# README, section "0.23.2".

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile(), gp_methods()

mode <- commandArgs(trailingOnly = TRUE)[1]
steps <- c("emit", "grad", "replay", "fit", "ab", "diag")
if (is.na(mode) || !mode %in% steps) {
  stop("usage: Rscript benchmarks/stanli/corpus_probe.R ", paste(steps, collapse = "|"),
       " [--options]", call. = FALSE)
}
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})
opt <- gp_args(list(out = "benchmarks/results/stanli/corpus", family = "invgaussian",
                    variants = "", variant = "sv_st0", set = "post", reps = 5L, warmup = 5L,
                    a = "sv_st0", b = "q3_while", tag = ""))
res_dir <- "benchmarks/stanli/results"
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
loops <- Sys.getenv("STANLI_STRUCTURED_LOOPS")
suffix <- if (nzchar(opt$tag)) paste0("_", opt$tag) else ""
csv <- function(x, what) {
  f <- file.path(res_dir, paste0("corpus_", what, suffix, ".csv"))
  utils::write.csv(x, f, row.names = FALSE)
  cat("wrote", f, "\n")
}

# ---- Programs ----------------------------------------------------------------

# The commit that added tests/cogmod; the four files are unchanged on main
# as of v0.23.2.
STANLI_COMMIT <- "80da58c9cc21ccbe0a3612eed4e8f8521072a4f9"
corpus_url <- function(f) {
  paste0("https://raw.githubusercontent.com/seantalts/stanli/", STANLI_COMMIT, "/tests/cogmod/", f)
}

# The data-generating values of stanli's generator (tools/gen_cogmod_models.R),
# on the unconstrained scale of the corpus programs. Variants with a dpar
# fixed lack its entry, and the points are built by name.
truth <- list(
  invgaussian = list(b = 0, Intercept = log(expm1(3)), Intercept_boundary = log(expm1(0.5)),
                     Intercept_sigmadrift = log(expm1(1)), Intercept_sigmandt = log(0.1),
                     Intercept_ndt = log(0.2), Intercept_poutlier = qlogis(0.02)),
  ddm = list(b = 0, Intercept = 0.5, Intercept_boundary = log(expm1(1)), Intercept_bias = 0,
             Intercept_sigmadrift = log(expm1(0.2)), Intercept_sigmabias = qlogis(0.05),
             Intercept_sigmandt = log(0.03), Intercept_ndt = log(0.2),
             Intercept_poutlier = qlogis(0.02)))

variant_names <- list(
  invgaussian = c("corpus", "sv_st0", "sv", "st0", "plain", "q1_meancrit", "q3_while", "q4_runlse"),
  ddm = c("corpus", "ddm7", "ddm5", "ddm4"))

subs1 <- function(code, old, new) {
  n <- lengths(regmatches(code, gregexpr(old, code, fixed = TRUE)))
  if (n != 1) stop("pattern found ", n, " times: ", substr(old, 1, 60), call. = FALSE)
  sub(old, new, code, fixed = TRUE)
}

emit <- function() {
  gp_load(".")
  for (fam in names(variant_names)) {
    dir <- file.path(opt$out, fam)
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)
    for (ext in c("stan", "json")) {
      utils::download.file(corpus_url(paste0("cm_", fam, ".", ext)),
                           file.path(dir, paste0("corpus.", ext)), quiet = TRUE, mode = "wb")
    }
    cj <- jsonlite::read_json(file.path(dir, "corpus.json"), simplifyVector = TRUE)
    df <- data.frame(RT = cj$Y, x = cj$X[, 2])
    if (fam == "ddm") df$Error <- cj$dec
    write_prog <- function(name, f) {
      prior <- suppressMessages(cogmod_priors(f, df))
      sv <- cogmod_stanvars(f)
      writeLines(as.character(brms::make_stancode(f, data = df, prior = prior, stanvars = sv)),
                 file.path(dir, paste0(name, ".stan")))
      sdat <- lapply(unclass(brms::make_standata(f, data = df, prior = prior, stanvars = sv)), identity)
      cmdstanr::write_stan_json(sdat, file.path(dir, paste0(name, ".json")))
    }
    if (fam == "invgaussian") {
      F <- function(...) brms::bf(RT ~ x, boundary ~ 1, ndt ~ 1, poutlier ~ 1, ...,
                                  family = cogmod_invgaussian())
      write_prog("sv_st0", F(sigmadrift ~ 1, sigmandt ~ 1))  # the corpus formula
      write_prog("sv", F(sigmadrift ~ 1, sigmandt = 0))      # the articles' formula
      write_prog("st0", F(sigmadrift = 0, sigmandt ~ 1))
      write_prog("plain", F(sigmadrift = 0, sigmandt = 0))
      emit_rewrites(dir)
    } else {
      F <- function(...) brms::bf(RT | dec(Error) ~ x, boundary ~ 1, bias ~ 1, ndt ~ 1,
                                  poutlier ~ 1, ..., family = cogmod_ddm())
      write_prog("ddm7", F(sigmadrift ~ 1, sigmabias ~ 1, sigmandt ~ 1))  # the corpus formula
      write_prog("ddm5", F(sigmadrift ~ 1, sigmabias = 0, sigmandt = 0))
      write_prog("ddm4", F(sigmadrift = 0, sigmabias = 0, sigmandt = 0))
    }
    cat("emitted", fam, "to", dir, "\n")
  }
}

# Three rewrites of the quadrature, applied to sv_st0:
#   q1_meancrit  choose the CDF pair or the survival pair by the Wald mean
#                (m < a / v) instead of by the midpoint's CDF: two fewer
#                cogmod_log_Phi() per node
#   q3_while     the node loop as a while loop, so that it is not unrolled
#   q4_runlse    a running log_sum_exp() in place of the vector of terms
emit_rewrites <- function(dir) {
  code <- paste(readLines(file.path(dir, "sv_st0.stan")), collapse = "\n")
  file.copy(file.path(dir, "sv_st0.json"), file.path(dir, paste0(c("q1_meancrit", "q3_while",
            "q4_runlse"), ".json")), overwrite = TRUE)
  put <- function(name, x) writeLines(x, file.path(dir, paste0(name, ".stan")))
  put("q1_meancrit", subs1(code,
    "if (cogmod_wald_logcdf(t - 0.5 * st0, v, a) < -0.69314718055994529) {",
    "if ((t - 0.5 * st0) * v < a) {"))
  body <- paste0("    real v = mid + half * gx[j];\n",
                 "    real lw = lgw[j] + log(half) + normal_lpdf(v | mu, sigmadrift) - lnorm;\n")
  loop <- paste0("  for (j in 1:64) {\n", body,
                 "    terms[j] = lw + cogmod_wald_st0(t, v, boundary, st0, what);\n  }")
  put("q3_while", subs1(code, loop, paste0("  int j = 1;\n  while (j <= 64) {\n", body,
    "    terms[j] = lw + cogmod_wald_st0(t, v, boundary, st0, what);\n    j += 1;\n  }")))
  put("q4_runlse", subs1(code, paste0(loop, "\n  return log_sum_exp(terms);"),
    paste0("  real acc = negative_infinity();\n  for (j in 1:64) {\n", body,
           "    acc = log_sum_exp(acc, lw + cogmod_wald_st0(t, v, boundary, st0, what));\n",
           "  }\n  return acc;")))
}

read_prog <- function(fam, v) {
  f <- file.path(opt$out, fam, paste0(v, c(".stan", ".json")))
  if (!all(file.exists(f))) stop("no ", f[1], "; run `emit` first", call. = FALSE)
  list(stan = f[1], json = f[2], code = paste(readLines(f[1]), collapse = "\n"),
       sdat = jsonlite::read_json(f[2], simplifyVector = TRUE))
}

# ---- Points ------------------------------------------------------------------

# stanli's evaluation_point() (tools/benchmark_timer.hpp): point 0 a fixed
# pattern, points 1-9 uniform on (-2, 2) from splitmix64. R has no unsigned
# 64-bit integers, so they are four 16-bit limbs, least significant first.
u64 <- function(x) { out <- numeric(4); for (i in 1:4) { out[i] <- x %% 65536; x <- x %/% 65536 }; out }
u64_hex <- function(h) rev(strtoi(substring(h, seq(1, 13, 4), seq(4, 16, 4)), 16L))
u64_norm <- function(s) { for (i in 1:3) { s[i + 1] <- s[i + 1] + s[i] %/% 65536; s[i] <- s[i] %% 65536 }; s[4] <- s[4] %% 65536; s }
u64_add <- function(a, b) u64_norm(a + b)
u64_mul <- function(a, b) {
  s <- numeric(4)
  for (i in 1:4) for (j in 1:(5 - i)) s[i + j - 1] <- s[i + j - 1] + a[i] * b[j]
  u64_norm(s)
}
u64_xor <- function(a, b) as.numeric(bitwXor(as.integer(a), as.integer(b)))
u64_bits <- function(a) unlist(lapply(a, function(l) (l %/% 2^(0:15)) %% 2))
u64_shr <- function(a, k) {
  bits <- c(u64_bits(a)[-(1:k)], rep(0, k))
  vapply(1:4, function(i) sum(bits[(16 * i - 15):(16 * i)] * 2^(0:15)), 0)
}
splitmix64 <- function(x) {
  x <- u64_add(x, u64_hex("9E3779B97F4A7C15"))
  x <- u64_mul(u64_xor(x, u64_shr(x, 30)), u64_hex("BF58476D1CE4E5B9"))
  x <- u64_mul(u64_xor(x, u64_shr(x, 27)), u64_hex("94D049BB133111EB"))
  u64_xor(x, u64_shr(x, 31))
}
evaluation_point <- function(point, n) {
  i <- seq_len(n) - 1
  if (point == 0) return(0.1 + 0.05 * (i %% 7) - 0.15 * (i %% 3))
  base <- splitmix64(u64_mul(u64(point), u64_hex("9E3779B97F4A7C15")))
  vapply(i, function(k) {
    r <- u64_shr(splitmix64(u64_xor(base, u64(k))), 11)
    -2 + 4 * sum(r * 65536^(0:3)) * 2^-53
  }, 0)
}

# Point sets for a program, in the order of its unconstrained parameters
# (all of them unconstrained reals here, so the CmdStan fit only names them).
point_sets <- function(fc, fam) {
  pars <- fc$metadata()$stan_variables
  tr <- truth[[fam]][intersect(pars, names(truth[[fam]]))]
  tr$b <- array(tr$b, 1)
  q0 <- as.numeric(fc$unconstrain_variables(tr))
  set.seed(3)
  list(fixed = lapply(0:9, evaluation_point, n = length(q0)),
       post = lapply(1:40, function(i) q0 + stats::rnorm(length(q0), 0, 0.1)))
}

# The same "post" points without CmdStan: the parameters stanli's program
# declares, in its order, are what the corpus names.
stanli_points <- function(fam, v, n = 10) {
  code <- read_prog(fam, v)$code
  block <- regmatches(code, regexpr("(?s)\nparameters \\{.*?\n\\}", code, perl = TRUE))
  pars <- regmatches(block, gregexpr("(?<=real |vector\\[Kc\\] )[A-Za-z_]+", block, perl = TRUE))[[1]]
  q0 <- unlist(truth[[fam]][pars])
  set.seed(3)
  lapply(seq_len(n), function(i) q0 + stats::rnorm(length(q0), 0, 0.1))
}

# ---- Modes -------------------------------------------------------------------

scaled_err <- function(a, b) max(abs(a - b) / pmax(abs(a), abs(b), 1))

# Passes of the point set per block: about 0.4 s for the faster timer, since
# Windows' system.time() ticks at ~10 ms.
calibrate <- function(timers) {
  k <- 1L
  repeat {
    t <- min(vapply(timers, function(f) system.time(f(k))[["elapsed"]], 0))
    if (t >= 0.05 || k >= 2L^20) break
    k <- 2L * k
  }
  max(1L, as.integer(ceiling(k * 0.4 / max(t, 1e-3))))
}

grad <- function() {
  gp_load(".")
  suppressPackageStartupMessages(library(stanli))
  fam <- opt$family
  vs <- if (nzchar(opt$variants)) strsplit(opt$variants, ",")[[1]] else variant_names[[fam]]
  res <- list(); chk <- list()
  for (v in vs) {
    x <- read_prog(fam, v)
    fc <- gp_methods(gp_compile(x$stan, dir = tempdir()), x$json)$fit
    sm <- stanli::stanli_model(code = x$code, data = x$sdat, threads_per_chain = 1)
    sets <- point_sets(fc, fam)
    for (s in names(sets)) {
      P <- sets[[s]]
      err <- vapply(P, function(q) {
        a <- fc$grad_log_prob(q); b <- stanli::log_prob_grad(sm, q)
        if (!is.finite(attr(a, "log_prob")) && !is.finite(b$lp)) return(NA_real_)
        scaled_err(c(attr(a, "log_prob"), as.numeric(a)), c(b$lp, as.numeric(b[[2]])))
      }, 0)
      chk[[length(chk) + 1]] <- data.frame(family = fam, variant = v, set = s, points = length(P),
                                           nonfinite = sum(is.na(err)),
                                           max_scaled_err = max(err, na.rm = TRUE))
      timers <- list(cmdstan = function(k) for (i in seq_len(k)) for (q in P) fc$grad_log_prob(q),
                     stanli = function(k) for (i in seq_len(k)) for (q in P) stanli::log_prob_grad(sm, q))
      passes <- calibrate(timers)
      for (r in seq_len(opt$reps)) for (e in sample(names(timers))) {
        st <- system.time(timers[[e]](passes))[["elapsed"]]
        res[[length(res) + 1]] <- data.frame(family = fam, variant = v, set = s, rep = r, engine = e,
                                             structured_loops = loops,
                                             us = st / (passes * length(P)) * 1e6)
      }
      cat(sprintf("%-12s %-5s %d passes per block\n", v, s, passes))
    }
  }
  res <- do.call(rbind, res); chk <- do.call(rbind, chk)
  csv(res, "grad"); csv(chk, "check")
  w <- stats::reshape(res[, c("variant", "set", "rep", "engine", "us")],
                      idvar = c("variant", "set", "rep"), timevar = "engine", direction = "wide")
  tab <- do.call(rbind, lapply(split(w, list(w$variant, w$set), drop = TRUE), function(x) data.frame(
    variant = x$variant[1], set = x$set[1],
    cmdstan_us = signif(stats::median(x$us.cmdstan), 3), stanli_us = signif(stats::median(x$us.stanli), 3),
    stanli_over_cmdstan = signif(stats::median(x$us.stanli / x$us.cmdstan), 2))))
  tab <- merge(tab, chk[, c("variant", "set", "nonfinite", "max_scaled_err")], by = c("variant", "set"))
  cat("STANLI_STRUCTURED_LOOPS =", if (nzchar(loops)) loops else "(unset)", "\n")
  print(tab[order(match(tab$variant, vs), tab$set), ], row.names = FALSE)
}

replay <- function() {
  suppressPackageStartupMessages(library(stanli))
  x <- read_prog(opt$family, opt$variant)
  sm <- stanli::stanli_model(code = x$code, data = x$sdat, threads_per_chain = 1)
  P <- stanli_points(opt$family, opt$variant)
  if (opt$set == "fixed") P <- lapply(0:9, evaluation_point, n = length(P[[1]]))
  step <- P[[1]]; step[2] <- step[2] + 1e-3
  g <- function(q) stanli::log_prob_grad(sm, q)
  for (q in P) g(q)
  # Each pattern is one pass of 20 to 200 gradients, repeated to about 0.4 s.
  patterns <- list(
    "one point repeated" = function() for (i in 1:20) g(P[[1]]),
    "two points alternated" = function() for (i in 1:10) { g(P[[1]]); g(P[[2]]) },
    "a 1e-3 step alternated" = function() for (i in 1:10) { g(P[[1]]); g(step) },
    "ten points cycled" = function() for (i in 1:3) for (q in P) g(q),
    "ten points, each 3 times in a row" = function() for (q in P) for (i in 1:3) g(q),
    "ten points, each 20 times in a row" = function() for (q in P) for (i in 1:20) g(q))
  n <- c(20, 20, 20, 30, 30, 200)
  ms <- mapply(function(f, n) {
    k <- max(1L, as.integer(ceiling(0.4 / max(system.time(f())[["elapsed"]], 0.005))))
    system.time(for (i in seq_len(k)) f())[["elapsed"]] / (k * n) * 1e3
  }, patterns, n)
  out <- data.frame(family = opt$family, variant = opt$variant, set = opt$set,
                    structured_loops = loops, pattern = names(patterns), ms = signif(unname(ms), 3))
  csv(out, "replay")
  print(out, row.names = FALSE)
}

fit <- function() {
  suppressPackageStartupMessages(library(stanli))
  x <- read_prog(opt$family, opt$variant)
  sm <- stanli::stanli_model(code = x$code, data = x$sdat, threads_per_chain = 1)
  q0 <- stanli_points(opt$family, opt$variant, 1)[[1]]
  wall <- system.time(f <- stanli::sample_model(sm, chains = 1, seed = 1, warmup = opt$warmup,
                                                samples = 0, save_warmup = TRUE, init = q0,
                                                refresh = 0))[["elapsed"]]
  sp <- stanli::get_sampler_params(f, inc_warmup = TRUE)
  if (is.list(sp)) sp <- sp[[1]]
  out <- data.frame(family = opt$family, variant = opt$variant, structured_loops = loops,
                    iterations = nrow(sp), leapfrogs = sum(sp[, "n_leapfrog__"]), wall_s = wall,
                    ms_per_leapfrog = wall / sum(sp[, "n_leapfrog__"]) * 1e3)
  csv(out, "fit")
  print(out, row.names = FALSE)
}

ab <- function() {
  gp_load(".")
  fam <- opt$family
  m <- lapply(c(a = opt$a, b = opt$b), function(v) {
    x <- read_prog(fam, v)
    gp_methods(gp_compile(x$stan, dir = tempdir()), x$json)$fit
  })
  if (!identical(m$a$metadata()$stan_variables, m$b$metadata()$stan_variables)) {
    stop("`ab` needs two programs with the same parameters", call. = FALSE)
  }
  P <- point_sets(m$a, fam)$post
  timers <- lapply(m, function(fc) function(k) for (i in seq_len(k)) for (q in P) fc$grad_log_prob(q))
  passes <- calibrate(timers)
  res <- list()
  for (r in seq_len(opt$reps)) for (e in sample(names(timers))) {
    st <- system.time(timers[[e]](passes))[["elapsed"]]
    res[[length(res) + 1]] <- data.frame(rep = r, arm = e, us = st / (passes * length(P)) * 1e6)
  }
  res <- do.call(rbind, res)
  res$variant <- c(a = opt$a, b = opt$b)[res$arm]
  csv(res, "ab")
  w <- stats::reshape(res[, c("rep", "arm", "us")], idvar = "rep", timevar = "arm", direction = "wide")
  cat(sprintf("CmdStan: %s %.0f us, %s %.0f us, ratio b / a %.3f [%.3f, %.3f] over %d blocks\n",
              opt$a, stats::median(w$us.a), opt$b, stats::median(w$us.b),
              stats::median(w$us.b / w$us.a), min(w$us.b / w$us.a), max(w$us.b / w$us.a), nrow(w)))
}

diag <- function() {
  Sys.setenv(STANLI_REGION_MAP_DIAGNOSTICS = "1")
  suppressPackageStartupMessages(library(stanli))
  for (fam in names(variant_names)) for (v in variant_names[[fam]]) {
    x <- read_prog(fam, v)
    cat("=====", fam, v, "\n")
    invisible(stanli::stanli_model(code = x$code, data = x$sdat, threads_per_chain = 1))
  }
}

cat("stanli", as.character(utils::packageVersion("stanli")), "\n")
switch(mode, emit = emit(), grad = grad(), replay = replay(), fit = fit(), ab = ab(), diag = diag())
