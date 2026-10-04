# stanli vs CmdStan on the LNR: does it run, and is it any faster?
# =================================================================
#
# Question (2026-10-02): paul-buerkner/brms#1911 proposes stanli - an
# interpreter for Stan programs, no C++ toolchain - as a brms backend, and
# shows a "light" route that works with brms as released: brm(empty = TRUE),
# sample its Stan code and data with stanli, as_stanfit(), rename_pars(). Does
# that work for a cogmod family, and how does it compare with brms's cmdstanr
# backend? README.md has the answer; this script reproduces it.
#
# Model: the decision_making vignette's LNR on its data (speed_acc,
# participants 1-3, RT <= 2 s; 4620 trials), as in benchmarks/lnr_vectorize.
# The program brms writes for it is `orig`; the variants, each checked
# against it:
#
#   lse     `orig` with log_mix() written out as log_sum_exp(); the package's
#           density otherwise, branches and all. stanli 0.18.1 refused
#           log_mix() in a parameter branch, and this got round it
#   bf      a branch-free copy of the likelihood (lnr_branchfree_functions.stan,
#           sigmabias = 0 only) that drops the far-tail series of log Phi
#   sel     the same with the package's numerics kept, every branch written
#           as a select of clamped arms (lnr_select_functions.stan,
#           sigmabias = 0)
#   orig_s  `orig` with both sigmas as reals in the model block - one
#           softplus each instead of one per trial - and nothing else
#           changed: what brms would write if it hoisted intercept-only
#           dpars. The same log density as `orig`, bit for bit
#   rw      the hand-written program from cogmod#5 (lnr_rewrite.stan): rows
#           split by `dec` in transformed data, vector arithmetic per group,
#           each parameter branch a 0/1 blend of clamped arms - and scalar
#           sigmas, as `orig_s`
#   rw_v    `rw` with the sigmas as brms writes them, vector[N] through the
#           softplus: what a vectorised family (custom_family(loop = FALSE))
#           would get from the formula people write, `sigmazero ~ 1`
#
# orig -> orig_s and rw_v -> rw are the sigmas alone; orig -> rw_v and
# orig_s -> rw the likelihood alone.
#
# Run from the package root, in steps; each reads what the previous wrote.
#
#   Rscript benchmarks/stanli/bench.R emit [--out DIR]
#   Rscript benchmarks/stanli/bench.R grad [--out DIR] [--reps 21] [--tag NAME]
#                                          [--cmdstan v,...] [--stanli v,...]
#   Rscript benchmarks/stanli/bench.R fit  [--out DIR] [--rounds 1] [--iter 500] [--seed 11]
#                                          [--arms a,b,...] [--tag NAME]
#
#   emit  `orig`, the variants, data and an init
#   grad  log density and gradient of the --cmdstan variants (`orig` always)
#         in CmdStan and of the --stanli variants in stanli, at the init, 35
#         points around it and 9 edge points (results/check[_TAG].csv); the
#         same at a point and a response that put log Phi far into its tail
#         (results/tail[_TAG].csv); then the cost of one gradient,
#         alternating blocks (results/time[_TAG].csv)
#   fit   4 chains x (ITER warmup + ITER sampling) from the same inits; by
#         default {orig, sel} x {cmdstanr, stanli}: brm(backend = "cmdstanr")
#         on `orig`, cmdstanr on the program alone for any other variant, and
#         the light route on each. --arms picks others, e.g. stanli with one
#         process per chain (see the fit section). Writes
#         results/fit_summary[_TAG].csv, fit_chains, fit_pars
#
# DIR (generated programs, executables, fits) defaults to
# benchmarks/results/stanli. Needs stanli (seantalts.r-universe.dev, then
# stanli::stanli_install() once), rtdists (the data), cmdstanr and CmdStan.
# `grad` compiles three programs with model methods, a minute or two each on
# Windows; `fit` compiles `orig` and `sel` for sampling. Run nothing else on
# the machine while `grad` times or `fit` samples.
#
# How the cost is measured: R loops calling grad_log_prob() (cmdstanr model
# methods) and stanli's log_prob_grad() over 200 points, a different one each
# call, so nothing can be cached; programs alternate block by block, as in
# benchmarks/gradient_cost.R, and the headline is the ratio of medians within
# one run. Both go through R and a C++ boundary per call; at ~1-4 ms per
# gradient that overhead does not register. CPU time is recorded beside
# elapsed time, to show each uses one core.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load(), gp_compile(), gp_methods()

mode <- commandArgs(trailingOnly = TRUE)[1]
steps <- c("emit", "grad", "fit")
if (is.na(mode) || !mode %in% steps) {
  stop("usage: Rscript benchmarks/stanli/bench.R ", paste(steps, collapse = "|"),
       " [--options]", call. = FALSE)
}
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})
opt <- gp_args(list(out = "benchmarks/results/stanli", reps = 21L, rounds = 1L, iter = 500L,
                   seed = 11L, arms = "cmdstanr_orig,cmdstanr_sel,stanli_orig,stanli_sel",
                   tag = "", cmdstan = "orig,orig_s,rw,rw_v",
                   stanli = "orig,orig_s,rw,rw_v,sel"))
out <- opt$out
res_dir <- "benchmarks/stanli/results"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
dir.create(res_dir, recursive = TRUE, showWarnings = FALSE)
path <- function(...) file.path(out, paste0(...))

gp_load(".")
suppressPackageStartupMessages(library(stanli))

# ---- Data and model ----------------------------------------------------------

# benchmarks/lnr_vectorize/bench.R's lnr_data() and its `vignette` model.
lnr_data <- function() {
  data(speed_acc, package = "rtdists", envir = environment())
  df <- data.frame(
    Participant = as.integer(as.character(speed_acc$id)),
    Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[
      as.character(speed_acc$condition)]),
    RT = speed_acc$rt,
    Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat))
  )
  df[df$Participant %in% c(1, 2, 3) & df$RT <= 2, ]
}
lnr_formula <- function() {
  brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition, sigmazero ~ 1,
           sigmaone ~ 1, sigmabias = 0, ndt ~ Condition, family = cogmod_lnr())
}

# The brms program with its functions block swapped for one of ours, and
# sigmabias dropped from the likelihood call. Parameters, priors and
# generated quantities are untouched, so rename_pars() and every brms
# post-processing function see the program they expect.
swap_functions <- function(code, file) {
  x <- strsplit(code, "\n")[[1]]
  f0 <- grep("^functions [{]", x)
  d0 <- grep("^data [{]", x)
  y <- c(x[seq_len(f0 - 1)], readLines(file), x[d0:length(x)])
  call0 <- "sigmaone[n], sigmabias, ndt[n]"
  if (sum(grepl(call0, y, fixed = TRUE)) != 1) stop("likelihood call not found", call. = FALSE)
  paste(sub(call0, "sigmaone[n], ndt[n]", y, fixed = TRUE), collapse = "\n")
}

# The brms program with log_mix() - the one function stanli 0.18.1 refuses
# inside a branch on a parameter - written out. Nothing else changes.
logsumexp <- function(code) {
  old <- "return log_mix(poutlier, lp_out, lp_dec);"
  if (!grepl(old, code, fixed = TRUE)) stop("log_mix() call not found", call. = FALSE)
  sub(old, "return log_sum_exp(log(poutlier) + lp_out, log1m(poutlier) + lp_dec);",
      code, fixed = TRUE)
}

# Each `old` must occur exactly once in `code`; replaced by `new`, in order.
subs <- function(code, ...) {
  s <- list(...)
  for (i in seq(1, length(s), by = 2)) {
    n <- lengths(regmatches(code, gregexpr(s[[i]], code, fixed = TRUE)))
    if (n != 1) stop("expected one `", s[[i]], "`, found ", n, call. = FALSE)
    code <- sub(s[[i]], s[[i + 1]], code, fixed = TRUE)
  }
  code
}

# A top-level block of a Stan program ("parameters", "generated quantities"),
# from its header to the next header.
stan_block <- function(code, name) {
  x <- strsplit(code, "\n")[[1]]
  heads <- grep("^[a-z][a-z ]* [{]\\s*$", x)
  i <- heads[x[heads] == paste(name, "{")]
  if (length(i) != 1) stop("no block `", name, "`", call. = FALSE)
  j <- c(heads[heads > i], length(x) + 1)[1] - 1
  trimws(paste(x[i:j], collapse = "\n"))
}

# The sigmas as reals: brms writes vector[N] sigmazero = rep_vector(0.0, N),
# adds the intercept and takes the softplus of all N. As reals the arithmetic
# is the same, once, so the log density is unchanged bit for bit.
scalar_sigmas <- function(code) {
  subs(code,
       "vector[N] sigmazero = rep_vector(0.0, N);", "real sigmazero = 0;",
       "vector[N] sigmaone = rep_vector(0.0, N);", "real sigmaone = 0;",
       "sigmazero[n], sigmaone[n]", "sigmazero, sigmaone")
}

# cogmod#5's program, read as posted. It must keep `orig`'s data,
# parameters, priors and generated quantities, or neither the log densities
# nor rename_pars() would line up.
rewrite <- function(code) {
  rw <- paste(readLines("benchmarks/stanli/lnr_rewrite.stan"), collapse = "\n")
  for (b in c("data", "parameters", "transformed parameters", "generated quantities"))
    if (!identical(stan_block(rw, b), stan_block(code, b)))
      stop("lnr_rewrite.stan's `", b, "` block differs from brms's", call. = FALSE)
  rw
}

# `rw` with the sigmas as brms writes them for `sigmazero ~ 1`: a vector of
# N softplus values, sliced per group. The values are `rw`'s, so is the log
# density; the cost is N softplus and N extra entries on the tape per sigma,
# and the check-blend in vector form.
rewrite_vec <- function(code) {
  subs(rewrite(code),
       "vector nu_w, real s_w, vector nu_l,", "vector nu_w, vector s_w, vector nu_l,",
       "real s_l, vector ndt,", "vector s_l, vector ndt,",
       "vector[n] sw = s_w * ones;", "vector[n] sw = s_w;",
       "./ (s_l * ones), n)", "./ s_l, n)",
       "real sigmazero = log1p_exp(Intercept_sigmazero);",
       "vector[N] sigmazero = log1p_exp(rep_vector(0.0, N) + Intercept_sigmazero);",
       "real sigmaone = log1p_exp(Intercept_sigmaone);",
       "vector[N] sigmaone = log1p_exp(rep_vector(0.0, N) + Intercept_sigmaone);",
       "int bad = sigmazero <= 0 || sigmaone <= 0 ||",
       "int bad = min(sigmazero) <= 0 || min(sigmaone) <= 0 ||",
       "real sz = sigmazero * ok + bad;", "vector[N] sz = sigmazero * ok + bad;",
       "real so = sigmaone * ok + bad;", "vector[N] so = sigmaone * ok + bad;",
       "mu[idx0], sz, nuone[idx0], so,", "mu[idx0], sz[idx0], nuone[idx0], so[idx0],",
       "nuone[idx1], so, mu[idx1], sz,", "nuone[idx1], so[idx1], mu[idx1], sz[idx1],")
}

.VARIANTS <- list(
  orig = identity,
  lse = logsumexp,
  bf = function(code) swap_functions(code, "benchmarks/stanli/lnr_branchfree_functions.stan"),
  sel = function(code) swap_functions(code, "benchmarks/stanli/lnr_select_functions.stan"),
  orig_s = scalar_sigmas,
  rw = rewrite,
  rw_v = rewrite_vec
)

inputs <- function() readRDS(path("inputs.rds"))
code_of <- function(v) paste(readLines(path("lnr_", v, ".stan")), collapse = "\n")
stanli_of <- function(v, data) stanli_model(code = code_of(v), data = data, threads_per_chain = 1)
init_json <- function(init, file) {
  # write_stan_json() writes a length-1 R vector as a scalar, which CmdStan
  # rejects for the vector[1] coefficients (`b`, `b_nuone`, `b_ndt`) here; a
  # 1-d array keeps the brackets.
  vec <- grepl("^b(_|$)", names(init))
  init[vec] <- lapply(init[vec], function(v) array(v, dim = length(v)))
  cmdstanr::write_stan_json(init, file)
}
relerr <- function(g, ref) max(abs(g - ref) / pmax(1, abs(ref)))

# ---- emit --------------------------------------------------------------------
if (mode == "emit") {
  df <- lnr_data()
  f <- lnr_formula()
  prior <- suppressMessages(cogmod_priors(f, df))
  sv <- cogmod_stanvars(f)
  code <- as.character(brms::make_stancode(f, data = df, prior = prior, stanvars = sv))
  sdat <- lapply(unclass(brms::make_standata(f, data = df, prior = prior, stanvars = sv)), identity)
  for (v in names(.VARIANTS)) writeLines(.VARIANTS[[v]](code), path("lnr_", v, ".stan"))
  cmdstanr::write_stan_json(sdat, path("lnr.data.json"))
  init <- cogmod_inits(f, df, jitter = 0)(1)
  init_json(init, path("lnr.init.json"))
  saveRDS(list(df = df, f = f, prior = prior, sv = sv, sdat = sdat, init = init), path("inputs.rds"))
  cat("N =", sdat$N, "- wrote", out, "\n")
  e <- tryCatch({ stanli_model(code = code, data = sdat); "built" },
                error = function(e) conditionMessage(e))
  cat("stanli on the original program:", e, "\n")
}

# ---- grad --------------------------------------------------------------------
if (mode == "grad") {
  inp <- inputs()
  cm_v <- unique(c("orig", strsplit(opt$cmdstan, ",", fixed = TRUE)[[1]]))
  sl_v <- strsplit(opt$stanli, ",", fixed = TRUE)[[1]]
  unknown <- setdiff(c(cm_v, sl_v), names(.VARIANTS))
  if (length(unknown)) stop("unknown variant: ", unknown[1], call. = FALSE)
  res_file <- function(stem) file.path(res_dir, paste0(stem, if (nzchar(opt$tag)) paste0("_", opt$tag), ".csv"))
  cat("stanli", as.character(utils::packageVersion("stanli")), "\n")

  # CmdStan, with model methods. `orig` is what everything is compared with,
  # so it has to build; any other may go missing - Windows Defender
  # quarantined lnr_bf.exe straight after compiling, twice on 2026-10-03 (a
  # false positive it also makes on other CmdStan executables) - and is then
  # reported and left out.
  cat("compiling with model methods:", paste(cm_v, collapse = ", "), "\n")
  mods <- list()
  cm <- list()
  for (v in cm_v) {
    r <- tryCatch({
      mods[[v]] <- gp_compile(path("lnr_", v, ".stan"))
      gp_methods(mods[[v]], path("lnr.data.json"), path("lnr.init.json"))
    }, error = function(e) {
      if (v == "orig") stop(e)
      cat("CmdStan `", v, "` unavailable: ", conditionMessage(e), "\n", sep = "")
      NULL
    })
    if (!is.null(r)) cm[[v]] <- r
  }
  t0 <- Sys.time()
  sl <- list()
  for (v in sl_v) {
    m <- tryCatch(stanli_of(v, inp$sdat), error = function(e) conditionMessage(e))
    if (is.character(m)) cat("stanli on `", v, "`: ", m, "\n", sep = "") else sl[[v]] <- m
  }
  cat("stanli build,", length(sl), "programs:", round(as.numeric(Sys.time() - t0, units = "secs"), 1), "s\n")

  # One interface for both: q (unconstrained) -> list(lp, grad). stanli's lp
  # carries the Jacobian and the constants, as CmdStan's log_prob(jacobian =
  # TRUE) does. A program that throws gives NA and the message.
  lpg_cmdstan <- function(fit) function(q) {
    g <- fit$grad_log_prob(q)
    list(lp = attr(g, "log_prob"), grad = as.numeric(g))
  }
  lpg_stanli <- function(m) function(q) log_prob_grad(m, q)
  engines <- c(setNames(lapply(cm, function(x) lpg_cmdstan(x$fit)), paste0("cmdstan_", names(cm))),
               setNames(lapply(sl, lpg_stanli), paste0("stanli_", names(sl))))
  safely <- function(f, q) tryCatch(c(f(q), err = ""), error = function(e)
    list(lp = NA_real_, grad = rep(NA_real_, length(q)), err = conditionMessage(e)))

  # How far apart two log densities are, in units in the last place of the
  # larger: 0 is the same double. Infinities are 0 apart only from themselves.
  ulp <- function(a, b) {
    if (is.na(a) || is.na(b)) return(NA_real_)
    if (!is.finite(a) || !is.finite(b)) return(if (identical(a, b)) 0 else Inf)
    if (a == b) return(0)
    abs(a - b) / 2^(floor(log2(max(abs(a), abs(b)))) - 52)
  }
  relerr <- function(g, ref) max(abs(g - ref) / pmax(1, abs(ref)))
  compare <- function(pts, ref_name = "cmdstan_orig") {
    do.call(rbind, lapply(names(pts), function(p) {
      r <- lapply(engines, safely, q = pts[[p]])
      ref <- r[[ref_name]]
      do.call(rbind, lapply(names(r), function(e) {
        x <- r[[e]]
        both <- all(is.finite(x$grad)) && all(is.finite(ref$grad))
        data.frame(point = p, engine = e, lp = x$lp, lp_ulp = ulp(x$lp, ref$lp),
                   grad_finite = all(is.finite(x$grad)),
                   grad_relerr = if (both) relerr(x$grad, ref$grad) else NA_real_,
                   err = substr(x$err, 1, 120))
      }))
    }))
  }

  # Check: the init, 5 points near it, 30 further out, and edge points. Every
  # parameter is a scalar or a 1-vector here, so the skeleton's names index the
  # unconstrained vector.
  q0 <- cm$orig$up0
  sk <- cm$orig$fit$variable_skeleton(transformed_parameters = FALSE, generated_quantities = FALSE)
  stopifnot(all(lengths(sk) == 1), length(sk) == length(q0))
  ix <- setNames(seq_along(sk), names(sk))
  edge <- function(...) { q <- q0; a <- c(...); q[ix[names(a)]] <- a; q }
  ndt_hi <- q0[ix[["Intercept_ndt"]]] + 0.4
  set.seed(1)
  pts <- c(list(init = q0),
           setNames(lapply(1:5, function(i) q0 + stats::rnorm(length(q0), 0, 0.15)), paste0("near", 1:5)),
           setNames(lapply(1:30, function(i) q0 + stats::rnorm(length(q0), 0, 0.5)), paste0("wide", 1:30)),
           list(
             pout_tiny = edge(poutlier = -30),        # poutlier 1e-13
             pout_zero = edge(poutlier = -750),       # inv_logit() underflows: exactly 0
             pout_one = edge(poutlier = 40),          # rounds to exactly 1
             ndt_hi = edge(Intercept_ndt = ndt_hi),   # ndt 1.5x: more responses below it
             ndt_hi_pout_zero = edge(Intercept_ndt = ndt_hi, poutlier = -750),  # -inf
             ndt_lo = edge(Intercept_ndt = -6),       # ndt 2.5 ms
             sig_lo = edge(Intercept_sigmazero = -3, Intercept_sigmaone = -3),  # far tails
             sig_hi = edge(Intercept_sigmazero = 3, Intercept_sigmaone = 3),
             nu_far = edge(Intercept = 3, Intercept_nuone = -3)))
  check <- compare(pts)
  utils::write.csv(check, res_file("check"), row.names = FALSE)
  edges <- !grepl("^(init|near|wide)", check$point)
  smry <- function(x) data.frame(
    points = length(unique(x$point)), lp_identical = sum(x$lp_ulp == 0, na.rm = TRUE),
    lp_ulp_max = suppressWarnings(max(x$lp_ulp[is.finite(x$lp_ulp)])),
    lp_mismatch = sum(!is.finite(x$lp_ulp) | is.na(x$lp_ulp)),
    grad_relerr_max = suppressWarnings(max(x$grad_relerr, na.rm = TRUE)),
    grad_nonfinite = sum(!x$grad_finite), errors = sum(nzchar(x$err)))
  cat("\nInit, near and wide points, against CmdStan `orig`:\n")
  print(do.call(rbind, lapply(split(check[!edges, ], check$engine[!edges]), smry)), digits = 3)
  cat("\nEdge points:\n")
  print(check[edges, c("point", "engine", "lp", "lp_ulp", "grad_finite", "grad_relerr", "err")],
        digits = 6, row.names = FALSE)

  # The far tail, which cogmod_log_Phi()'s x < -25 series is there for. Both
  # sigmas at softplus(-3) = 0.049 and the first trial replaced by a 12 s
  # response (gradient_check.R's slowest): the loser's survival then sits
  # at x ~ -60, past the x = -38 where the erfc route underflows.
  sd12 <- inp$sdat
  sd12$Y[1] <- 12
  cmdstanr::write_stan_json(sd12, path("lnr12.data.json"))
  init12 <- inp$init
  init12$Intercept_sigmazero <- init12$Intercept_sigmaone <- -3
  init_json(init12, path("lnr12.init.json"))
  cm12 <- lapply(mods[names(cm)], gp_methods, data = path("lnr12.data.json"),
                 init = path("lnr12.init.json"))
  q12 <- cm12$orig$up0
  s_l <- log1p(exp(-3))
  nu_l <- if (sd12$dec[1] == 0) init12$Intercept_nuone else init12$Intercept
  cat(sprintf("\ntail check: loser's x on the 12 s trial = %.1f\n",
              (-nu_l - log(12 - exp(init12$Intercept_ndt))) / s_l))
  engines <- c(setNames(lapply(cm12, function(x) lpg_cmdstan(x$fit)), paste0("cmdstan_", names(cm12))),
               setNames(lapply(names(sl), function(v) lpg_stanli(stanli_of(v, sd12))), paste0("stanli_", names(sl))))
  tail <- compare(list(tail = q12))
  print(tail[, c("engine", "lp", "lp_ulp", "grad_finite", "grad_relerr", "err")], digits = 10, row.names = FALSE)
  utils::write.csv(tail, res_file("tail"), row.names = FALSE)

  # Cost.
  K <- 200
  set.seed(3)
  P <- lapply(seq_len(K), function(k) q0 + stats::rnorm(length(q0), 0, 0.1))
  timer <- function(f) function() for (k in seq_len(K)) f(P[[k]])
  timers <- c(setNames(lapply(cm, function(x) { fit <- x$fit; timer(fit$grad_log_prob) }), paste0("cmdstan_", names(cm))),
              setNames(lapply(sl, function(m) timer(function(q) log_prob_grad(m, q))), paste0("stanli_", names(sl))))
  for (nm in names(timers)) timers[[nm]]()  # warm up
  rows <- list()
  for (r in seq_len(opt$reps)) {
    for (nm in sample(names(timers))) {
      st <- system.time(timers[[nm]]())
      rows[[length(rows) + 1]] <- data.frame(rep = r, program = nm,
        us_per_grad = st[["elapsed"]] / K * 1e6,
        cpu_over_elapsed = (st[["user.self"]] + st[["sys.self"]]) / st[["elapsed"]])
    }
  }
  tm <- do.call(rbind, rows)
  utils::write.csv(tm, res_file("time"), row.names = FALSE)
  # us is the median over blocks; the ratio to CmdStan `orig` is taken within
  # each block, then its median and 10th to 90th percentile: steadier than a
  # ratio of medians when the machine's speed drifts during the run.
  med <- tapply(tm$us_per_grad, tm$program, stats::median)
  base <- tm[tm$program == "cmdstan_orig", ]
  rel <- tm$us_per_grad / base$us_per_grad[match(tm$rep, base$rep)]
  qs <- tapply(rel, tm$program, stats::quantile, probs = c(0.1, 0.5, 0.9))
  tab <- data.frame(program = names(med), us = round(med),
                    ratio = round(vapply(qs[names(med)], `[`, numeric(1), 2), 2),
                    p10 = round(vapply(qs[names(med)], `[`, numeric(1), 1), 2),
                    p90 = round(vapply(qs[names(med)], `[`, numeric(1), 3), 2),
                    cpu = round(tapply(tm$cpu_over_elapsed, tm$program, stats::median)[names(med)], 2))
  print(tab[order(tab$us), ], row.names = FALSE)
}

# ---- fit ---------------------------------------------------------------------
# Sampling efficiency. An arm is <engine>_<variant>[_proc | _rand], e.g.
# (--arms, comma-separated; the default is {orig, sel} x {cmdstanr, stanli}):
#
#   cmdstanr_orig     brm(backend = "cmdstanr"), the standard route
#   cmdstanr_sel      `sel` compiled by cmdstanr, sampled, read back as a
#                     stanfit; likewise cmdstanr_<any other variant>
#   stanli_orig       the light route from brms#1911 on the program brms writes
#   stanli_sel        the light route on `sel`, and so on
#   stanli_orig_proc  as stanli_orig, but one R process per chain
#   stanli_orig_rand  as stanli_orig, from stanli's own random inits (radius 2)
#
# Every arm but the _rand ones starts its four chains from the same
# cogmod_inits() and the same seed (--seed, which also draws the inits);
# stanli takes them on the unconstrained scale (unconstrain(), one row per
# chain). Given inits, stanli 0.19.0 and 0.19.1 run the chains one after
# another whatever parallel_chains says (repro_parallel_init.R), so the
# plain stanli arms pay the sum of their chains. The _proc arms get round it
# as CmdStan does, one process per chain: four PSOCK workers each sample
# chains = 1 from their own init row, at seed 100 * --seed + chain, and
# rstan::sflist2stanfit() merges the four. The _rand arms get round it the
# other way, with no inits, so their chains run on stanli's own threads.
#
# Compiling and building (the workers' too) happen before the clock starts, so
# `wall_s` is sampling to brmsfit. Each round runs the arms in turn, starting
# from a different one each round and each seed, so that across a set of seeds
# each arm goes first equally often; the draws repeat exactly across rounds,
# so extra rounds only show how much the timings move (--rounds, default 1).
# --tag names the output files (fit_summary_<tag>.csv, ...), so that runs with
# other arms do not overwrite these.
if (mode == "fit") {
  inp <- inputs()
  f <- inp$f; df <- inp$df; prior <- inp$prior; sv <- inp$sv
  arm_names <- strsplit(opt$arms, ",", fixed = TRUE)[[1]]
  # "stanli_rw_v_proc" -> engine "stanli", variant "rw_v", how "_proc".
  arm_parts <- lapply(setNames(arm_names, arm_names), function(a) {
    m <- regmatches(a, regexec("^(cmdstanr|stanli)_(.+?)(_proc|_rand)?$", a, perl = TRUE))[[1]]
    if (!length(m) || !m[3] %in% names(.VARIANTS) || (m[2] == "cmdstanr" && nzchar(m[4])))
      stop("unknown arm: ", a, call. = FALSE)
    list(engine = m[2], v = m[3], how = m[4])
  })
  variant_of <- function(a) arm_parts[[a]]$v
  set.seed(opt$seed)
  inits <- cogmod_inits(f, df)
  init_list <- lapply(1:4, function(i) inits(i))
  init_files <- vapply(1:4, function(i) {
    p <- path("fit_init", i, ".json"); init_json(init_list[[i]], p); p
  }, character(1))

  # backend = "cmdstanr" so that stancode() is the code brm() compiles below.
  dummy <- brms::brm(f, data = df, prior = prior, stanvars = sv, backend = "cmdstanr",
                     empty = TRUE)
  code_orig <- as.character(brms::stancode(dummy))
  vs <- unique(vapply(arm_names, variant_of, character(1)))
  code <- lapply(setNames(vs, vs), function(v) .VARIANTS[[v]](code_orig))
  sdat <- lapply(unclass(brms::standata(dummy)), identity)
  as_brmsfit <- function(stanfit) { d <- dummy; d$fit <- stanfit; brms:::rename_pars(d) }
  timed <- function(expr) { t0 <- Sys.time(); force(expr); as.numeric(Sys.time() - t0, units = "secs") }

  # brm() writes its program with cmdstanr::write_stan_file(), named by hash,
  # into this directory: compiling the same text here first leaves it an
  # executable that is up to date.
  options(cmdstanr_write_stan_file_dir = normalizePath(out))
  build <- list()
  cmod <- list()
  for (a in grep("^cmdstanr_", arm_names, value = TRUE)) {
    v <- variant_of(a)
    build[[a]] <- timed(cmod[[v]] <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code[[v]]), quiet = TRUE))
  }
  sm <- list()
  sl_arms <- grep("^stanli_", arm_names, value = TRUE)
  for (v in unique(vapply(sl_arms, variant_of, character(1)))) {
    b <- timed(sm[[v]] <- stanli_model(code = code[[v]], data = sdat))
    for (a in sl_arms[vapply(sl_arms, variant_of, character(1)) == v]) build[[a]] <- b
  }
  init_u <- lapply(sm, function(m) do.call(rbind, lapply(init_list, unconstrain, model = m)))

  proc <- grep("_proc$", arm_names, value = TRUE)
  if (length(proc)) {
    cl <- parallel::makePSOCKcluster(4)
    pv <- unique(vapply(proc, variant_of, character(1)))
    parallel::clusterExport(cl, c("code", "sdat", "pv"), envir = environment())
    b <- timed(parallel::clusterEvalQ(cl, {
      suppressPackageStartupMessages(library(stanli))
      wm <- lapply(setNames(pv, pv), function(v) stanli_model(code = code[[v]], data = sdat))
      NULL
    }))
    for (a in proc) build[[a]] <- b
  }

  n <- opt$iter
  stanli_arm <- function(v, init = init_u[[v]]) function() as_brmsfit(as_stanfit(sample_model(
    sm[[v]], chains = 4, seed = opt$seed, warmup = n, samples = n, init = init,
    parallel_chains = 4, refresh = 0)))
  proc_arm <- function(v) function() {
    sfs <- parallel::parLapply(cl, 1:4, function(i, v, init, seed, n) {
      stanli::as_stanfit(stanli::sample_model(wm[[v]], chains = 1, seed = 100 * seed + i,
                                              warmup = n, samples = n, init = init[i, ], refresh = 0))
    }, v = v, init = init_u[[v]], seed = opt$seed, n = n)
    as_brmsfit(rstan::sflist2stanfit(sfs))
  }
  cmdstanr_arm <- function(v) if (v == "orig") function() brms::brm(
    f, data = df, prior = prior, stanvars = sv, backend = "cmdstanr", chains = 4,
    cores = 4, warmup = n, iter = 2 * n, init = init_list, seed = opt$seed, refresh = 0,
    silent = 2) else function() {
      fc <- cmod[[v]]$sample(data = sdat, chains = 4, parallel_chains = 4, iter_warmup = n,
                             iter_sampling = n, init = init_files, seed = opt$seed, refresh = 0,
                             show_messages = FALSE)
      # brms's own reader, as brm(backend = "cmdstanr") uses: rstan::read_stan_csv()
      # does not parse CmdStan 2.38's CSV header ("object 'n_kept' not found").
      as_brmsfit(brms:::read_csv_as_stanfit(fc$output_files(), model = cmod[[v]]))
    }
  arms <- lapply(arm_parts, function(p) {
    if (p$engine == "cmdstanr") return(cmdstanr_arm(p$v))
    if (p$how == "_proc") proc_arm(p$v)
    else if (p$how == "_rand") stanli_arm(p$v, init = NULL)
    else stanli_arm(p$v)
  })

  fits <- list()
  chains <- list()
  for (r in seq_len(opt$rounds)) {
    k <- (opt$seed + r - 2) %% length(arms)
    arm_order <- names(arms)[(seq_along(arms) + k - 1) %% length(arms) + 1]
    for (a in arm_order) {
      cat("seed", opt$seed, "round", r, "-", a, "\n")
      t0 <- Sys.time()
      fit <- arms[[a]]()
      wall <- as.numeric(Sys.time() - t0, units = "secs")
      el <- rstan::get_elapsed_time(fit$fit)
      np <- brms::nuts_params(fit)
      per <- function(p) as.numeric(tapply(np$Value[np$Parameter == p], np$Chain[np$Parameter == p], sum))
      chains[[length(chains) + 1]] <- data.frame(
        seed = opt$seed, round = r, arm = a, position = match(a, arm_order),
        chain = seq_len(nrow(el)), warmup_s = el[, "warmup"],
        sample_s = el[, "sample"], leapfrog = per("n_leapfrog__"),
        divergent = per("divergent__"), build_s = build[[a]], wall_s = wall)
      fits[[a]] <- fit
    }
  }
  chains <- do.call(rbind, chains)
  if (length(proc)) parallel::stopCluster(cl)
  out_name <- function(stem, ext) paste0(stem, if (nzchar(opt$tag)) paste0("_", opt$tag), ext)
  saveRDS(fits, path(out_name("fits", ".rds")))

  pars <- do.call(rbind, lapply(names(fits), function(a) {
    d <- posterior::summarise_draws(
      brms::as_draws_df(fits[[a]], variable = "^b_|^poutlier", regex = TRUE),
      "mean", "sd", "rhat", "ess_bulk", "ess_tail")
    data.frame(seed = opt$seed, arm = a, as.data.frame(d))
  }))

  # Per arm and round. ESS per second is the smallest bulk ESS over the nine
  # parameters divided by the slowest chain's warmup + sampling: what the fit
  # would take if its chains ran in parallel. ESS per wall second divides by
  # what it did take, sampling to brmsfit.
  summ <- do.call(rbind, lapply(split(chains, list(chains$round, chains$arm), drop = TRUE), function(x) {
    p <- pars[pars$arm == x$arm[1], ]
    chain_s <- x$warmup_s + x$sample_s
    data.frame(seed = opt$seed, round = x$round[1], arm = x$arm[1], position = x$position[1],
               build_s = x$build_s[1],
               chain_s_mean = mean(chain_s), chain_s_max = max(chain_s), wall_s = x$wall_s[1],
               leapfrog = sum(x$leapfrog), ms_per_leapfrog = 1000 * sum(chain_s) / sum(x$leapfrog),
               ess_bulk_min = min(p$ess_bulk), ess_tail_min = min(p$ess_tail),
               ess_bulk_per_s = min(p$ess_bulk) / max(chain_s),
               ess_bulk_per_wall_s = min(p$ess_bulk) / x$wall_s[1],
               ess_bulk_per_1k_leapfrog = 1000 * min(p$ess_bulk) / sum(x$leapfrog),
               rhat_max = max(p$rhat), divergent = sum(x$divergent))
  }))
  summ <- summ[order(summ$round, match(summ$arm, names(arms))), ]
  if ("cmdstanr_orig" %in% summ$arm) {
    i0 <- match(summ$round, summ$round[summ$arm == "cmdstanr_orig"])
    summ$ess_per_s_vs_cmdstanr_orig <- summ$ess_bulk_per_s / summ$ess_bulk_per_s[summ$arm == "cmdstanr_orig"][i0]
    summ$ess_per_wall_s_vs_cmdstanr_orig <- summ$ess_bulk_per_wall_s / summ$ess_bulk_per_wall_s[summ$arm == "cmdstanr_orig"][i0]
  }

  utils::write.csv(chains, file.path(res_dir, out_name("fit_chains", ".csv")), row.names = FALSE)
  utils::write.csv(pars, file.path(res_dir, out_name("fit_pars", ".csv")), row.names = FALSE)
  utils::write.csv(summ, file.path(res_dir, out_name("fit_summary", ".csv")), row.names = FALSE)
  print(summ, digits = 3, row.names = FALSE)
  print(pars, digits = 4)
  # The same posterior every way: elpd differences should all be ~0.
  print(loo::loo_compare(lapply(fits, brms::loo)))
}
