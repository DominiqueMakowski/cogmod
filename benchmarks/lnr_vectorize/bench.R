# Vectorising cogmod_lnr(): is it feasible, and is it any faster?
# =================================================================
#
# Question (2026-09-24): brms can write a custom family's likelihood as ONE
# vectorised call instead of a loop over observations (custom_family(loop =
# FALSE), there since brms 2.16). Is that worth doing for the LNR, which is
# simple but slow on real data sets? README.md has the answer; this script
# reproduces it.
#
# Variants of the likelihood, all on the same brms program otherwise:
#
#   base   the package as it is: loop = TRUE, one scalar call per trial
#   v0     loop = FALSE, the loop moved inside the function (scalar calls)
#   v1     loop = FALSE, one fused loop, no nested user-function calls
#   v2     loop = FALSE, vector arithmetic over index sets (winner 0, winner
#          1, t <= ndt), cogmod_log_Phi() vectorised by branch
#   v3     loop = TRUE again - NOT vectorised - with the loser's survival as
#          one lognormal_lccdf() (guarded by cogmod_log_Phi()'s own switch at
#          z = 25) instead of ~9 separate operations
#
# on four models of the vignette's data (speed_acc, participants 1-3,
# RT <= 2 s; 4620 trials):
#
#   vignette    the decision_making vignette's formula: sigmas `~ 1`
#   scalars     the same with the sigmas left out of bf() (brms aux reals)
#   hier        random intercepts and slopes on both nus, intercepts on ndt
#   startpoint  sigmabias ~ 1: the general (sigmabias > 0) density
#
# Run from the package root, in steps; each reads what the previous wrote.
#
#   Rscript benchmarks/lnr_vectorize/bench.R emit    [--out DIR]
#   Rscript benchmarks/lnr_vectorize/bench.R compile [--out DIR] [--workers 4]
#   Rscript benchmarks/lnr_vectorize/bench.R check   [--out DIR]
#   Rscript benchmarks/lnr_vectorize/bench.R profile [--out DIR]
#   Rscript benchmarks/lnr_vectorize/bench.R time    [--out DIR] [--reps 21]
#            [--models vignette,scalars,hier,startpoint] [--variants base,v0,v1,v2,v3]
#   Rscript benchmarks/lnr_vectorize/bench.R threads [--out DIR] [--reps 11]
#
# DIR defaults to benchmarks/results/lnr_vectorize. Needs rtdists (the data),
# cmdstanr and CmdStan. `compile` builds 26 programs, about 1-2 minutes each
# on Windows; `time` takes about an hour at --reps 21, and `check` spends most
# of its time in the startpoint model's warmup. Run nothing else on the
# machine while `time` or `threads` runs.
#
# How the cost is measured: NUTS with adaptation off, at the step size and
# diagonal metric adapted by `check`'s short base fit, from its last draw.
# CmdStan's own sampling timer divided by the leapfrog steps taken is the
# cost of one gradient, as a fit pays it. Variants alternate block by block,
# as in benchmarks/gradient_cost.R; the headline is the ratio of medians
# within one run. Absolute times moved 15% between runs on the laptop the
# README's numbers come from; ratios did not.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()
source("benchmarks/lnr_vectorize/proto.R")

mode <- commandArgs(trailingOnly = TRUE)[1]
steps <- c("emit", "compile", "check", "profile", "time", "threads")
if (is.na(mode) || !mode %in% steps) {
  stop("usage: Rscript benchmarks/lnr_vectorize/bench.R ",
       paste(steps, collapse = "|"), " [--options]", call. = FALSE)
}
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})
opt <- gp_args(list(out = "benchmarks/results/lnr_vectorize", workers = 4L, reps = 0L,
                    models = "vignette,scalars,hier,startpoint",
                    variants = "base,v0,v1,v2,v3"))
out <- opt$out
dir.create(out, recursive = TRUE, showWarnings = FALSE)
path <- function(...) file.path(out, paste0(...))

# ---- Data and models ---------------------------------------------------------

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

# The same, salted with the responses gradient_check.R adds: the fastest
# possible trial and 3, 6 and 12 s ones, both responses, both conditions.
lnr_tails <- function(df) {
  tails <- expand.grid(Participant = 1L, Condition = c("Accuracy", "Speed"),
                       RT = c(min(df$RT) + 1e-4, 3, 6, 12), Error = 0:1,
                       stringsAsFactors = FALSE)
  rbind(df, tails)
}

.MODELS <- list(
  vignette = function(fam) brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition,
    sigmazero ~ 1, sigmaone ~ 1, sigmabias = 0, ndt ~ Condition, family = fam),
  scalars = function(fam) brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition,
    sigmabias = 0, ndt ~ Condition, family = fam),
  hier = function(fam) brms::bf(RT | dec(Error) ~ Condition + (Condition | Participant),
    nuone ~ Condition + (Condition | Participant), sigmazero ~ 1, sigmaone ~ 1,
    sigmabias = 0, ndt ~ Condition + (1 | Participant), family = fam),
  startpoint = function(fam) brms::bf(RT | dec(Error) ~ Condition, nuone ~ Condition,
    sigmazero ~ 1, sigmaone ~ 1, sigmabias ~ 1, ndt ~ Condition, family = fam)
)

# NUTS iterations per timing block. The two models whose base fit adapts a
# small step size in 200 warmup iterations run to treedepth 10, about a
# minute per 30 iterations; two iterations still give ~1300 gradients.
.ITERS <- c(vignette = 30L, scalars = 30L, hier = 2L, startpoint = 2L)

# profile() sections around brms's linear predictors, the likelihood and the
# priors, in a model block brms wrote. Works for the looped likelihood and for
# the single vectorised call.
add_profiles <- function(x) {
  m0 <- grep("^model [{]", x)
  decl <- max(grep("^    vector[[]N[]] [a-z]+ = rep_vector[(]0.0, N[)];", x))
  lik <- grep("target [+]= cogmod_lnr_lpdf", x)
  looped <- grepl("for [(]n in 1:N[)]", x[lik - 1])
  a <- if (looped) lik - 1 else lik
  b <- if (looped) lik + 1 else lik
  pri <- grep("^  target [+]= lprior;", x)
  end <- m0 + grep("^[}]", x[(m0 + 1):length(x)])[1]
  c(x[1:decl], '    profile("linpred") {', x[(decl + 1):(a - 1)], "    }",
    '    profile("lik") {', x[a:b], "    }",
    x[(b + 1):(pri - 1)], '  profile("prior") {', x[pri:(end - 1)], "  }", x[end:length(x)])
}

# ---- emit --------------------------------------------------------------------
if (mode == "emit") {
  gp_load(".")
  df <- lnr_data()
  df_t <- lnr_tails(df)
  for (nm in names(.MODELS)) {
    f0 <- .MODELS[[nm]](cogmod_lnr())
    prior <- suppressMessages(cogmod_priors(f0, df))
    sd0 <- brms::make_standata(f0, data = df, prior = prior)
    cmdstanr::write_stan_json(lapply(unclass(sd0), identity), path(nm, ".data.json"))
    sdt <- brms::make_standata(f0, data = df_t, prior = prior)
    cmdstanr::write_stan_json(lapply(unclass(sdt), identity), path(nm, ".tails.data.json"))
    code0 <- brms::make_stancode(f0, data = df, prior = prior, stanvars = cogmod_stanvars(f0))
    # write_stan_json() writes a length-1 vector as a scalar; keep declared vectors 1-D
    init <- cogmod_inits(f0, df, jitter = 0)(1)
    vecs <- regmatches(code0, gregexpr("\n  vector(<[^>]*>)?\\[[^]]+\\] ([A-Za-z0-9_]+);", code0))[[1]]
    vecs <- sub(";$", "", sub("^.* ", "", vecs))
    for (p in intersect(vecs, names(init))) init[[p]] <- array(init[[p]], dim = length(init[[p]]))
    cmdstanr::write_stan_json(init, path(nm, ".init.json"))

    writeLines(code0, path(nm, "__base.stan"))
    f1 <- .MODELS[[nm]](lnr_vec_family())
    for (st in c("v0", "v1", "v2")) {
      writeLines(brms::make_stancode(f1, data = df, prior = prior, stanvars = lnr_vec_stanvars(f1, st)),
                 path(nm, "__", st, ".stan"))
    }
    writeLines(brms::make_stancode(f0, data = df, prior = prior, stanvars = lnr_v3_stanvars()),
               path(nm, "__v3.stan"))
    if (nm != "startpoint") writeLines(add_profiles(strsplit(code0, "\n")[[1]]), path(nm, "__prof.stan"))
    cat(nm, ": N =", sd0$N, "; vectors in the loop = FALSE call:",
        paste(names(which(lnr_arg_types(f1))), collapse = ", "), "\n")
  }
  for (st in c("v2", "v3")) {
    writeLines(add_profiles(readLines(path("vignette__", st, ".stan"))), path("vignette__", st, "prof.stan"))
  }
  # The vignette model, base family, with brms within-chain threading.
  f0 <- .MODELS$vignette(cogmod_lnr())
  prior <- suppressMessages(cogmod_priors(f0, df))
  thr <- brms::threading(4, grainsize = 250)
  writeLines(brms::make_stancode(f0, data = df, prior = prior, stanvars = cogmod_stanvars(f0), threads = thr),
             path("vignette__thr.stan"))
  sd <- brms::make_standata(f0, data = df, prior = prior, threads = thr)
  cmdstanr::write_stan_json(lapply(unclass(sd), identity), path("vignette.thr.data.json"))
  quit(status = 0)
}

# ---- compile -----------------------------------------------------------------
if (mode == "compile") {
  stan <- list.files(out, pattern = "[.]stan$", full.names = TRUE)
  cl <- parallel::makeCluster(opt$workers)
  on.exit(parallel::stopCluster(cl))
  res <- parallel::parLapply(cl, stan, function(f) {
    t0 <- Sys.time()
    cpp <- if (grepl("__thr[.]stan$", f)) list(stan_threads = TRUE) else list()
    cmdstanr::cmdstan_model(f, cpp_options = cpp, force_recompile = TRUE, quiet = TRUE)
    sprintf("%-28s %4.0f s", basename(f), as.numeric(Sys.time() - t0, units = "secs"))
  })
  cat(unlist(res), sep = "\n")
  quit(status = 0)
}

suppressPackageStartupMessages(library(cmdstanr))
# The executables need tbb.dll, which cmdstanr puts on PATH only for its own calls.
if (.Platform$OS.type == "windows") {
  Sys.setenv(PATH = paste(normalizePath(file.path(cmdstan_path(), "stan/lib/stan_math/lib/tbb")),
                          Sys.getenv("PATH"), sep = ";"))
}
exe <- function(nm, v) path(nm, "__", v, ".exe")
models <- strsplit(opt$models, ",")[[1]]
variants <- strsplit(opt$variants, ",")[[1]]

# ---- check -------------------------------------------------------------------
# Every variant against the base: log density and gradient, from CmdStan's own
# log_prob method at 20 post-warmup draws of the base model, on the real data
# and on the tail-salted data. Also leaves behind the adapted step size,
# metric and last draw that `profile`, `time` and `threads` start from.
if (mode == "check") {
  lp_grad <- function(nm, v, draws_csv, data_json) {
    o <- tempfile(fileext = ".csv")
    res <- system2(exe(nm, v), c("log_prob", paste0("constrained_params=", draws_csv), "jacobian=1",
                                 "data", paste0("file=", data_json), "output", paste0("file=", o),
                                 "sig_figs=18"), stdout = TRUE, stderr = TRUE)
    if (!file.exists(o)) stop(nm, "/", v, ": ", paste(utils::tail(res, 5), collapse = "\n"))
    as.matrix(utils::read.csv(o, comment.char = "#"))
  }
  rows <- list()
  for (nm in models) {
    mod <- cmdstan_model(exe_file = exe(nm, "base"))
    fit <- mod$sample(data = path(nm, ".data.json"), init = path(nm, ".init.json"),
                      chains = 1, iter_warmup = 200, iter_sampling = 20, seed = 1,
                      refresh = 0, show_messages = FALSE, output_dir = out,
                      output_basename = paste0(nm, "__draws"))
    draws_csv <- fit$output_files()
    last <- fit$draws(format = "draws_df")[20, ]
    saveRDS(list(step = fit$metadata()$step_size_adaptation,
                 inv_metric = fit$inv_metric()[[1]], last = last), path(nm, ".adapt.rds"))
    # the last draw as an init, in the shapes the JSON init had (draws are column-major)
    raw <- jsonlite::fromJSON(path(nm, ".init.json"), simplifyVector = FALSE)
    lst <- posterior::as_draws_list(last)[[1]]
    init <- lapply(names(raw), function(p) {
      v <- unname(unlist(lst[grep(paste0("^", p, "(\\[|$)"), names(lst), value = TRUE)]))
      if (!is.list(raw[[p]])) v
      else if (is.list(raw[[p]][[1]])) matrix(v, nrow = length(raw[[p]]))
      else array(v, dim = length(v))
    })
    names(init) <- names(raw)
    write_stan_json(init, path(nm, ".last.json"))

    for (dat in c("data", "tails")) {
      dj <- path(nm, if (dat == "data") ".data.json" else ".tails.data.json")
      b <- lp_grad(nm, "base", draws_csv, dj)
      g <- grep("^g_", colnames(b))
      for (v in setdiff(variants, "base")) {
        x <- lp_grad(nm, v, draws_csv, dj)
        rows[[length(rows) + 1]] <- data.frame(model = nm, data = dat, variant = v,
          lp_base = b[1, "lp__"], max_abs_dlp = max(abs(x[, "lp__"] - b[, "lp__"])),
          max_rel_dgrad = max(abs(x[, g] - b[, g]) / pmax(1, abs(b[, g]))),
          finite = all(is.finite(x)))
      }
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(out, "check.csv"), row.names = FALSE)
  print(res, digits = 3, row.names = FALSE)
  quit(status = 0)
}

# One NUTS run with adaptation off; returns ms per gradient and the leapfrog count.
run_fixed <- function(mod, nm, data_json, iters, seed, threads = NULL) {
  ad <- readRDS(path(nm, ".adapt.rds"))
  fit <- mod$sample(data = data_json, init = path(nm, ".last.json"), chains = 1,
    threads_per_chain = threads, iter_warmup = 0, iter_sampling = iters,
    adapt_engaged = FALSE, step_size = ad$step, inv_metric = diag(ad$inv_metric),
    seed = seed, refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
  nl <- sum(fit$sampler_diagnostics(format = "draws_matrix")[, "n_leapfrog__"])
  list(fit = fit, ms = 1000 * fit$time()$chains$sampling / nl, leapfrog = nl)
}

# ---- profile -----------------------------------------------------------------
# Share of the gradient per section, and autodiff tape entries per trial.
if (mode == "profile") {
  progs <- c(vignette = "prof", scalars = "prof", hier = "prof",
             vignette = "v2prof", vignette = "v3prof")
  rows <- list()
  for (k in seq_along(progs)) {
    nm <- names(progs)[k]
    N <- jsonlite::fromJSON(path(nm, ".data.json"))$N
    r <- run_fixed(cmdstan_model(exe_file = exe(nm, progs[[k]])), nm, path(nm, ".data.json"),
                   .ITERS[[nm]], 1)
    pr <- r$fit$profiles()[[1]]
    rows[[k]] <- data.frame(model = nm, program = progs[[k]], section = pr$name,
      share = pr$total_time / sum(pr$total_time), forward_share = pr$forward_time / pr$total_time,
      tape_per_trial = (pr$chain_stack + pr$no_chain_stack) / pr$autodiff_calls / N,
      ns_per_trial = 1e9 * pr$total_time / pr$autodiff_calls / N)
  }
  res <- do.call(rbind, rows)
  res <- res[res$section != "prior", ]
  utils::write.csv(res, file.path(out, "profile.csv"), row.names = FALSE)
  print(res, digits = 3, row.names = FALSE)
  quit(status = 0)
}

# ---- time --------------------------------------------------------------------
if (mode == "time") {
  reps <- if (opt$reps > 0) opt$reps else 21L
  rows <- list()
  for (nm in models) {
    vs <- intersect(variants, c("base", "v0", "v1", "v2", "v3"))
    mods <- stats::setNames(lapply(vs, function(v) cmdstan_model(exe_file = exe(nm, v))), vs)
    one <- function(v, seed) run_fixed(mods[[v]], nm, path(nm, ".data.json"), .ITERS[[nm]], seed)
    for (v in vs) one(v, 1)  # warm the file cache
    for (r in seq_len(reps)) {
      for (v in if (r %% 2) vs else rev(vs)) {
        x <- one(v, 100 + r)
        rows[[length(rows) + 1]] <- data.frame(model = nm, variant = v, rep = r,
                                               ms = x$ms, leapfrog = x$leapfrog)
      }
    }
    res <- do.call(rbind, rows[vapply(rows, function(x) x$model == nm, TRUE)])
    med <- tapply(res$ms, res$variant, stats::median)[vs]
    cat(sprintf("%-11s %-5s %.3f ms/gradient  ratio to base %.3f\n", nm, vs, med,
                med / med[["base"]]), sep = "")
    utils::write.csv(do.call(rbind, rows), file.path(out, "time.csv"), row.names = FALSE)
  }
  quit(status = 0)
}

# ---- threads -----------------------------------------------------------------
# The vignette model with the package as it is, plain vs brms-threaded
# (reduce_sum, grainsize 250) at 1, 2 and 4 threads: wall time per gradient.
if (mode == "threads") {
  reps <- if (opt$reps > 0) opt$reps else 11L
  m_base <- cmdstan_model(exe_file = exe("vignette", "base"))
  # Built from the .stan file (which finds the compiled executable up to date)
  # so that cmdstanr knows it has stan_threads; from exe_file alone it does
  # not, and silently drops threads_per_chain.
  m_thr <- cmdstan_model(path("vignette__thr.stan"), cpp_options = list(stan_threads = TRUE), quiet = TRUE)
  cfg <- list(base = list(m_base, ".data.json", NULL), thr1 = list(m_thr, ".thr.data.json", 1),
              thr2 = list(m_thr, ".thr.data.json", 2), thr4 = list(m_thr, ".thr.data.json", 4))
  one <- function(k, seed) {
    c <- cfg[[k]]
    run_fixed(c[[1]], "vignette", path("vignette", c[[2]]), .ITERS[["vignette"]], seed, c[[3]])$ms
  }
  for (k in names(cfg)) one(k, 1)
  rows <- list()
  for (r in seq_len(reps)) {
    for (k in if (r %% 2) names(cfg) else rev(names(cfg))) {
      rows[[length(rows) + 1]] <- data.frame(config = k, rep = r, ms = one(k, 100 + r))
    }
  }
  res <- do.call(rbind, rows)
  utils::write.csv(res, file.path(out, "threads.csv"), row.names = FALSE)
  med <- tapply(res$ms, res$config, stats::median)[names(cfg)]
  cat(sprintf("%-5s %.3f ms/gradient  speed-up vs base %.2f\n", names(cfg), med,
              med[["base"]] / med), sep = "")
}
