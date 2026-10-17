# Step 2 of the LNR start check (README.md), part 2: one cell, model x arm x
# seed. Arm A starts `nuone` at 0.7 (the 0.3.5 value, a 50/50 race); arm B
# sets the registry entry to 0 in-session, and nothing else differs: same
# data, same program, same priors, same seed - and because the jitter is
# additive on the link scale and drawn in the same order, the same jitter on
# every parameter, so the two arms' chains differ only in the centre of
# `Intercept_nuone`. The starts actually used go to
# `results/<model>__<arm>__s<seed>__inits.csv`.
#
# Chains are drawn with CmdStan directly, warmup kept, 250 + 250 iterations,
# 4 chains, diag_e, as the driftone prior check did.
#
# On the cluster (task.slurm) this is a script:
#
#   Rscript fit_cell.R <row of cells.csv>
#   Env: LS_DIR (the pushed directory: cells.csv, spec.rds, tree/, results/),
#        LS_WORK (this task's scratch directory: the program is compiled and
#        CmdStan writes there)
#
# Locally, fit_arm.R sources it and calls fit_cell() per cell.

set_arm <- function(arm) {
  ns <- asNamespace("cogmod")
  ch <- get(".CHOICE", envir = ns)
  ch$cogmod_lnr$init$nuone <- c(A = 0.7, B = 0)[[arm]]
  if (bindingIsLocked(".CHOICE", ns)) unlockBinding(".CHOICE", ns)
  assign(".CHOICE", ch, envir = ns)
  get(".CHOICE", envir = ns)$cogmod_lnr$init$nuone
}

fit_cell <- function(name, arm, seed, spec, stan_dir, out_dir, resdir,
                     n_warm = 250L, n_samp = 250L) {
  log_msg <- function(...) {
    cat(format(Sys.time(), "%H:%M:%S"), name, arm, seed, ..., "\n")
  }
  stem <- file.path(resdir, sprintf("%s__%s__s%d", name, arm, seed))
  if (file.exists(paste0(stem, "__chains.csv"))) {
    log_msg("exists, skipping")
    return(invisible(NULL))
  }
  nuone_start <- set_arm(arm)
  log_msg("cogmod", as.character(utils::packageVersion("cogmod")),
          "brms", as.character(utils::packageVersion("brms")),
          "cmdstan", cmdstanr::cmdstan_version(), "nuone start", nuone_start)

  sp <- spec[[name]]
  f <- sp$f
  d <- sp$data
  pri <- suppressMessages(suppressWarnings(cogmod::cogmod_priors(f, d)))
  sv <- cogmod::cogmod_stanvars(f)
  code <- suppressWarnings(brms::make_stancode(f, data = d, prior = pri,
                                               stanvars = sv))
  sdata <- suppressWarnings(brms::make_standata(f, data = d, prior = pri,
                                                stanvars = sv))
  # Same text as prep.R wrote locally, so an up-to-date executable is reused;
  # on the cluster each cell compiles in its own scratch directory.
  file <- cmdstanr::write_stan_file(code, dir = stan_dir, basename = name)
  mod <- cmdstanr::cmdstan_model(file, quiet = TRUE)

  # The starts, one list per chain, drawn before sampling so that they can
  # be recorded. cmdstanr accepts a list of per-chain lists.
  set.seed(seed)
  initf <- cogmod::cogmod_inits(f, d)
  inits <- lapply(1:4, initf)
  st <- do.call(rbind, lapply(seq_along(inits), function(i) {
    v <- inits[[i]]
    v <- v[lengths(v) == 1]
    data.frame(chain = i, par = names(v), value = unlist(v, use.names = FALSE))
  }))
  utils::write.csv(st, paste0(stem, "__inits.csv"), row.names = FALSE)
  # Off the seed again before sampling. processx < 3.9 tags each chain
  # process with an id drawn by sample() plus the current second, and
  # cmdstanr's cleanup (on.exit, kill_tree()) SIGKILLs every process on the
  # machine carrying that tag. The two arms of a seed consume the RNG
  # identically, so when they started sampling on one node in the same
  # second the first to finish killed the other's running chains: 16 of 96
  # chains in array 11596401, 20 of 96 in 11596580, all exit status -9, all
  # in same-node pairs. CmdStan's own seed is passed explicitly, so this
  # changes nothing in the draws.
  set.seed(NULL)

  # CmdStan writes into a directory of this cell's own, kept if a chain is
  # incomplete.
  cell_dir <- file.path(out_dir, basename(stem))
  unlink(cell_dir, recursive = TRUE)
  dir.create(cell_dir, recursive = TRUE)
  log_msg("sampling")
  fit <- tryCatch(mod$sample(
    data = unclass(sdata), seed = seed, chains = 4, parallel_chains = 4,
    iter_warmup = n_warm, iter_sampling = n_samp, init = inits,
    refresh = 0, show_messages = FALSE, show_exceptions = FALSE,
    save_warmup = TRUE, metric = "diag_e", output_dir = cell_dir
  ), error = function(e) {
    log_msg("sample() failed:", conditionMessage(e))
    NULL
  })

  # Which chains finished is read off CmdStan's CSVs, not off cmdstanr, so
  # that a chain lost for any reason shows up rather than vanishing: cmdstanr
  # drops a chain it did not see finish (it goes by the "Elapsed Time" lines
  # it parsed, not by the exit status) and renumbers the rest 1..k, which in
  # array 11596401 misaligned the per-chain times and hid the chains killed
  # there (see set.seed(NULL) above). So each chain's exit status,
  # cmdstanr's verdict and the draws in its CSV are recorded, the times come
  # from the CSV's footer, and the fit is rebuilt from the complete CSVs.
  files <- if (!is.null(fit)) {
    tryCatch(fit$output_files(include_failed = TRUE), error = function(e) NULL)
  }
  if (length(files) == 0L) files <- list.files(cell_dir, pattern = "\\.csv$", full.names = TRUE)
  files <- files[file.exists(files) & !grepl("-profile-", basename(files))]
  # cmdstanr names them <model>-<timestamp>-<chain id>-<hash>.csv
  ids <- as.integer(sub("^.*-([0-9]+)-[0-9a-zA-Z]+\\.csv$", "\\1", basename(files)))
  footer <- function(l, section) {
    x <- grep(sprintf("seconds (%s)", section), l, fixed = TRUE, value = TRUE)
    if (length(x)) as.numeric(sub("^.*[^0-9.]([0-9.]+) seconds.*$", "\\1", x[1])) else NA_real_
  }
  st <- data.frame(chain = 1:4)
  st$rc <- if (!is.null(fit)) fit$return_codes()[1:4] else NA_integer_
  st$cmdstanr_ok <- if (!is.null(fit)) {
    st$chain %in% ids[files %in% fit$output_files()]
  } else FALSE
  st$draws <- 0L
  st$time_warm <- st$time_samp <- st$time_total <- NA_real_
  for (k in seq_along(files)) {
    l <- readLines(files[k])
    j <- match(ids[k], st$chain)
    st$draws[j] <- sum(grepl("^[-0-9]", l))
    st$time_warm[j] <- footer(l, "Warm-up")
    st$time_samp[j] <- footer(l, "Sampling")
    st$time_total[j] <- footer(l, "Total")
  }
  st$complete <- st$draws == n_warm + n_samp
  bad <- !st$complete | (!is.na(st$rc) & st$rc != 0)

  # The evidence for any chain that is incomplete, or that cmdstanr did not
  # count: its console output as cmdstanr read it, and an incomplete CSV.
  fd <- file.path(resdir, "failed")
  for (i in st$chain[bad | !st$cmdstanr_ok]) {
    dir.create(fd, showWarnings = FALSE)
    f <- sprintf("%s__chain%d", file.path(fd, basename(stem)), i)
    if (!is.null(fit)) writeLines(c(sprintf("exit status %s", st$rc[i]), fit$output(i)),
                                  paste0(f, ".txt"))
    if (bad[i] && any(ids == i)) file.copy(files[ids == i], paste0(f, ".csv"), overwrite = TRUE)
  }
  ok <- files[match(st$chain[st$complete], ids)]
  if (!length(ok)) {
    log_msg("no complete chain; rc", st$rc, "draws", st$draws)
    utils::write.csv(st, paste0(stem, "__chains.csv"), row.names = FALSE)
    return(invisible(TRUE))
  }
  cf <- cmdstanr::as_cmdstan_fit(ok, check_diagnostics = FALSE)
  # .chain k of the rebuilt fit is the k-th complete chain
  idk <- st$chain[st$complete]

  sdg <- as.data.frame(cf$sampler_diagnostics(inc_warmup = TRUE,
                                              format = "draws_df"))
  lp <- as.data.frame(cf$draws("lp__", inc_warmup = TRUE,
                               format = "draws_df"))
  per <- do.call(rbind, lapply(seq_along(idk), function(k) {
    x <- sdg[sdg$.chain == k, ]
    i <- idk[k]
    w <- x$.iteration <= n_warm
    l <- lp$lp__[lp$.chain == k]
    lw <- l[w]
    ls <- l[!w]
    # The cold-start transient, as benchmarks/inits_ablation measured it:
    # the first warmup iteration at which lp__ reaches the chain's own
    # post-warmup 5th percentile.
    thr <- stats::quantile(ls, 0.05, names = FALSE)
    transient <- which(lw >= thr)[1]
    # CmdStan's console output for the chain, in which a rejected start
    # shows up as "Rejecting initial value".
    out <- if (!is.null(fit)) fit$output(i) else character()
    data.frame(
      chain = i,
      leap_warm = sum(x$n_leapfrog__[w]), leap_samp = sum(x$n_leapfrog__[!w]),
      div_warm = sum(x$divergent__[w]), div = sum(x$divergent__[!w]),
      td10_warm = sum(x$treedepth__[w] >= 10), td10 = sum(x$treedepth__[!w] >= 10),
      stepsize = x$stepsize__[!w][1],
      lp_start = lw[1], lp_mean = mean(ls),
      transient = if (is.na(transient)) n_warm + 1L else transient,
      rejected = sum(grepl("Rejecting initial value", out))
    )
  }))
  per$ebfmi <- tryCatch(cf$diagnostic_summary(quiet = TRUE)$ebfmi,
                        error = function(e) NA_real_)
  chains <- merge(st, per, by = "chain", all.x = TRUE)
  chains$frozen <- chains$stepsize < 0.01 * stats::median(chains$stepsize, na.rm = TRUE)
  chains$node <- Sys.info()[["nodename"]]

  vars <- cf$metadata()$model_params
  vars <- vars[grepl("^(b_|b\\[|Intercept|poutlier|sd_)", vars) &
                 !grepl("^(lprior|lp__|z_|r_)", vars)]
  dd <- as.data.frame(cf$draws(vars, format = "draws_df"))
  pars <- as.data.frame(posterior::summarise_draws(
    posterior::as_draws_df(dd[, vars, drop = FALSE]), median,
    ~stats::quantile(.x, probs = c(0.05, 0.95)), "rhat", "ess_bulk",
    "ess_tail"))
  spread <- vapply(vars, function(v) {
    m <- tapply(dd[[v]], dd$.chain, stats::median)
    diff(range(m)) / stats::sd(dd[[v]])
  }, numeric(1))
  pars$chain_spread <- spread[pars$variable]
  utils::write.csv(pars, paste0(stem, "__pars.csv"), row.names = FALSE)
  utils::write.csv(chains, paste0(stem, "__chains.csv"), row.names = FALSE)
  log_msg(sprintf(
    "done in %.0fs: chains complete %d (exit status %s; cmdstanr counted %d), div %d, frozen %d, rejected %d, max rhat %.3f, min ess_bulk %.0f, leapfrog warm %d samp %d",
    max(chains$time_total, na.rm = TRUE), sum(chains$complete),
    paste(chains$rc, collapse = " "), sum(chains$cmdstanr_ok),
    sum(chains$div, na.rm = TRUE), sum(chains$frozen, na.rm = TRUE),
    sum(chains$rejected, na.rm = TRUE), max(pars$rhat, na.rm = TRUE),
    min(pars$ess_bulk, na.rm = TRUE), sum(chains$leap_warm, na.rm = TRUE),
    sum(chains$leap_samp, na.rm = TRUE)))
  if (!any(bad)) unlink(cell_dir, recursive = TRUE)
  # TRUE when a chain is incomplete, so that task.slurm keeps the scratch
  # directory and the task shows as failed.
  invisible(any(bad))
}

if (sys.nframe() == 0L) {
  idx <- as.integer(commandArgs(trailingOnly = TRUE)[1])
  DIR <- Sys.getenv("LS_DIR")
  WORK <- Sys.getenv("LS_WORK")
  cells <- utils::read.csv(file.path(DIR, "cells.csv"), stringsAsFactors = FALSE)
  suppressMessages(pkgload::load_all(file.path(DIR, "tree"), quiet = TRUE))
  suppressMessages(library(brms))
  spec <- readRDS(file.path(DIR, "spec.rds"))
  resdir <- file.path(DIR, "results")
  dir.create(resdir, showWarnings = FALSE)
  bad <- fit_cell(cells$model[idx], cells$arm[idx], cells$seed[idx], spec,
                  stan_dir = WORK, out_dir = WORK, resdir = resdir)
  quit(status = if (isTRUE(bad)) 3L else 0L)
}
