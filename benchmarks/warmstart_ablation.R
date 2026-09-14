# Benchmark: what does each half of a warm start buy - the metric, or the step
# size - and what happens to the whole of it when the source no longer matches
# the target?
#
# Companion to the "Improving Sampling Efficiency and Performance" article
# (vignettes/articles/performance.qmd), "Warm Starts" section. Not part of the
# R package; run it from the repository root:
#
#   Rscript benchmarks/warmstart_ablation.R
#
# or from an R session, overriding any of the settings in `defaults` below
# without editing the file:
#
#   bench_overrides <- list(families = "lnr", seeds = 1); source("benchmarks/warmstart_ablation.R")
#
# The first question. cogmod_warmstart() hands a new fit three things at once:
# the previous fit's adapted inverse metric, its step size, and its posterior
# means as starting values. warm_start_subset.R shows the three together are
# worth several times a cold short warmup, and that the starting values alone
# are worth nothing. That leaves the interesting comparison unmade: of the two
# arguments that do the work, cogmod_inv_metric() and cogmod_step_size(), does
# one carry the effect on its own?
#
# They are not obviously separable. The step size Stan can use is a property
# of the metric it is using - a metric that rescales the awkward directions
# lets much longer steps stay stable - so a step size taken from a fit whose
# metric is not supplied may be far too large (divergences, then re-adaptation
# down) and a well-scaled metric with Stan's default step size of 1 has to
# find its own way down anyway. The plausible outcomes are that the metric
# carries nearly all of it, or that neither alone does much and the pair is
# what matters.
#
# The second question. Every warm start assumes the source resembles the
# target. The usual case is the one above, a pilot on some of the participants
# of the final data set, and there the assumption is safe by construction. But
# a warm start may also be taken from last year's study, or from a fit of the
# same model to a different population, and then the metric, step size,
# starting values and (if used) priors may be anywhere from slightly to badly
# off. How gracefully does each degrade?
#
# Design. Per family: ONE pilot fit on `n_subset` participants with a full
# `warmup_pilot` warmup, which is the warm-start source for everything below.
# Then THREE target data sets, each of `n_participants` participants:
#
#   pilot_included    the pilot's participants plus as many others - the
#                     standard setup, where the pilot is a subset of the target
#   new_participants  the same number of participants, none of whom the pilot
#                     saw: same task, same population, so the population-level
#                     parameters are about the same, but not one group-level
#                     entry of the metric carries over by level
#   shifted           new_participants with every RT transformed to
#                     rt_shift + rt_scale * RT, mostly a shift: the
#                     non-decision time of a slower population, roughly
#                     doubled and many posterior SDs from the pilot's, while
#                     the decision-time parameters barely move (the scale is
#                     1.1). The formula is the same; one parameter is not, and
#                     it is the one every RT model has
#
# On the first target, eight runs, all with `sample_iter` sampling iterations
# and everything else held fixed:
#
#   reference           cold inits, warmup_full - what a normal fit of this
#                       model costs
#   base                cold inits, warmup_short - the same short warmup as
#                       every cell below
#   step_size           base + cogmod_step_size()
#   inv_metric          base + cogmod_inv_metric()
#   both                base + inv_metric + step_size
#   both + inits        both, plus cogmod_inits(warmstart = ) - the full warm
#                       start, everything cogmod_warmstart() extracts
#   all + pilot priors  both + inits, plus cogmod_priors(warmstart = <the
#                       pilot>) - the full warm start with the pilot's
#                       posterior as the prior on top: everything the pilot
#                       can possibly give the new fit
#   oracle priors       base + cogmod_priors(warmstart = <the reference run>)
#
# `step_size`, `inv_metric` and `oracle priors` each add exactly one thing to
# `base` at a fixed short warmup, with inits held at the cold ones, so that
# each helper is measured on its own. `reference` says what the whole exercise
# is competing against; `both` and `both + inits` assemble the helpers, and say
# whether the warm inits add anything once the sampler is already well set up
# (warm_start_subset.R says they should not); `all + pilot priors` asks whether
# priors add anything on top of the assembled warm start.
#
# On the two transfer targets, the runs in `runs_transfer`: the same minus the
# oracle, whose source (a full fit of the same target) has nothing to do with
# the pilot. The question there is how each piece of the warm start, and the
# assembled whole with and without the pilot's starting values and priors,
# fares against a target the pilot does not describe. The single-helper cells
# earn their place here because of the RDM: a transferred metric with cold
# starting values can leave a chain stuck, and whether that also happens on
# participants the pilot never saw is exactly what these columns are for.
#
# The two priors cells are not like the others, and their numbers are not
# speedups you can have. cogmod_priors(warmstart = ) re-centres every
# population-level prior on `normal(median, prior_scale * sd)` of a previous
# fit's posterior, which changes the model rather than the sampler.
#
#   oracle priors       is centred on the *reference run of the same target* -
#                       a full cold fit of the very data the cell then fits
#                       again. The prior therefore already knows the answer,
#                       from those same data: maximal double counting (the
#                       data enter once as prior and once as likelihood), the
#                       intervals are too narrow, and nothing in the sampler's
#                       diagnostics reveals it. Nobody can do this in practice,
#                       which is why it is called an oracle. What the cell
#                       measures is an upper bound: how fast these models
#                       sample when the priors are as informative and as
#                       well-placed as they could possibly be. If that bound
#                       is not far above `inv_metric`, informative priors are
#                       not where the speed is. Only on the first target.
#   all + pilot priors  is centred on the pilot, a different fit. On
#                       `pilot_included` that still double-counts the pilot's
#                       participants (half the data, once as prior, once as
#                       likelihood). On `new_participants` it is the one
#                       legitimate use of the argument - an independent sample
#                       of the same population - and the cell says what it is
#                       worth. On `shifted` it is a wrong prior stated with
#                       whatever confidence the pilot had, and the cell says
#                       how much damage that does: watch `z_shift`.
#
# Two consequences for reading the table. A priors cell samples a different
# posterior from every other cell, so its `ess_per_s` is not comparable like
# for like: an effective draw from a narrower posterior is a different unit.
# The `sd_ratio` and `z_shift` columns of the summary are there to say how
# different - see below. And because a prior is part of the Stan program,
# those cells compile their own executables, which every other cell of the
# family avoids by reusing the pilot's.
#
# Families: LNR, DDM (the 4-parameter one - sigmadrift, sigmabias and sigmandt
# fixed at 0, see helpers.R), RDM and LBA2, each in the "mixed" structure -
# participant random intercepts on two or three parameters. With
# n_participants = 10 that is 10 standardized effects per term, so most of the
# metric's entries are `z` entries, which is the case cogmod_warmstart() has
# to get right.
#
# What it records, per fit (one row of the CSV):
#   dataset, family, run, seed
#   warmup_s, sample_s, wall_s   Stan's own per-chain times; wall is the
#                                slowest chain, warmup + sampling
#   ess_bulk_min, ess_tail_min   smallest ESS over all parameters, `z_` aside
#   ess_per_s                    ess_bulk_min / wall_s: the headline figure
#   ess_per_grad                 ess_bulk_min / post-warmup leapfrog steps -
#                                hardware-independent, and the one to read if
#                                the machine was doing something else at 3am
#   leapfrog_mean, treedepth_mean
#   stepsize_mean                mean post-warmup step size. With `step_size`
#                                supplied but not the metric, watch this fall
#                                back during warmup - that is the mechanism
#   divergent, rhat_max, n_pars
#
# and, in the summary rather than per fit, two columns that compare the
# posterior each run produced against the reference run's *of the same
# target*, over the population-level parameters (`b_*` and `sd_*`, written per
# fit to warmstart_ablation_postsd.csv):
#   sd_ratio        median over parameters of this run's posterior SD over the
#                   reference's. About 1 for every run that only changed the
#                   sampler - that is the check that they did only change the
#                   sampler - and below 1 for the priors cells, by the amount
#                   the informative prior sharpened the posterior
#   z_shift         largest move of a posterior mean away from the reference's,
#                   in units of the reference's posterior SD. Small everywhere
#                   if the runs agree on the answer; a large value in a
#                   sampler-only run means that run did not converge (check
#                   rhat_max: a chain stuck at max treedepth shows up here),
#                   and a large value in `all + pilot priors` on `shifted` is
#                   the wrong prior pulling the answer
#
# Runtime. Per family: 1 pilot + 8 runs on the first target + 7 on each of the
# other two, times length(seeds) - 44 fits with two seeds - plus the Stan
# compilations: the pilot's, which every default-prior run reuses through
# update(), and one per distinct prior specification (the oracle priors differ
# by seed, the pilot priors by target). A first pass at 150 trials per
# participant and 500 sampling iterations took 45 min (LNR) and 2 h (DDM) for
# the first target alone, and the RDM about 8 h for all three; the defaults
# below - 80 trials, 300 sampling iterations, 600 warmup for the reference and
# the pilot - cost roughly a third of that, so the four families should fit in
# one night, the LBA2 being much the slowest. `datasets = "pilot_included"`
# recovers the original, single-target benchmark; `families` and `seeds` are
# the other knobs. Changing any of the size settings makes the existing
# results CSV incomparable - move it out of `out_dir` first, or resume will
# skip cells that were fitted under the old settings.
#
# If it dies partway, `resume = TRUE` (the default) reads the CSV back and
# skips the (dataset, family, run, seed) rows already in it. A CSV written by
# the single-target version of this script, before the `dataset` column
# existed, is read as `pilot_included`. The pilot fit of an unfinished family
# is redone - it is the source the helpers need and is not stored between
# sessions - so nothing is lost but its own wall time. The reference run's
# warm-start table, which the `oracle priors` cell needs, IS stored (one file
# per target, family and seed), so a resumed session does not have to refit
# it; if the file is missing that cell is skipped with a note rather than
# silently falling back to the default priors.
#
# Output: benchmarks/results/warmstart_ablation.csv (one row per fit),
# _summary.csv (means over seeds, each cell relative to its target's `base`)
# and a figure of the ESS per second of every run - one row of panels per
# family, one column per target - saved to man/figures/ for the performance
# article. The figure is drawn by warmstart_ablation_figure.R, which can also
# be run on its own to redraw it from the CSV. The log, the pilot and
# reference tables and the postsd CSV are gitignored; the two CSVs above and
# the figure are what gets committed.

source("benchmarks/helpers.R")

# ---- Settings -----------------------------------------------------------------

defaults <- list(
  families = c("lnr", "ddm", "rdm", "lba2"),
  n_participants = 10,                # participants per target, drawn from rtdists::speed_acc
  n_subset = 5,                       # of the first target's, the pilot's share
  n_trials = 80,                      # trials per participant, balanced across conditions.
                                      # Few for a real study, but the cost of a fit is
                                      # linear in them and the geometry is not
  chains = 4,
  cores = 4,
  sample_iter = 300,                  # post-warmup iterations, the same for every run.
                                      # 1200 draws: enough to estimate a minimum ESS,
                                      # and ESS/s is a rate, so more only costs time
  warmup_short = 200,                 # the warmup a warm start is meant to make enough
  warmup_full = 600,                  # the reference run's: a properly adapted cold fit
  warmup_pilot = 600,                 # the source's. Its own knob, since the pilot is
                                      # not one of the models being compared - lower it
                                      # to ask what a cheap pilot is worth instead
  seeds = 1:2,                        # more seeds = less noisy ratios, proportionally longer
  datasets = c("pilot_included", "new_participants", "shifted"),
  runs = c("reference", "base", "step_size", "inv_metric", "both", "both + inits",
           "all + pilot priors", "oracle priors"),                 # on pilot_included
  runs_transfer = c("reference", "base", "step_size", "inv_metric", "both", "both + inits",
                    "all + pilot priors"),                            # on the other two
  rt_shift = 0.3,                     # the `shifted` target: RT -> rt_shift + rt_scale * RT,
  rt_scale = 1.1,                     # in seconds. Moves every ndt up by ~0.3 s and
                                      # barely touches the decision-time scale (10%)
  prior_scale = 3,                    # the priors cells' prior SD, as a multiple of
                                      # the source posterior's. 1 would use the
                                      # source posterior outright
  resume = TRUE,                      # skip rows already in the results CSV
  out_dir = file.path("benchmarks", "results")
)
settings <- if (exists("bench_overrides")) modifyList(defaults, bench_overrides) else defaults
list2env(settings, envir = environment())
dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)
results_file <- file.path(out_dir, "warmstart_ablation.csv")
post_file <- file.path(out_dir, "warmstart_ablation_postsd.csv")

runs_for <- function(dataset) if (identical(dataset, "pilot_included")) runs else runs_transfer

# The reference run's adaptation and posterior summaries, kept per target,
# family and seed so that the `oracle priors` cell can be centred on its own
# seed's reference and a resumed session need not refit it. The first target
# keeps the file name the single-target version of this script used.
ref_file <- function(dataset, family, seed) {
  tag <- if (identical(dataset, "pilot_included")) "" else paste0("_", dataset)
  file.path(out_dir, sprintf("warmstart_ablation_%s%s_reference_seed%s.csv", family, tag, seed))
}

# ---- Data -------------------------------------------------------------------------
# The pilot's participants are the first n_subset of the first target's. The
# second target is drawn from the participants the pilot did not see (same
# seed, different pool), the third is the second with its RTs transformed.

df <- bench_data(n_participants, n_trials)
sub_ids <- levels(df$Participant)[seq_len(n_subset)]
df_sub <- droplevels(df[df$Participant %in% sub_ids, ])
cat(sprintf("Pilot subset: %d participants, %d trials\n", n_subset, nrow(df_sub)))

targets <- list(pilot_included = df)
if (any(datasets != "pilot_included")) {
  df_new <- bench_data(n_participants, n_trials, exclude = as.integer(sub_ids))
  stopifnot(!any(levels(df_new$Participant) %in% sub_ids))
  targets$new_participants <- df_new
  df_shift <- df_new
  df_shift$RT <- rt_shift + rt_scale * df_shift$RT
  targets$shifted <- df_shift
  cat(sprintf("Transfer targets: %d participants none of which are the pilot's; ",
              n_participants),
      sprintf("shifted RTs = %.2f + %.2f * RT (median %.2f s, was %.2f s)\n",
              rt_shift, rt_scale, median(df_shift$RT), median(df_new$RT)))
}
targets <- targets[datasets]

# Rows from a previous, interrupted run. Keyed on the four things that
# identify a fit; the pilot rows are keyed too, but never skipped, since the
# pilot is refitted whenever its family has anything left to do. A CSV from
# before the `dataset` column existed is the first target's.
key <- function(dataset, family, run, seed) paste(dataset, family, run, seed, sep = " | ")
with_dataset <- function(d) {
  if (!is.null(d) && is.null(d$dataset)) d <- cbind(dataset = "pilot_included", d)
  d
}
previous <- if (isTRUE(resume) && file.exists(results_file)) {
  p <- with_dataset(read.csv(results_file, stringsAsFactors = FALSE))
  cat(sprintf("Resuming: %d rows already in %s\n", nrow(p), results_file))
  # The same cell under different size settings is a different measurement.
  stale <- p$run != "pilot" & (p$n_obs != nrow(df) | p$sample_iter != sample_iter)
  if (any(stale)) {
    stop(sprintf(paste0("%d rows of %s were fitted with other settings (n_obs or ",
                        "sample_iter differ). Move the old results out of `out_dir` ",
                        "before rerunning, or set resume = FALSE to overwrite them."),
                 sum(stale), results_file), call. = FALSE)
  }
  p
} else NULL
done <- if (is.null(previous)) character() else key(previous$dataset, previous$family, previous$run, previous$seed)
results <- if (is.null(previous)) list() else list(previous)
posts <- if (isTRUE(resume) && file.exists(post_file)) {
  list(with_dataset(read.csv(post_file, stringsAsFactors = FALSE)))
} else list()

# ---- Measurement ----------------------------------------------------------------

measure <- function(fit, fallback_wall) {
  vars <- variables(fit)
  vars <- vars[!vars %in% c("lp__", "lprior") & !startsWith(vars, "z_")]
  s <- summarise_draws(as_draws_array(fit, variable = vars), "rhat", "ess_bulk", "ess_tail")

  np <- nuts_params(fit)  # post-warmup only
  by_par <- function(p) np$Value[np$Parameter == p]
  n_leapfrog <- sum(by_par("n_leapfrog__"))

  # Stan's per-chain timings (columns warmup, sample). Chains run in parallel,
  # so the wall time of a run is its slowest chain; the two components are
  # summed over chains, which is what tells warmup cost from sampling cost.
  t <- tryCatch(rstan::get_elapsed_time(fit$fit), error = function(e) NULL)
  if (!is.null(t) && !all(is.finite(t))) t <- NULL
  ess_min <- min(s$ess_bulk, na.rm = TRUE)
  wall_s <- if (is.null(t)) fallback_wall else max(rowSums(t))

  data.frame(
    warmup_s = if (is.null(t)) NA_real_ else sum(t[, "warmup"]),
    sample_s = if (is.null(t)) NA_real_ else sum(t[, "sample"]),
    wall_s = wall_s,
    ess_bulk_min = ess_min,
    ess_tail_min = min(s$ess_tail, na.rm = TRUE),
    ess_per_s = ess_min / wall_s,
    ess_per_grad = ess_min / n_leapfrog,
    leapfrog_mean = mean(by_par("n_leapfrog__")),
    treedepth_mean = mean(by_par("treedepth__")),
    stepsize_mean = mean(by_par("stepsize__")),
    divergent = sum(by_par("divergent__")),
    rhat_max = max(s$rhat, na.rm = TRUE),
    # Dimension of the space Stan samples in, i.e. the length the inverse
    # metric has to have. The b_*_Intercept columns are recomputed from the
    # centred `Intercept*` parameters, not sampled twice.
    n_pars = sum(!grepl("^b_(.*_)?Intercept$", vars))
  )
}

record <- function(dataset, family, run, seed, participants, n_obs, warmup, fit, fallback_wall) {
  row <- cbind(
    data.frame(dataset = dataset, family = family, run = run, seed = seed,
               participants = participants, n_obs = n_obs,
               chains = chains, warmup = warmup, sample_iter = sample_iter),
    measure(fit, fallback_wall)
  )
  print(row[, c("wall_s", "warmup_s", "sample_s", "ess_bulk_min", "ess_per_s",
                "ess_per_grad", "stepsize_mean", "divergent", "rhat_max")],
        row.names = FALSE, digits = 3)
  results[[length(results) + 1]] <<- row
  write.csv(do.call(rbind, results), results_file, row.names = FALSE)

  # The posterior itself, not just how fast it arrived: enough to say later
  # whether a run that was only supposed to change the sampler changed the
  # answer, and by how much a priors cell narrowed or moved it. Long format,
  # appended, one row per population-level parameter.
  pop <- grep("^(b|sd)_", variables(fit), value = TRUE)
  if (length(pop)) {
    ps <- summarise_draws(as_draws_array(fit, variable = pop), "mean", "sd")
    posts[[length(posts) + 1]] <<- data.frame(
      dataset = dataset, family = family, run = run, seed = seed,
      variable = ps$variable, mean = ps$mean, sd = ps$sd
    )
    write.csv(do.call(rbind, posts), post_file, row.names = FALSE)
  }
}

# ---- One cell -----------------------------------------------------------------------
# The arguments each cell passes to brm(), on top of the ones held fixed. Only
# `init`, `warmup`, `prior` and the presence of `inv_metric` / `step_size`
# vary. `tab` is the pilot's own warm-start table, and every helper maps it
# onto the target `dfx` describes - which is the point: the same source, three
# targets.

run_args <- function(run, dataset, family, f, tab, seed, dfx) {
  # cogmod_inits() returns a function, which brms calls once per chain, so one
  # object serves every cold cell. The helpers are wrapped in closures because
  # switch() evaluates only the arm it selects: the `reference` and `base`
  # cells then never touch the warm start at all.
  cold <- list(init = cogmod_inits(f, dfx), warmup = warmup_short)
  warm <- function() list(init = cogmod_inits(f, dfx, warmstart = tab), warmup = warmup_short)
  metric <- function() list(inv_metric = cogmod_inv_metric(f, dfx, warmstart = tab))
  step <- function() list(step_size = cogmod_step_size(f, dfx, warmstart = tab))
  # The cells that change the model. `stanvars` is passed because a new prior
  # means a new Stan program, so update() has to recompile, and it needs the
  # family's functions block to do that. The oracle's source is this seed's
  # reference run on this target, written to disk when that run finished;
  # without the file there is nothing to centre on, and NULL tells the loop to
  # skip rather than quietly run the cell on the default priors and label it
  # `oracle priors`.
  oracle_priors <- function() {
    rf <- ref_file(dataset, family, seed)
    if (!file.exists(rf)) return(NULL)
    list(prior = cogmod_priors(f, dfx, warmstart = rf, prior_scale = prior_scale),
         stanvars = cogmod_stanvars(f))
  }
  pilot_priors <- function() {
    list(prior = cogmod_priors(f, dfx, warmstart = tab, prior_scale = prior_scale),
         stanvars = cogmod_stanvars(f))
  }
  switch(run,
    "reference"          = list(init = cold$init, warmup = warmup_full),
    "base"               = cold,
    "step_size"          = c(cold, step()),
    "inv_metric"         = c(cold, metric()),
    "both"               = c(cold, metric(), step()),
    "both + inits"       = c(warm(), metric(), step()),
    "all + pilot priors" = c(warm(), metric(), step(), pilot_priors()),
    "oracle priors"      = { pr <- oracle_priors(); if (is.null(pr)) NULL else c(cold, pr) },
    stop("unknown run: ", run)
  )
}

# update() rather than brm(): the pilot's compiled program is reused, because
# the participants and their RTs are data and the Stan code does not change
# with them. Only a new `prior` forces a recompile.
fit_run <- function(m_pilot, args, seed, dfx) {
  do.call(stats::update, c(
    list(m_pilot, newdata = dfx, init = args$init,
         chains = chains, cores = cores,
         iter = args$warmup + sample_iter, warmup = args$warmup,
         seed = seed, refresh = 0, silent = 2),
    args[intersect(names(args), c("inv_metric", "step_size", "prior", "stanvars"))]
  ))
}

# ---- Run ------------------------------------------------------------------------

stamp <- function(...) cat(sprintf("[%s] ", format(Sys.time(), "%H:%M:%S")), ..., "\n", sep = "")

for (family in families) {
  # target-major, then run, then seed
  todo <- do.call(rbind, lapply(names(targets), function(ds) {
    g <- expand.grid(seed = seeds, run = runs_for(ds), stringsAsFactors = FALSE)
    data.frame(dataset = ds, run = g$run, seed = g$seed, stringsAsFactors = FALSE)
  }))
  todo <- todo[!key(todo$dataset, family, todo$run, todo$seed) %in% done, , drop = FALSE]
  if (nrow(todo) == 0) {
    stamp("==== ", family, ": already complete, skipping ====")
    next
  }

  stamp("==== ", family, ": pilot fit on ", n_subset, " participants (compiles) ====")
  f <- make_formula(family, "mixed")
  # cogmod_priors() on the FULL first target, even for the pilot. Almost all
  # of what it sets is a constant, but brms's own default for the response
  # intercept is read off the response, so cogmod_priors(f, df_sub) and
  # cogmod_priors(f, df) differ by a few decimals - which would be a different
  # Stan program, so update() would recompile for every cell instead of
  # reusing the pilot's executable, and the cells would not all be stating the
  # same priors. This is the one prior specification every default-prior run
  # shares, on every target (update() keeps the fitted model's priors), and
  # the only thing the priors cells change.
  prior <- cogmod_priors(f, df)
  m_pilot <- tryCatch(
    suppressMessages(
      brm(f, data = df_sub, prior = prior, stanvars = cogmod_stanvars(f),
          init = cogmod_inits(f, df_sub), backend = "cmdstanr",
          chains = chains, cores = cores,
          iter = warmup_pilot + sample_iter, warmup = warmup_pilot,
          seed = 1, refresh = 0, silent = 2)
    ),
    error = function(e) { message("pilot failed: ", conditionMessage(e)); NULL }
  )
  if (is.null(m_pilot)) next
  if (!key("pilot_included", family, "pilot", 1) %in% done) {
    record("pilot_included", family, "pilot", 1, n_subset, nrow(df_sub), warmup_pilot,
           m_pilot, NA_real_)
  }

  # The pilot's own table, which each helper maps onto each target in turn.
  # Going through as.data.frame() is the same path a cluster job would take (a
  # few kilobytes of CSV instead of a brmsfit), and keeps the brmsfit out of
  # the helpers, so what is measured is what a user would actually pass. The
  # mapping onto each target is printed once, so the log says how many entries
  # carried over by level and how many are new participants.
  tab <- as.data.frame(cogmod_warmstart(m_pilot))
  write.csv(tab, file.path(out_dir, sprintf("warmstart_ablation_%s_pilot.csv", family)),
            row.names = FALSE)
  for (ds in unique(todo$dataset)) {
    cat("Mapped onto", ds, ": ")
    print(cogmod_warmstart(tab, f, targets[[ds]]))
  }

  for (i in seq_len(nrow(todo))) {
    ds <- todo$dataset[i]
    run <- todo$run[i]
    seed <- todo$seed[i]
    dfx <- targets[[ds]]
    stamp("---- ", family, " | ", ds, " | ", run, " | seed ", seed)
    args <- run_args(run, ds, family, f, tab, seed, dfx)
    if (is.null(args)) {
      message("no reference table for ", family, " on ", ds, " seed ", seed,
              " (", ref_file(ds, family, seed), "); skipping the `oracle priors` cell. ",
              "Run `reference` for this seed first.")
      next
    }
    timing <- system.time(
      fit <- tryCatch(fit_run(m_pilot, args, seed, dfx),
                      error = function(e) { message("failed: ", conditionMessage(e)); NULL })
    )
    if (is.null(fit)) next
    record(ds, family, run, seed, n_participants, nrow(dfx), args$warmup, fit,
           unname(timing["elapsed"]))

    # The reference is the `oracle priors` cell's source. Its own model on its
    # own data, so cogmod_warmstart() reproduces its adaptation and carries the
    # posterior median and SD of every parameter that has a prior.
    if (identical(run, "reference")) {
      write.csv(as.data.frame(cogmod_warmstart(fit)), ref_file(ds, family, seed),
                row.names = FALSE)
    }
    rm(fit); gc(verbose = FALSE)
  }

  rm(m_pilot); gc(verbose = FALSE)
}

# ---- Summary ------------------------------------------------------------------------

results <- do.call(rbind, results)
full <- results[results$run != "pilot", ]
if (!any(full$run == "base")) {
  cat("\nNothing to compare (need the `base` run).\n")
} else {
  agg <- aggregate(
    cbind(wall_s, warmup_s, sample_s, ess_bulk_min, ess_per_s, ess_per_grad,
          leapfrog_mean, stepsize_mean, divergent, rhat_max) ~ dataset + family + run,
    data = full, FUN = function(x) mean(x, na.rm = TRUE), na.action = na.pass
  )
  base <- agg[agg$run == "base", ]
  cmp <- merge(agg, base, by = c("dataset", "family"), suffixes = c("", "_base"))

  # Did the run change the answer, or only the route to it? Every run's
  # population-level posterior against the reference run's on the same
  # target, parameter by parameter and seed by seed: the ratio of the SDs says
  # whether the posterior is the same width, and the shift of the means in
  # units of the reference SD says whether it is in the same place.
  post <- do.call(rbind, posts)
  shifts <- data.frame(dataset = character(0), family = character(0), run = character(0),
                       sd_ratio = numeric(0), z_shift = numeric(0))
  if (!is.null(post) && any(post$run == "reference")) {
    ref <- post[post$run == "reference", c("dataset", "family", "seed", "variable", "mean", "sd")]
    pc <- merge(post, ref, by = c("dataset", "family", "seed", "variable"), suffixes = c("", "_ref"))
    pc <- pc[is.finite(pc$sd_ref) & pc$sd_ref > 0, , drop = FALSE]
    pc$sd_ratio <- pc$sd / pc$sd_ref
    pc$z_shift <- abs(pc$mean - pc$mean_ref) / pc$sd_ref
    shifts <- merge(
      aggregate(sd_ratio ~ dataset + family + run, pc, stats::median),
      aggregate(z_shift ~ dataset + family + run, pc, max),
      by = c("dataset", "family", "run")
    )
  }

  cell <- function(d) paste(d$dataset, d$family, sep = " | ")
  summary_tbl <- data.frame(
    dataset = cmp$dataset,
    family = cmp$family,
    run = cmp$run,
    n_pars = full$n_pars[match(cell(cmp), cell(full))],
    wall_s = round(cmp$wall_s),
    warmup_s = round(cmp$warmup_s),
    sample_s = round(cmp$sample_s),
    stepsize = round(cmp$stepsize_mean, 3),
    leapfrog_ratio = round(cmp$leapfrog_mean / cmp$leapfrog_mean_base, 2),
    ess_bulk_min = round(cmp$ess_bulk_min),
    ess_per_grad_ratio = round(cmp$ess_per_grad / cmp$ess_per_grad_base, 2),
    ess_per_s = round(cmp$ess_per_s, 2),
    ess_per_s_ratio = round(cmp$ess_per_s / cmp$ess_per_s_base, 2),
    divergent = round(cmp$divergent, 1),
    rhat_max = round(cmp$rhat_max, 3)
  )
  i <- match(paste(cell(summary_tbl), summary_tbl$run), paste(cell(shifts), shifts$run))
  summary_tbl$sd_ratio <- round(shifts$sd_ratio[i], 2)
  summary_tbl$z_shift <- round(shifts$z_shift[i], 2)
  run_order <- unique(c(runs, runs_transfer))
  summary_tbl <- summary_tbl[order(match(summary_tbl$dataset, datasets),
                                   match(summary_tbl$family, families),
                                   match(summary_tbl$run, run_order)), ]
  cat("\n==== relative to the same target's `base` (cold inits, short warmup) ====\n")
  cat("ess_per_s_ratio > 1: more effective draws per second, the number that matters.\n",
      "ess_per_grad_ratio > 1 / leapfrog_ratio < 1: fewer gradient evaluations per\n",
      "effective draw, i.e. the geometry is better scaled - that is the metric's job.\n",
      "`stepsize` well below the pilot's, in the `step_size` row, means the supplied\n",
      "step size did not survive warmup without the matching metric.\n",
      "rhat_max > 1.05: a chain got stuck; that cell's ess_per_s is a failure, not a speed.\n",
      "sd_ratio ~ 1 and z_shift small: the run sampled the same posterior as the\n",
      "reference, only faster or slower - which is what every row except the two\n",
      "priors rows must show for its ess_per_s_ratio to mean anything. `oracle priors`\n",
      "sits below 1 by the amount centring the priors on the reference's own posterior\n",
      "sharpened it, so its ess_per_s is an upper bound, not a like-for-like speedup.\n",
      "`all + pilot priors` on `shifted` with a large z_shift is a wrong prior moving the answer.\n\n",
      sep = "")
  print(summary_tbl, row.names = FALSE)
  write.csv(summary_tbl, file.path(out_dir, "warmstart_ablation_summary.csv"), row.names = FALSE)

  # The figure lives in its own file so that it can be redrawn from the CSV
  # without refitting anything: `Rscript benchmarks/warmstart_ablation_figure.R`.
  if (requireNamespace("ggplot2", quietly = TRUE)) {
    source("benchmarks/warmstart_ablation_figure.R")
    save_warmstart_ablation_figure(full, families = families, datasets = datasets)
  }
  stamp("Results written to ", normalizePath(out_dir))
}
