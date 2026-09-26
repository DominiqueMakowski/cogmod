# The design, as one table per experiment. fit_cell.R fits row
# SLURM_ARRAY_TASK_ID of it; summarise.R reads the rows back. Kept as a
# function of the experiment name so that the driver (run.sh) and the fitting
# script cannot disagree about how many cells there are.
#
# Experiment A - real lexical-decision RTs (speed_acc, rtdists), four families
# chosen for the four kinds of rule: log moments (LogNormal), the Wald
# inversion, the skewness-based ex-Gaussian, and root-finding on the CV
# (Weibull). Three data sets: 10 participants x 100 trials (as the laptop
# runs), the same shifted to 0.3 + 1.1 RT (the warmstart ablation's slow
# population, where the constants miss by a factor of four), and every trial
# of every participant. Two warmup lengths: 500, and 150, which is inside the
# cold-start transient AGENT.md 4.4.1 measured, so a start that shortens the
# transient has room to show.
#
# Experiment B - the cluster's own use case: the Muller-Lyer data of the
# Illusion Game with the production formula pattern (a 2-D tensor smooth of
# difficulty x illusion strength plus a participant intercept on every
# distributional parameter), at 120 and 480 participants, for the three
# families whose rules differ most. The choice column is ignored: these are
# RT-only families.
ablation_cells <- function(exp_id) {
  g <- switch(
    exp_id,
    A = expand.grid(
      seed = 1:8,
      scheme = c("constant", "data"),
      warmup = c(150L, 500L),
      family = c("lognormal", "invgaussian", "exgaussian", "weibull"),
      data = c("sa10", "sa10_shift", "sa_all"),
      stringsAsFactors = FALSE
    ),
    B = expand.grid(
      seed = 1:3,
      scheme = c("constant", "data"),
      warmup = c(300L, 1000L),
      family = c("lognormal", "invgaussian", "exgaussian"),
      data = c("igc120", "igc480"),
      stringsAsFactors = FALSE
    ),
    stop("unknown experiment '", exp_id, "'")
  )
  g$id <- seq_len(nrow(g))
  g$exp <- exp_id
  g[, c("exp", "id", "data", "family", "warmup", "scheme", "seed")]
}

# Chains and threads per cell. Four chains everywhere, so that Rhat means
# something; threads sized to the data.
ablation_layout <- function(exp_id) {
  switch(exp_id,
    A = list(chains = 4L, threads = 2L, samples = 500L),
    B = list(chains = 4L, threads = 4L, samples = 500L)
  )
}
