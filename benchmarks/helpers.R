# Shared by the benchmark scripts in this directory: the dataset and the model
# formulas. Sourced from the repository root by metric_dense_e.R and
# warm_start.R; not meant to be run on its own.

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  library(cogmod)
})

# ---- Data -----------------------------------------------------------------------
# Experiment 1 of Wagenmakers et al. (2008), as in the decision-making article:
# lexical decision under speed vs accuracy instructions, choice coded as
# correct (0) vs error (1). `n_participants` are drawn at random (fixed seed)
# and `n_trials` kept per participant, balanced across the two conditions.

bench_data <- function(n_participants, n_trials, seed = 2026) {
  data(speed_acc, package = "rtdists")
  df_all <- data.frame(
    Participant = as.integer(as.character(speed_acc$id)),
    Condition = unname(c(accuracy = "Accuracy", speed = "Speed")[as.character(speed_acc$condition)]),
    RT = speed_acc$rt,
    Error = as.integer(as.character(speed_acc$response) != as.character(speed_acc$stim_cat))
  )
  df_all <- df_all[df_all$RT <= 2, ]

  set.seed(seed)
  ids <- sort(sample(unique(df_all$Participant), n_participants))
  df <- do.call(rbind, lapply(ids, function(id) {
    d <- df_all[df_all$Participant == id, ]
    do.call(rbind, lapply(split(d, d$Condition), function(dc) {
      dc[sample(nrow(dc), min(nrow(dc), n_trials %/% 2)), ]
    }))
  }))
  df$Participant <- factor(df$Participant)
  rownames(df) <- NULL
  cat(sprintf("Data: %d trials, %d participants, %.1f%% errors\n",
              nrow(df), n_participants, 100 * mean(df$Error)))
  df
}

# ---- Models -------------------------------------------------------------------
# Which dpars get the Condition effect, which are intercept-only, which carry a
# participant random intercept under the "mixed" structure, and which are fixed
# constants. Condition sits on the drifts and (where the family has one) the
# boundary, the classic speed-accuracy trade-off account.

specs <- list(
  ddm = list(
    cond = c("mu", "boundary"),
    int = c("bias", "ndt", "poutlier"),
    re = c("mu", "boundary", "ndt"),
    fixed = list(sigmadrift = 0, sigmabias = 0, sigmandt = 0)
  ),
  lba2 = list(
    cond = c("mu", "driftone", "boundary"),
    int = c("sigmaone", "sigmabias", "ndt", "poutlier"),
    re = c("mu", "boundary", "ndt"),
    fixed = list(sigmazero = 1)          # LBA evidence scale is arbitrary; see ?cogmod_lba2
  ),
  lnr = list(
    cond = c("mu", "nuone"),
    int = c("sigmazero", "sigmaone", "ndt", "poutlier"),
    re = c("mu", "ndt"),
    fixed = list()
  ),
  rdm = list(
    cond = c("mu", "driftone", "boundary"),
    int = c("sigmabias", "ndt", "poutlier"),
    re = c("mu", "boundary", "ndt"),
    fixed = list()
  )
)

make_formula <- function(family, structure = c("fixed", "mixed")) {
  structure <- match.arg(structure)
  spec <- specs[[family]]
  rhs <- function(dpar) {
    fe <- if (dpar %in% spec$cond) "Condition" else "1"
    re <- if (structure == "mixed" && dpar %in% spec$re) " + (1 | Participant)" else ""
    paste0(fe, re)
  }
  lhs <- function(dpar) if (dpar == "mu") "RT | dec(Error)" else dpar
  dpars <- c(spec$cond, spec$int)
  forms <- lapply(dpars, function(d) as.formula(paste(lhs(d), "~", rhs(d))))
  fam <- get(paste0("cogmod_", family), envir = asNamespace("cogmod"))()
  do.call(brms::bf, c(forms, spec$fixed, list(family = fam)))
}

# ---- Diagnostics ------------------------------------------------------------------

# Smallest bulk ESS and largest Rhat over the parameters of a brmsfit (or of a
# posterior draws object), ignoring constants and the internal `z_` scores.
ess_rhat <- function(x) {
  d <- if (inherits(x, "brmsfit")) {
    vars <- variables(x)
    vars <- vars[!vars %in% c("lp__", "lprior") & !startsWith(vars, "z_")]
    as_draws_array(x, variable = vars)
  } else x
  s <- summarise_draws(d, "ess_bulk", "rhat")
  c(ess_bulk_min = min(s$ess_bulk, na.rm = TRUE), rhat_max = max(s$rhat, na.rm = TRUE))
}
