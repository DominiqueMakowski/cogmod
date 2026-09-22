# Fit ONE cell of the inits ablation and write one CSV row.
#
# Environment: ABL_EXPERIMENT (A or B), SLURM_ARRAY_TASK_ID (the row of
# ablation_cells()), ABL_OUT (where the row goes), ABL_DATA (where
# speed_acc.csv is), SLURM_CPUS_PER_TASK.
#
# The two schemes differ ONLY in the `init` list handed to CmdStan. "data" is
# cogmod_inits() as built in this tree, with the data-aware layer; "constant"
# is cogmod_inits() as it was before it - the registry constants, with only
# ndt read off the response - rebuilt here from the same internals so that a
# single installed cogmod serves both arms.
#
# brms is used for what it is good at - the Stan program, the data list, the
# priors - and CmdStan is driven directly for the fit, because the metrics
# that matter here live in the WARMUP draws (lp__ and the sampler diagnostics
# iteration by iteration), which brms discards.

suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  # ABL_LOAD_ALL=<package dir> runs the script against a source tree (a local
  # smoke test); on the cluster the tarball is installed and this is unset.
  if (nzchar(Sys.getenv("ABL_LOAD_ALL"))) {
    pkgload::load_all(Sys.getenv("ABL_LOAD_ALL"), quiet = TRUE)
  } else {
    library(cogmod)
  }
})

exp_id <- Sys.getenv("ABL_EXPERIMENT", unset = "A")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
out_dir <- Sys.getenv("ABL_OUT", unset = "results")
data_dir <- Sys.getenv("ABL_DATA", unset = "data")
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "8"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("cells.R")
cells <- ablation_cells(exp_id)
layout <- ablation_layout(exp_id)
cell <- cells[task, ]
out_file <- file.path(out_dir, sprintf("%s_%04d.csv", exp_id, task))
if (file.exists(out_file)) {
  cat("cell", task, "already done:", out_file, "\n")
  quit(status = 0)
}
cat("cell:", paste(names(cell), unlist(cell), sep = "=", collapse = " "), "\n")
cat("cogmod", as.character(packageVersion("cogmod")), "from", find.package("cogmod"),
    "| data-aware layer:", exists(".data_start", envir = asNamespace("cogmod")), "\n")
cat("brms", as.character(packageVersion("brms")), "| cmdstanr",
    as.character(packageVersion("cmdstanr")), "| CmdStan", as.character(cmdstan_version()), "\n")
if (!exists(".data_start", envir = asNamespace("cogmod"))) {
  stop("the installed cogmod has no data-aware layer; the two arms would be identical")
}


# Data ----------------------------------------------------------------------

speed_acc_data <- function(which) {
  sa <- read.csv(file.path(data_dir, "speed_acc.csv"))
  df <- data.frame(Participant = sa$id, Condition = sa$condition, RT = sa$rt)
  df <- df[is.finite(df$RT) & df$RT <= 2, ]
  if (which %in% c("sa10", "sa10_shift")) {
    set.seed(2026)
    ids <- sort(sample(unique(df$Participant), 10))
    df <- do.call(rbind, lapply(ids, function(id) {
      d <- df[df$Participant == id, ]
      do.call(rbind, lapply(split(d, d$Condition), function(dc) {
        dc[sample(nrow(dc), 50), ]
      }))
    }))
  }
  if (which == "sa10_shift") df$RT <- 0.3 + 1.1 * df$RT
  df$Participant <- factor(df$Participant)
  df$Condition <- factor(df$Condition)
  rownames(df) <- NULL
  df
}

# The Illusion Game data, prepared as analysis/server/fit_model.R prepares it
# (Muller-Lyer only; difficulty rescaled to [-1, 1], strength to [-1, 1] with
# its sign), without the datawizard dependency.
igc_data <- function(n_participants) {
  base <- "https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/"
  df <- do.call(rbind, lapply(1:3, function(i) {
    read.csv(paste0(base, "illusion_part", i, ".csv"))
  }))
  df <- df[df$Illusion_Type == "MullerLyer", ]
  df$Illusion_Difference <- abs(df$Illusion_Difference)
  nrm <- function(x) (x - min(x)) / (max(x) - min(x))
  df$Illusion_DifferenceZ <- 2 * nrm(df$Illusion_Difference) - 1
  df$Illusion_StrengthZ <- sign(df$Illusion_Strength) * nrm(abs(df$Illusion_Strength))
  keep <- unique(df$Participant)[seq_len(n_participants)]
  df <- df[df$Participant %in% keep, ]
  df$Participant <- factor(df$Participant)
  rownames(df) <- NULL
  df
}

data <- switch(cell$data,
  sa10 = , sa10_shift = , sa_all = speed_acc_data(cell$data),
  igc120 = igc_data(120L),
  igc480 = igc_data(480L)
)
cat("rows:", nrow(data), " participants:", nlevels(data$Participant),
    " RT median:", round(median(data$RT), 3), "\n")


# Formula -------------------------------------------------------------------

t2f <- function(lhs) {
  stats::as.formula(paste0(
    lhs, " ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, ",
    "k = c(5, 5), bs = c('cr', 'cr')) + (1 | Participant)"
  ), env = globalenv())
}

f <- if (exp_id == "A") {
  switch(cell$family,
    lognormal = bf(RT ~ Condition + (1 | Participant), sigma ~ 1, ndt ~ 1,
                   sigmabias = 0, family = cogmod_lognormal()),
    invgaussian = bf(RT ~ Condition + (1 | Participant), boundary ~ Condition,
                     ndt ~ 1, sigmadrift = 0, sigmandt = 0,
                     family = cogmod_invgaussian()),
    exgaussian = bf(RT ~ Condition + (1 | Participant), sigma ~ 1,
                    tau ~ Condition, family = cogmod_exgaussian()),
    # shape on mu, scale on sigma: the participant sits on the scale
    weibull = bf(RT ~ 1, sigma ~ Condition + (1 | Participant), ndt ~ 1,
                 family = cogmod_weibull())
  )
} else {
  switch(cell$family,
    lognormal = bf(t2f("RT"), t2f("sigma"), sigmabias = 0, t2f("ndt"),
                   poutlier ~ 1 + (1 | Participant), family = cogmod_lognormal()),
    invgaussian = bf(t2f("RT"), t2f("boundary"), sigmadrift = 0, sigmandt = 0,
                     t2f("ndt"), poutlier ~ 1 + (1 | Participant),
                     family = cogmod_invgaussian()),
    exgaussian = bf(t2f("RT"), t2f("sigma"), t2f("tau"),
                    family = cogmod_exgaussian())
  )
}

# cogmod's priors, plus normal(0, 1) on the slopes brms leaves flat - the
# rule analysis/server/fit_model.R uses, so that B is the production model.
priors <- cogmod_priors(f, data)
for (par in c("", setdiff(unique(priors$dpar), c("poutlier", "")))) {
  blanket <- priors$class == "b" & priors$dpar == par &
    !nzchar(priors$coef) & !nzchar(priors$group)
  # A dpar with no population-level coefficients (sigma ~ 1) has no `b` row
  # at all, and a prior on one is an error; skip it as well as the ones
  # cogmod has already set.
  if (!any(blanket) || any(blanket & nzchar(priors$prior))) next
  priors <- c(priors, brms::prior_string("normal(0, 1)", class = "b", dpar = par),
              replace = TRUE)
}


# Starting values -----------------------------------------------------------

# cogmod_inits() before the data-aware layer: the constants, ndt from the
# response. Same plan machinery, same jitter, same tiering.
constant_inits <- function(formula, data) {
  ns <- asNamespace("cogmod")
  family <- ns$.cogmod_family(formula)
  targets <- ns$.init_targets(family)
  links <- ns$.family_links(family)
  code <- suppressWarnings(make_stancode(formula, data = data, family = family))
  sdata <- suppressWarnings(make_standata(formula, data = data, family = family))
  if (!is.null(targets$ndt)) targets$ndt <- ns$.ndt_start(sdata$Y, targets$ndt)
  plan <- ns$.init_plan(ns$.stan_param_decls(code), as.list(sdata), targets, links)
  out <- ns$.init_fun(plan, NULL)
  attr(out, "targets") <- targets
  out
}

init_fun <- if (cell$scheme == "data") cogmod_inits(f, data) else constant_inits(f, data)
targets <- attr(init_fun, "targets")
from_data <- attr(targets, "from_data")
cat("targets:", paste(names(targets), signif(unlist(targets), 4), sep = "=", collapse = " "),
    "| from data:", if (length(from_data)) paste(from_data, collapse = ",") else "-", "\n")


# Fit -----------------------------------------------------------------------

chains <- layout$chains
threads <- max(1L, min(layout$threads, cpus %/% chains))
sv <- cogmod_stanvars(f)
code <- make_stancode(f, data = data, prior = priors, stanvars = sv,
                      threads = threading(threads), backend = "cmdstanr")
sdata <- make_standata(f, data = data, prior = priors, stanvars = sv,
                       threads = threading(threads))
sdata <- as.list(sdata)

# Same cpp_options as analysis/server/fit_model.R and precompile.R, so the
# precompiled header the production runs built is the one used here and no
# array task ever tries to build one (AGENT.md 3.3).
t_compile <- system.time(
  mod <- cmdstan_model(
    write_stan_file(code, dir = tempdir()),
    cpp_options = list(stan_threads = TRUE, STAN_CPP_OPTIMS = TRUE,
                       STAN_NO_RANGE_CHECKS = TRUE),
    stanc_options = list("O1")
  )
)[["elapsed"]]
cat(sprintf("compiled in %.0f s\n", t_compile))

set.seed(cell$seed)
init_list <- lapply(seq_len(chains), function(i) init_fun(i))

row <- data.frame(
  cell, n_rows = nrow(data), n_participants = nlevels(data$Participant),
  chains = chains, threads = threads, node = Sys.info()[["nodename"]],
  from_data = paste(from_data, collapse = ","), status = "ok",
  stringsAsFactors = FALSE
)

t_fit <- Sys.time()
fit <- tryCatch(
  mod$sample(
    data = sdata, chains = chains, parallel_chains = chains,
    threads_per_chain = threads, iter_warmup = cell$warmup,
    iter_sampling = layout$samples, seed = cell$seed, init = init_list,
    save_warmup = TRUE, refresh = 0, show_messages = FALSE,
    output_dir = tempdir()
  ),
  error = function(e) e
)
row$wall_s <- round(as.numeric(Sys.time() - t_fit, units = "secs"), 1)

if (inherits(fit, "error")) {
  row$status <- paste("error:", conditionMessage(fit))
  write.csv(row, out_file, row.names = FALSE)
  cat("FAILED:", conditionMessage(fit), "\n")
  quit(status = 0)
}


# Metrics -------------------------------------------------------------------

n_warm <- cell$warmup
n_samp <- layout$samples
tm <- fit$time()$chains
row$n_chains_ok <- nrow(tm)
row$warmup_s <- round(sum(tm$warmup), 1)
row$sample_s <- round(sum(tm$sampling), 1)

diag <- fit$sampler_diagnostics(inc_warmup = TRUE) # iter x chain x variable
lf <- diag[, , "n_leapfrog__", drop = TRUE]
td <- diag[, , "treedepth__", drop = TRUE]
dv <- diag[, , "divergent__", drop = TRUE]
dim(lf) <- dim(td) <- dim(dv) <- c(n_warm + n_samp, row$n_chains_ok)
wi <- seq_len(n_warm)
si <- n_warm + seq_len(n_samp)
row$warmup_leapfrog <- sum(lf[wi, ])
row$warmup_maxtree_frac <- round(mean(td[wi, ] >= 10), 3)
row$warmup_divergences <- sum(dv[wi, ])
row$sample_leapfrog_mean <- round(mean(lf[si, ]), 1)
row$sample_treedepth_mean <- round(mean(td[si, ]), 2)
row$divergences <- sum(dv[si, ])
row$step_size <- signif(mean(unlist(fit$metadata()$step_size_adaptation)), 3)

# The cold-start transient: for each chain, the first warmup iteration at
# which lp__ reaches the 5th percentile of that chain's own post-warmup lp__.
# The worst chain is what a fit waits for.
lp <- fit$draws("lp__", inc_warmup = TRUE)
lp <- matrix(as.numeric(lp), nrow = n_warm + n_samp)
trans <- apply(lp, 2, function(x) {
  q <- quantile(x[si], 0.05, names = FALSE)
  hit <- which(x[wi] >= q)
  if (length(hit)) hit[1] else NA_integer_
})
row$transient_max <- if (all(is.na(trans))) NA else max(trans, na.rm = TRUE)
row$transient_mean <- round(mean(trans, na.rm = TRUE), 1)
# Leapfrog steps spent before every chain had arrived: the cost of the start.
row$transient_leapfrog <- if (is.na(row$transient_max)) NA else sum(lf[seq_len(row$transient_max), ])

# Population-level parameters only: intercepts, slopes, smooth coefficients
# and every group / smooth SD. Participant-level draws are excluded because a
# single participant's ndt / poutlier trade-off dominates the global maximum
# (AGENT.md 4.4.1) and would hide what the start changes.
vars <- grep("^(b|Intercept|bs|sd|sds)(_[A-Za-z0-9_]*)?(\\[[0-9]+\\])?$|^(poutlier|ndt|sigma|tau|boundary|shape|sigmabias)$",
             fit$metadata()$stan_variables, value = TRUE)
vars <- setdiff(vars, c("b", "bs")[!c("b", "bs") %in% fit$metadata()$stan_variables])
sm <- summarise_draws(fit$draws(variables = vars), "rhat", "ess_bulk", "ess_tail")
row$n_pop_params <- nrow(sm)
row$max_rhat <- round(max(sm$rhat, na.rm = TRUE), 4)
row$median_rhat <- round(median(sm$rhat, na.rm = TRUE), 4)
row$min_ess_bulk <- round(min(sm$ess_bulk, na.rm = TRUE))
row$median_ess_bulk <- round(median(sm$ess_bulk, na.rm = TRUE))
row$min_ess_tail <- round(min(sm$ess_tail, na.rm = TRUE))
row$min_ess_per_s <- round(row$min_ess_bulk / row$sample_s, 2)
row$median_ess_per_s <- round(row$median_ess_bulk / row$sample_s, 2)
# ESS per second of the WHOLE chain, warmup included: what a user pays.
row$min_ess_per_total_s <- round(row$min_ess_bulk / (row$warmup_s + row$sample_s), 3)

write.csv(row, out_file, row.names = FALSE)
cat("REPORT", paste(names(row), unlist(row), sep = "=", collapse = " "), "\n")
