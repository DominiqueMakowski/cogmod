# Fit ONE chain of the LBA screen and write one CSV row plus its draws.
#
# Environment: SCREEN_EXPERIMENT (lba or rdm), SLURM_ARRAY_TASK_ID (the row of
# screen_cells()), SCREEN_OUT (where results go), SCREEN_LIB (this tree's
# cogmod), SLURM_CPUS_PER_TASK. For a local smoke test: SCREEN_LOAD_ALL=<this
# tree> for the 0.3.4 arms, and SCREEN_PARTICIPANTS / SCREEN_WARMUP /
# SCREEN_SAMPLES to shrink it.
#
# The arms differ in the formula and in which cogmod builds the priors and the
# starting values; everything else is analysis/server/fit_model.R of the
# Illusion Game project, so that lba_old is the production gam_lba on fewer
# participants. CmdStan is driven directly rather than through brm() because
# the warmup draws are part of what is measured, and brms discards them.

exp_id <- Sys.getenv("SCREEN_EXPERIMENT", unset = "lba")
task <- as.integer(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
out_dir <- Sys.getenv("SCREEN_OUT", unset = "results")
cpus <- as.integer(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "4"))
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source("cells.R")
cell <- screen_cells(exp_id)[task, ]
env_int <- function(name, default) {
  v <- Sys.getenv(name)
  if (nzchar(v)) as.integer(v) else default
}
n_part <- env_int("SCREEN_PARTICIPANTS", screen_layout$participants)
n_warm <- env_int("SCREEN_WARMUP", screen_layout$warmup)
n_samp <- env_int("SCREEN_SAMPLES", screen_layout$samples)
stem <- sprintf("%s_%02d_%s_chain%d", exp_id, task, cell$arm, cell$chain)
out_file <- file.path(out_dir, paste0(stem, ".csv"))
if (file.exists(out_file)) {
  cat("cell", task, "already done:", out_file, "\n")
  quit(status = 0)
}
cat("cell:", paste(names(cell), unlist(cell), sep = "=", collapse = " "), "\n")


# Which cogmod -------------------------------------------------------------

# The 0.3.3 arm must not see this tree's cogmod: drop its library from the
# search path so library(cogmod) falls through to the production one. The
# check afterwards is what keeps a mixed-up library from passing silently.
screen_lib <- Sys.getenv("SCREEN_LIB")
if (cell$defaults == "0.3.3" && nzchar(screen_lib)) {
  .libPaths(setdiff(.libPaths(), normalizePath(screen_lib, mustWork = FALSE)))
}
suppressPackageStartupMessages({
  library(brms)
  library(cmdstanr)
  library(posterior)
  if (cell$defaults != "0.3.3" && nzchar(Sys.getenv("SCREEN_LOAD_ALL"))) {
    pkgload::load_all(Sys.getenv("SCREEN_LOAD_ALL"), quiet = TRUE)
  } else {
    library(cogmod)
  }
})
ns <- asNamespace("cogmod")
new_defaults <- exists(".zs_scales", envir = ns)
cat("cogmod", as.character(packageVersion("cogmod")), "from",
    tryCatch(find.package("cogmod"), error = function(e) "load_all"),
    "| 0.3.4 defaults:", new_defaults, "\n")
cat("brms", as.character(packageVersion("brms")), "| cmdstanr",
    as.character(packageVersion("cmdstanr")), "| CmdStan",
    as.character(cmdstan_version()), "\n")
if (new_defaults != (cell$defaults == "0.3.4")) {
  stop("arm ", cell$arm, " wants the ", cell$defaults, " defaults but got a ",
       "cogmod ", if (new_defaults) "with" else "without", " them")
}


# Data, as analysis/server/fit_model.R prepares it ------------------------

base <- "https://raw.githubusercontent.com/RealityBending/IllusionGameComputational/refs/heads/main/data/"
df <- do.call(rbind, lapply(1:3, function(i) {
  read.csv(paste0(base, "illusion_part", i, ".csv"))
}))
df$Illusion_Difference <- abs(df$Illusion_Difference)
# datawizard::normalize() within illusion, without the dependency
nrm <- function(x) (x - min(x)) / (max(x) - min(x))
df$Illusion_DifferenceZ <- NA_real_
df$Illusion_StrengthZ <- NA_real_
for (ty in unique(df$Illusion_Type)) {
  r <- df$Illusion_Type == ty
  df$Illusion_DifferenceZ[r] <- 2 * nrm(df$Illusion_Difference[r]) - 1
  df$Illusion_StrengthZ[r] <- sign(df$Illusion_Strength[r]) *
    nrm(abs(df$Illusion_Strength[r]))
}
keep <- unique(df$Participant)[seq_len(n_part)]
data <- df[df$Participant %in% keep & df$Illusion_Type == "MullerLyer", ]
rownames(data) <- NULL
cat("rows:", nrow(data), " participants:", length(unique(data$Participant)),
    " error rate:", round(mean(data$Error), 3), "\n")


# Formula, as analysis/server/models.R writes it ---------------------------

t2f <- function(lhs) {
  stats::as.formula(paste0(
    lhs, " ~ t2(Illusion_DifferenceZ, Illusion_StrengthZ, ",
    "k = c(5, 5), bs = c('cr', 'cr')) + (1 | Participant)"
  ), env = globalenv())
}
p_only <- function(lhs) {
  stats::as.formula(paste(lhs, "~ 1 + (1 | Participant)"), env = globalenv())
}
rt <- t2f("RT | dec(Error)")
out <- poutlier ~ 1 + (1 | Participant)
f <- switch(cell$model,
  gam_lba = bf(rt, t2f("driftone"), t2f("sigmaone"), t2f("sigmabias"),
               t2f("boundary"), t2f("ndt"), sigmazero = 1, out,
               family = cogmod_lba2()),
  gam_lba_a = bf(rt, t2f("driftone"), p_only("sigmaone"), t2f("sigmabias"),
                 t2f("boundary"), t2f("ndt"), sigmazero = 1, out,
                 family = cogmod_lba2()),
  gam_lba_b = bf(rt, t2f("driftone"), t2f("sigmaone"), p_only("sigmabias"),
                 t2f("boundary"), t2f("ndt"), sigmazero = 1, out,
                 family = cogmod_lba2()),
  gam_rdm5 = bf(rt, t2f("driftone"), t2f("sigmabias"), t2f("boundary"),
                t2f("ndt"), out, family = cogmod_rdm()),
  gam_rdm5_b = bf(rt, t2f("driftone"), p_only("sigmabias"), t2f("boundary"),
                  t2f("ndt"), out, family = cogmod_rdm()),
  stop("unknown model ", cell$model)
)

# cogmod's priors plus normal(0, 1) on the slopes brms leaves flat, as
# fit_model.R does. A dpar with no population-level slope (`~ 1 + (1 | P)`)
# has no blanket `b` row, and a prior on one is an error.
priors <- cogmod_priors(f, data)
for (par in c("", setdiff(unique(priors$dpar), c("poutlier", "")))) {
  blanket <- priors$class == "b" & priors$dpar == par &
    !nzchar(priors$coef) & !nzchar(priors$group)
  if (!any(blanket) || any(blanket & nzchar(priors$prior))) next
  priors <- c(priors, brms::prior_string("normal(0, 1)", class = "b", dpar = par),
              replace = TRUE)
}


# Fit ----------------------------------------------------------------------

threads <- max(1L, min(screen_layout$threads, cpus))
sv <- cogmod_stanvars(f)
code <- make_stancode(f, data = data, prior = priors, stanvars = sv,
                      threads = threading(threads), backend = "cmdstanr")
sdata <- as.list(make_standata(f, data = data, prior = priors, stanvars = sv,
                               threads = threading(threads)))
# The production cpp_options, so the precompiled header the production runs
# built is the one used here.
t_compile <- system.time(
  mod <- cmdstan_model(
    write_stan_file(code, dir = tempdir()),
    cpp_options = list(stan_threads = TRUE, STAN_CPP_OPTIMS = TRUE,
                       STAN_NO_RANGE_CHECKS = TRUE),
    stanc_options = list("O1")
  )
)[["elapsed"]]
cat(sprintf("compiled in %.0f s\n", t_compile))

# Cold start, as in production: cogmod_inits() with its default jitter. The
# init draw is seeded by the chain number, so chain k of every arm starts
# from the same random stream.
set.seed(cell$chain)
init <- cogmod_inits(f, data)(1)
row <- data.frame(cell, participants = n_part, rows = nrow(data),
                  warmup = n_warm, samples = n_samp, threads = threads,
                  node = Sys.info()[["nodename"]], status = "ok",
                  stringsAsFactors = FALSE)
t_fit <- Sys.time()
fit <- tryCatch(
  mod$sample(
    data = sdata, chains = 1, threads_per_chain = threads,
    iter_warmup = n_warm, iter_sampling = n_samp, seed = cell$chain,
    init = list(init), save_warmup = TRUE, refresh = 0,
    show_messages = FALSE, output_dir = tempdir()
  ),
  error = function(e) e
)
row$wall_h <- round(as.numeric(Sys.time() - t_fit, units = "hours"), 2)
if (inherits(fit, "error")) {
  row$status <- paste("error:", conditionMessage(fit))
  write.csv(row, out_file, row.names = FALSE)
  cat("FAILED:", conditionMessage(fit), "\n")
  quit(status = 0)
}


# Per-chain metrics ----------------------------------------------------------

diag <- fit$sampler_diagnostics(inc_warmup = TRUE)
lf <- as.numeric(diag[, 1, "n_leapfrog__"])
td <- as.numeric(diag[, 1, "treedepth__"])
dv <- as.numeric(diag[, 1, "divergent__"])
wi <- seq_len(n_warm)
si <- n_warm + seq_len(n_samp)
row$warmup_leapfrog <- sum(lf[wi])
row$warmup_maxtree_frac <- round(mean(td[wi] >= 10), 3)
row$sample_leapfrog_mean <- round(mean(lf[si]), 1)
row$sample_maxtree_frac <- round(mean(td[si] >= 10), 3)
row$divergences <- sum(dv[si])
row$step_size <- signif(unlist(fit$metadata()$step_size_adaptation)[1], 3)
lp <- as.numeric(fit$draws("lp__", inc_warmup = TRUE))
row$lp_mean <- round(mean(lp[si]), 1)
row$lp_sd <- round(stats::sd(lp[si]), 1)

# Population-level parameters for the pooled Rhat in summarise.R:
# intercepts, slopes, the smooths' linear parts, and every group and smooth
# SD. Participant-level draws stay out; they are many and say little here.
vars <- grep("^(b|Intercept|bs|sd|sds)(_[A-Za-z0-9_]*)?$",
             fit$metadata()$stan_variables, value = TRUE)
draws <- fit$draws(variables = c("lp__", vars))
saveRDS(list(cell = row, draws = draws, lp_warmup = lp[wi]),
        file.path(out_dir, paste0(stem, ".rds")))
write.csv(row, out_file, row.names = FALSE)
cat("REPORT", paste(names(row), unlist(row), sep = "=", collapse = " "), "\n")
