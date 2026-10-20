# One HMC run of the RT article's ex-Wald model (README.md), warmup saved.
#
#   Rscript benchmarks/exwald_modes/fit.R <prior> <inits> <name> [warmup] [sampling]
#
#   prior     0.3.4 | 0.3.5 | P1 (common.R)
#   inits     seed:<s>[,<s>...]  cogmod_inits() exactly as a 4-chain brm() call with
#                                set.seed(s) gets them: cmdstanr calls init(1) once
#                                before the per-chain calls, so chain k starts from
#                                the (k + 1)th draw. Four chains per seed.
#             modes:<file.rds>   one chain at each mode in modes.R's output (jitter
#                                0.02 on every unconstrained coordinate)
#   name      run directory under EXW_SCRATCH (default: tempdir())
#
# CmdStan seed = 1000 + the first init seed (or 11 for modes). Summarised by
# collect.R.
args <- commandArgs(trailingOnly = TRUE)
source("benchmarks/exwald_modes/common.R")
prior <- args[1]
inits_arg <- args[2]
name <- args[3]
W <- if (length(args) >= 4) as.integer(args[4]) else 250L
S <- if (length(args) >= 5) as.integer(args[5]) else 250L
scratch <- Sys.getenv("EXW_SCRATCH", tempdir())
rundir <- file.path(scratch, name)
dir.create(rundir, recursive = TRUE, showWarnings = FALSE)

prog <- make_program(prior)
exedir <- file.path(scratch, paste0("exe_", gsub("[^0-9A-Za-z]", "", prior)))
dir.create(exedir, showWarnings = FALSE)
mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(prog$code, dir = exedir), quiet = TRUE)

kind <- sub(":.*", "", inits_arg)
val <- sub("^[a-z]+:", "", inits_arg)
inits <- list()
origin <- character(0)
if (kind == "seed") {
  seeds <- as.integer(strsplit(val, ",")[[1]])
  for (s in seeds) {
    set.seed(s)
    ini <- cogmod_inits(FORMULA, D)
    invisible(ini(1))
    for (ch in 1:4) {
      inits[[length(inits) + 1]] <- ini(ch)
      origin <- c(origin, sprintf("seed %d chain %d", s, ch))
    }
  }
  sseed <- 1000 + seeds[1]
} else {
  md <- readRDS(val)
  set.seed(42)
  for (k in seq_along(md$par)) {
    inits[[k]] <- as_init(md$par[[k]] + stats::rnorm(9, 0, 0.02))
    origin <- c(origin, paste("mode", md$label[k]))
  }
  sseed <- 11
}
fit <- mod$sample(
  data = prog$sdata, seed = sseed, chains = length(inits), parallel_chains = length(inits),
  iter_warmup = W, iter_sampling = S, init = inits, refresh = 0, show_messages = FALSE,
  show_exceptions = FALSE, save_warmup = TRUE, metric = "diag_e", output_dir = rundir
)
saveRDS(list(
  prior = prior, origin = origin, warmup = W, inits = fit$init(),
  draws = as.data.frame(fit$draws(inc_warmup = TRUE, format = "draws_df")),
  sdg = as.data.frame(fit$sampler_diagnostics(inc_warmup = TRUE, format = "draws_df")),
  time = fit$time()
), file.path(rundir, "draws.rds"))
cat("done:", rundir, "\n")
