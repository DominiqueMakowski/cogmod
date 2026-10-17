# Why does stanli lose its per-gradient lead inside a 4-chain fit?
# ===============================================================
#
# multilevel.R's cluster run (2026-10-12) found stanli at 0.66-0.86x CmdStan
# per gradient, single-threaded, on five of its six cheap models, but at
# 0.82-1.55x per leapfrog step in the sampling phase of the same fits, four
# chains in parallel, and worse on the EPYC 7513 nodes than on the 9355s.
# Running inside the fit costs stanli 1.5-2.2x its single-threaded gradient,
# and CmdStan 1.3-1.7x. Two explanations: the four chains contend inside
# stanli's one process (they are threads of it: allocator, locks, shared
# state), or chains that are hungrier for memory contend for the node's
# caches and memory bandwidth. Five configurations, the same model, inits
# and seed:
#
#   cmdstan_1     CmdStan, one chain alone
#   cmdstan_4     CmdStan, four chains in parallel (four processes)
#   stanli_1      stanli, one chain alone
#   stanli_4      stanli, four chains in parallel (four threads of one process)
#   stanli_4proc  stanli, four R processes of one chain each, in parallel
#
# stanli_4 / stanli_1 against cmdstan_4 / cmdstan_1 is the extra cost of
# running stanli's chains together; stanli_4proc splits it: what it shares
# with stanli_4 is the hardware, what it does not is the process.
#
#   Rscript benchmarks/stanli/multilevel.R emit --model lnr
#   Rscript benchmarks/stanli/chains_probe.R --model lnr --seed 1 [--iter 300] [--tag NAME]
#
# The measure is the cost of a leapfrog step in the sampling phase, per
# chain: sampling time over sampling leapfrogs. Writes
# benchmarks/stanli/results/ml_chainsprobe_<model>_<seed>[_<tag>].csv.

source("benchmarks/gradient_programs.R")  # gp_args(), gp_load()
opt <- gp_args(list(out = "benchmarks/results/stanli_ml", model = "lnr", seed = 1L,
                    iter = 300L, tag = ""))
gp_load(".")
suppressPackageStartupMessages(library(stanli))
inp <- readRDS(file.path(opt$out, paste0(opt$model, "_inputs.rds")))
f <- inp$f; df <- inp$df
code <- as.character(brms::make_stancode(f, data = df, prior = inp$prior, stanvars = inp$sv))
sdat <- lapply(unclass(brms::make_standata(f, data = df, prior = inp$prior, stanvars = inp$sv)), identity)
set.seed(opt$seed)
inits <- cogmod_inits(f, df)
init_list <- lapply(1:4, function(i) inits(i))
# Off the seed again before sampling: see multilevel.R, `fit`.
set.seed(NULL)
cmod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(code, dir = opt$out), quiet = TRUE)
sm <- stanli_model(code = code, data = sdat)
init_u <- do.call(rbind, lapply(init_list, unconstrain, model = sm))
n <- opt$iter
cat(opt$model, "seed", opt$seed, "- stanli", as.character(utils::packageVersion("stanli")),
    "- CmdStan", as.character(cmdstanr::cmdstan_version()), "\n")

# Per chain: sampling seconds and sampling leapfrogs.
per_chain_cmdstan <- function(fc) {
  tt <- fc$time()$chains
  lf <- posterior::as_draws_df(fc$sampler_diagnostics())
  data.frame(sample_s = tt$sampling, leapfrog = as.numeric(tapply(lf$n_leapfrog__, lf$.chain, sum)))
}
per_chain_stanli <- function(fit) {
  sf <- as_stanfit(fit)
  sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
  data.frame(sample_s = rstan::get_elapsed_time(sf)[, "sample"],
             leapfrog = vapply(sp, function(x) sum(x[, "n_leapfrog__"]), numeric(1)))
}

# Four processes, each building the model and sampling one chain from its own
# init; the workers are started (and stanli loaded) before the clock matters.
cl <- parallel::makePSOCKcluster(4)
on.exit(parallel::stopCluster(cl), add = TRUE)
parallel::clusterExport(cl, c("code", "sdat", "init_u", "n", "opt"), envir = environment())
invisible(parallel::clusterEvalQ(cl, {
  suppressPackageStartupMessages(library(stanli))
  sm <- stanli_model(code = code, data = sdat)
  NULL
}))

configs <- list(
  cmdstan_1 = function() per_chain_cmdstan(cmod$sample(
    data = sdat, chains = 1, parallel_chains = 1, iter_warmup = n, iter_sampling = n,
    init = init_list[1], seed = opt$seed, refresh = 0, show_messages = FALSE)),
  cmdstan_4 = function() per_chain_cmdstan(cmod$sample(
    data = sdat, chains = 4, parallel_chains = 4, iter_warmup = n, iter_sampling = n,
    init = init_list, seed = opt$seed, refresh = 0, show_messages = FALSE)),
  stanli_1 = function() per_chain_stanli(sample_model(
    sm, chains = 1, parallel_chains = 1, seed = opt$seed, warmup = n, samples = n,
    init = init_u[1, , drop = FALSE], refresh = 0)),
  stanli_4 = function() per_chain_stanli(sample_model(
    sm, chains = 4, parallel_chains = 4, seed = opt$seed, warmup = n, samples = n,
    init = init_u, refresh = 0)),
  stanli_4proc = function() do.call(rbind, parallel::parLapply(cl, 1:4, function(i) {
    fit <- sample_model(sm, chains = 1, parallel_chains = 1, seed = 100 * opt$seed + i,
                        warmup = n, samples = n, init = init_u[i, , drop = FALSE], refresh = 0)
    sf <- as_stanfit(fit)
    sp <- rstan::get_sampler_params(sf, inc_warmup = FALSE)
    data.frame(sample_s = rstan::get_elapsed_time(sf)[, "sample"],
               leapfrog = sum(sp[[1]][, "n_leapfrog__"]))
  })))

rows <- list()
for (cfg in sample(names(configs))) {
  cat("config:", cfg, "\n")
  t0 <- Sys.time()
  x <- configs[[cfg]]()
  wall <- as.numeric(Sys.time() - t0, units = "secs")
  rows[[cfg]] <- data.frame(model = opt$model, seed = opt$seed, config = cfg,
                            position = length(rows) + 1, chain = seq_len(nrow(x)), x,
                            ms_per_leapfrog = 1000 * x$sample_s / x$leapfrog, wall_s = wall)
}
res <- do.call(rbind, rows)
key <- paste0(opt$model, "_", opt$seed, if (nzchar(opt$tag)) paste0("_", opt$tag))
dir.create("benchmarks/stanli/results", showWarnings = FALSE)
utils::write.csv(res, file.path("benchmarks/stanli/results", paste0("ml_chainsprobe_", key, ".csv")),
                 row.names = FALSE)
# Pooled over chains: sampling seconds over sampling leapfrogs.
s <- aggregate(cbind(sample_s, leapfrog) ~ config, res, sum)
s$ms_per_leapfrog <- 1000 * s$sample_s / s$leapfrog
base <- setNames(s$ms_per_leapfrog, s$config)
s$vs_alone <- s$ms_per_leapfrog / base[ifelse(grepl("^cmdstan", s$config), "cmdstan_1", "stanli_1")]
s$vs_cmdstan_same <- s$ms_per_leapfrog / base[sub("^cmdstan_4proc$", "cmdstan_4",
                                                  sub("^stanli", "cmdstan", s$config))]
print(s, digits = 3, row.names = FALSE)
