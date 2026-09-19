# Stage A: end to end with a warm start. The posterior mode is read off an
# optimiser CSV (--mode); both trees then run the same short warmup from that
# state, under the same hand-set diagonal metric, and the same number of
# sampling iterations. See README.md.
#
#   Rscript ess_warm.R --base <dir> --pr <dir> --out <dir> --mode <csv>
#           [--warmup 100 --sampling 200 --chains 4 --order base,pr --treedepth 7]
args <- commandArgs(trailingOnly = TRUE)
opt <- list(base = "", pr = "", out = ".", warmup = 150L, sampling = 400L,
            chains = 4L, order = "base,pr", seed = 3L, treedepth = 8L, mode = "")
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[i + 1]) else args[i + 1]
}
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

dirs <- c(base = opt$base, pr = opt$pr)
mods <- lapply(dirs, function(d) cmdstan_model(file.path(d, "cogmod_ddm.stan"), dir = file.path(d, "exe")))
data <- file.path(dirs[["base"]], "cogmod_ddm.data.json")
init <- file.path(dirs[["base"]], "cogmod_ddm.init.json")
pars <- c("Intercept", paste0("Intercept_", c("boundary", "bias", "sigmadrift",
                                             "sigmabias", "sigmandt", "ndt", "poutlier")))

# ---- warm start ----------------------------------------------------------------
# The mode from the optimiser run of 17:34 (347 s, lp = -39.35), read off its
# CSV; a Laplace metric was tried and abandoned - its draws land where the
# 7-parameter density costs seconds per evaluation. The starting metric is a
# hand-set diagonal in the posterior's rough units, identical for both trees;
# window adaptation replaces it after the first window either way.
mode_csv <- opt$mode
m <- read.csv(mode_csv, comment.char = "#")
mu <- as.numeric(m[1, pars])
inv_metric <- c(0.15, 0.05, 0.1, 0.3, 1, 0.5, 0.1, 1)^2
cat("mode:      ", paste(pars, signif(mu, 4), sep = "=", collapse = " "), "
")
cat("inv_metric:", signif(inv_metric, 3), "
")
# each chain starts at the mode nudged by a quarter of the metric's sd
set.seed(opt$seed)
inits <- lapply(seq_len(opt$chains), function(ch) {
  as.list(setNames(mu + 0.25 * sqrt(inv_metric) * rnorm(length(mu)), pars))
})
saveRDS(list(mode = mu, inv_metric = inv_metric, inits = inits), file.path(opt$out, "warm.rds"))

# ---- runs --------------------------------------------------------------------
order <- strsplit(opt$order, ",")[[1]]
rows <- list()
for (k in seq_along(order)) {
  which <- order[k]
  cat(sprintf("[%s] warm run %d/%d: %s\n", format(Sys.time(), "%H:%M:%S"), k, length(order), which))
  od <- file.path(opt$out, sprintf("warm%d_%s", k, which))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  fit <- mods[[which]]$sample(
    data = data, init = inits, seed = opt$seed, chains = opt$chains,
    parallel_chains = opt$chains, iter_warmup = opt$warmup, iter_sampling = opt$sampling,
    inv_metric = inv_metric, step_size = 0.1, max_treedepth = opt$treedepth, save_warmup = TRUE,
    refresh = 50, show_messages = TRUE, output_dir = od
  )
  tm <- fit$time()$chains
  sd_all <- fit$sampler_diagnostics(inc_warmup = TRUE)
  sd <- fit$sampler_diagnostics()
  nleap <- as.numeric(sd[, , "n_leapfrog__"])
  nleap_w <- as.numeric(sd_all[seq_len(opt$warmup), , "n_leapfrog__"])
  s <- summarise_draws(fit$draws(variables = pars), "rhat", "ess_bulk", "ess_tail")
  row <- data.frame(
    run = k, tree = which, warmup = opt$warmup, sampling = opt$sampling, chains = opt$chains,
    warmup_cpu_s = sum(tm$warmup), sampling_cpu_s = sum(tm$sampling), wall_s = max(tm$total),
    leapfrog_warmup_total = sum(nleap_w), leapfrog_sampling_total = sum(nleap),
    leapfrog_sampling_mean = mean(nleap),
    accept_mean = mean(as.numeric(sd[, , "accept_stat__"])),
    stepsize = mean(as.numeric(sd[1, , "stepsize__"])),
    divergent = sum(as.numeric(sd[, , "divergent__"])),
    treedepth_cap_hits = sum(as.numeric(sd[, , "treedepth__"]) >= opt$treedepth),
    rhat_max = max(s$rhat), ess_bulk_min = min(s$ess_bulk),
    ess_bulk_min_par = s$variable[which.min(s$ess_bulk)],
    ess_bulk_mean = mean(s$ess_bulk), ess_tail_min = min(s$ess_tail)
  )
  row$us_per_grad <- 1e6 * (row$warmup_cpu_s + row$sampling_cpu_s) /
    (row$leapfrog_warmup_total + row$leapfrog_sampling_total)
  row$ess_min_per_cpu_s <- row$ess_bulk_min / row$sampling_cpu_s
  row$ess_min_per_total_cpu_s <- row$ess_bulk_min / (row$warmup_cpu_s + row$sampling_cpu_s)
  row$ess_min_per_1000_grad <- 1000 * row$ess_bulk_min / row$leapfrog_sampling_total
  print(row, digits = 4)
  print(as.data.frame(s), digits = 4)
  rows[[k]] <- row
  write.csv(do.call(rbind, rows), file.path(opt$out, "ess_warm.csv"), row.names = FALSE)
}
cat("done\n")
