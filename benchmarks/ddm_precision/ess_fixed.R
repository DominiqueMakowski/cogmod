# Stage 2: sampling only, from one shared adapted state, adaptation off.
#
#   Rscript ess_fixed.R --from <dir of a finished run's CSVs> --base <dir> --pr <dir>
#                       --out <dir> [--iter 500 --chains 4 --order base,pr,pr,base]
#
# The step size (median over chains) and diagonal inverse metric (mean over
# chains) come from the finished run; each chain starts at that run's last
# draw of the same chain. Both trees then sample with identical settings, so
# every difference is the density and its gradient. Reports acceptance,
# leapfrogs, divergences, ESS per second and per gradient evaluation.
args <- commandArgs(trailingOnly = TRUE)
opt <- list(from = "", base = "", pr = "", out = ".", iter = 500L, chains = 4L,
            order = "base,pr,pr,base", seed = 7L, treedepth = 10L)
for (i in seq(1, length(args), by = 2)) {
  key <- sub("^--", "", args[i])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(args[i + 1]) else args[i + 1]
}
suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })

csvs <- list.files(opt$from, pattern = "^cogmod_ddm-\\d+-\\d-[0-9a-f]+\\.csv$", full.names = TRUE)
stopifnot(length(csvs) == opt$chains)
src <- read_cmdstan_csv(csvs)
step <- median(unlist(src$step_size))
inv_metric <- Reduce(`+`, src$inv_metric) / length(src$inv_metric)
pars <- c("Intercept", paste0("Intercept_", c("boundary", "bias", "sigmadrift",
                                             "sigmabias", "sigmandt", "ndt", "poutlier")))
n_it <- dim(src$post_warmup_draws)[1]
inits <- lapply(seq_len(opt$chains), function(ch) {
  as.list(setNames(as.numeric(src$post_warmup_draws[n_it, ch, pars]), pars))
})
cat(sprintf("adapted step size %.4g (chains: %s)\n", step, paste(signif(unlist(src$step_size), 3), collapse = ", ")))
cat("inverse metric:", signif(inv_metric, 3), "\n")

dirs <- c(base = opt$base, pr = opt$pr)
mods <- lapply(dirs, function(d) cmdstan_model(file.path(d, "cogmod_ddm.stan"), dir = file.path(d, "exe")))
data <- file.path(dirs[["base"]], "cogmod_ddm.data.json")
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

order <- strsplit(opt$order, ",")[[1]]
rows <- list()
for (k in seq_along(order)) {
  which <- order[k]
  cat(sprintf("[%s] fixed run %d/%d: %s\n", format(Sys.time(), "%H:%M:%S"), k, length(order), which))
  od <- file.path(opt$out, sprintf("fixed%d_%s", k, which))
  dir.create(od, recursive = TRUE, showWarnings = FALSE)
  fit <- mods[[which]]$sample(
    data = data, init = inits, seed = opt$seed, chains = opt$chains,
    parallel_chains = opt$chains, iter_warmup = 0, iter_sampling = opt$iter,
    adapt_engaged = FALSE, step_size = step, inv_metric = inv_metric,
    max_treedepth = opt$treedepth, refresh = 0, show_messages = FALSE, output_dir = od
  )
  tm <- fit$time()$chains
  sd <- fit$sampler_diagnostics()
  nleap <- as.numeric(sd[, , "n_leapfrog__"])
  s <- summarise_draws(fit$draws(variables = pars), "rhat", "ess_bulk", "ess_tail")
  row <- data.frame(
    run = k, tree = which, iter = opt$iter, chains = opt$chains,
    sampling_cpu_s = sum(tm$sampling), wall_s = max(tm$total),
    accept_mean = mean(as.numeric(sd[, , "accept_stat__"])),
    leapfrog_mean = mean(nleap), leapfrog_total = sum(nleap),
    treedepth_mean = mean(as.numeric(sd[, , "treedepth__"])),
    divergent = sum(as.numeric(sd[, , "divergent__"])),
    rhat_max = max(s$rhat), ess_bulk_min = min(s$ess_bulk),
    ess_bulk_min_par = s$variable[which.min(s$ess_bulk)],
    ess_bulk_mean = mean(s$ess_bulk), ess_tail_min = min(s$ess_tail)
  )
  row$us_per_grad <- 1e6 * row$sampling_cpu_s / row$leapfrog_total
  row$ess_min_per_cpu_s <- row$ess_bulk_min / row$sampling_cpu_s
  row$ess_mean_per_cpu_s <- row$ess_bulk_mean / row$sampling_cpu_s
  row$ess_min_per_1000_grad <- 1000 * row$ess_bulk_min / row$leapfrog_total
  print(row, digits = 4)
  print(as.data.frame(s), digits = 4)
  rows[[k]] <- row
  write.csv(do.call(rbind, rows), file.path(opt$out, "ess_fixed.csv"), row.names = FALSE)
}
cat("done\n")
