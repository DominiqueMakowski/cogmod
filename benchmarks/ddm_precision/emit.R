# Program emission (and an abandoned cold-start fit) for the 7-parameter DDM,
# from the same brms-generated program, data and start on two trees that
# differ only in .DDM_WIENER_PRECISION (base 1e-4, PR 1e-3).
#
#   Rscript ess_bench.R emit --pkg <tree> --out <dir> [--n 300]
#   Rscript ess_bench.R run --base <dir> --pr <dir> --out <dir>
#            [--warmup 400 --sampling 400 --chains 4 --order base,pr,pr,base]
#
# `run` samples the two programs in the given order (ABBA by default, so a
# drift in machine load over the hour hits both alike), each with `chains`
# chains in parallel, and writes <out>/ess_bench.csv with one row per run:
# wall time, leapfrogs, step size, divergences, and bulk / tail ESS of the
# eight intercepts, min over parameters, per second of sampling and per
# gradient evaluation.
args <- commandArgs(trailingOnly = TRUE)
mode <- args[1]
opt <- list(pkg = ".", out = ".", n = 300L, base = "", pr = "", warmup = 400L,
            sampling = 400L, chains = 4L, order = "base,pr,pr,base", seed = 1L,
            treedepth = 10L)
a <- args[-1]
for (i in seq(1, length(a), by = 2)) {
  key <- sub("^--", "", a[i])
  opt[[key]] <- if (is.numeric(opt[[key]])) as.numeric(a[i + 1]) else a[i + 1]
}

if (mode == "emit") {
  source(file.path(getwd(), "benchmarks/gradient_programs.R"))
  gp_load(opt$pkg)
  gp_emit(opt$out, "cogmod_ddm", n = as.integer(opt$n), extremes = FALSE)
  cat("precision in program:",
      regmatches(readLines(file.path(opt$out, "cogmod_ddm.stan")),
                 regexpr("sigmadrift, sw, sigmandt[^;]*;", readLines(file.path(opt$out, "cogmod_ddm.stan")))),
      "\n")
  quit(status = 0)
}

suppressPackageStartupMessages({ library(cmdstanr); library(posterior) })
dirs <- c(base = opt$base, pr = opt$pr)
mods <- lapply(dirs, function(d) {
  dir.create(file.path(d, "exe"), showWarnings = FALSE)
  cmdstan_model(file.path(d, "cogmod_ddm.stan"), dir = file.path(d, "exe"))
})
data <- file.path(dirs[["base"]], "cogmod_ddm.data.json")
init <- file.path(dirs[["base"]], "cogmod_ddm.init.json")
stopifnot(identical(readLines(data), readLines(file.path(dirs[["pr"]], "cogmod_ddm.data.json"))))

order <- strsplit(opt$order, ",")[[1]]
rows <- list()
for (k in seq_along(order)) {
  which <- order[k]
  cat(sprintf("[%s] run %d/%d: %s\n", format(Sys.time(), "%H:%M:%S"), k, length(order), which))
  dir.create(file.path(opt$out, sprintf("run%d_%s", k, which)), recursive = TRUE, showWarnings = FALSE)
  fit <- mods[[which]]$sample(
    data = data, init = init, seed = opt$seed, chains = opt$chains,
    parallel_chains = opt$chains, iter_warmup = opt$warmup,
    iter_sampling = opt$sampling, max_treedepth = opt$treedepth,
    refresh = 0, show_messages = FALSE,
    output_dir = file.path(opt$out, sprintf("run%d_%s", k, which))
  )
  tm <- fit$time()$chains
  sd <- fit$sampler_diagnostics()
  nleap <- as.numeric(sd[, , "n_leapfrog__"])
  drw <- fit$draws(variables = c("Intercept", paste0("Intercept_", c(
    "boundary", "bias", "sigmadrift", "sigmabias", "sigmandt", "ndt", "poutlier"))))
  s <- summarise_draws(drw, "rhat", "ess_bulk", "ess_tail")
  samp_time <- sum(tm$sampling)   # CPU seconds across chains
  wall <- max(tm$total)
  row <- data.frame(
    run = k, tree = which, chains = opt$chains, treedepth_cap = opt$treedepth,
    wall_s = wall, warmup_cpu_s = sum(tm$warmup), sampling_cpu_s = samp_time,
    leapfrog_mean = mean(nleap), leapfrog_total = sum(nleap),
    stepsize = mean(as.numeric(sd[1, , "stepsize__"])),
    divergent = sum(as.numeric(sd[, , "divergent__"])),
    max_treedepth_hits = sum(as.numeric(sd[, , "treedepth__"]) >= opt$treedepth),
    rhat_max = max(s$rhat), ess_bulk_min = min(s$ess_bulk), ess_tail_min = min(s$ess_tail),
    ess_bulk_min_par = s$variable[which.min(s$ess_bulk)],
    ess_bulk_mean = mean(s$ess_bulk)
  )
  row$ess_per_cpu_s <- row$ess_bulk_min / row$sampling_cpu_s
  row$ess_per_wall_s <- row$ess_bulk_min / sum(tm$sampling) * opt$chains  # parallel chains
  row$ess_per_1000_grad <- 1000 * row$ess_bulk_min / row$leapfrog_total
  row$us_per_grad <- 1e6 * (sum(tm$warmup) + samp_time) / row$leapfrog_total
  print(row, digits = 4)
  print(as.data.frame(s), digits = 4)
  rows[[k]] <- row
  write.csv(do.call(rbind, rows), file.path(opt$out, "ess_bench.csv"), row.names = FALSE)
}
cat("done\n")
