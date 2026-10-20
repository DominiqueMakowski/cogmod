# Aggregate chains_probe.R runs (`run.sh submit mlchains`, pulled into
# results/hpc/mlchains_<tag>_<model>_<seed>/).
#
#   Rscript benchmarks/stanli/hpc/summarise_chains.R [--dir benchmarks/stanli/results/hpc] [--tag 0.20.0]
#
# Per task, the sampling-phase cost of a leapfrog step, pooled over chains,
# relative to the same engine's single chain alone (`alone`) and to CmdStan
# in the same configuration (`vs_cmdstan`); then medians [min, max] over
# seeds, per model and CPU type.
source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(dir = "benchmarks/stanli/results/hpc", tag = "0.20.0"))
f <- Sys.glob(file.path(opt$dir, paste0("mlchains_", opt$tag, "_*"), "ml_chainsprobe_*.csv"))
if (!length(f)) stop("no mlchains results under ", opt$dir, call. = FALSE)
res <- do.call(rbind, lapply(f, function(x) {
  nt <- readLines(file.path(dirname(x), "node.txt"))
  cbind(cpu = sub(".*EPYC ([0-9]+).*", "\\1", grep("model name", nt, value = TRUE)[1]),
        host = sub("^node: ([^.]+).*", "\\1", nt[1]), utils::read.csv(x))
}))
s <- aggregate(cbind(sample_s, leapfrog) ~ cpu + host + model + seed + config, res, sum)
s$ms <- 1000 * s$sample_s / s$leapfrog
w <- stats::reshape(s[c("cpu", "host", "model", "seed", "config", "ms")],
                    idvar = c("cpu", "host", "model", "seed"), timevar = "config", direction = "wide")
names(w) <- sub("^ms\\.", "", names(w))
w$cm_4_alone <- w$cmdstan_4 / w$cmdstan_1
w$sl_4_alone <- w$stanli_4 / w$stanli_1
w$sl_4proc_alone <- w$stanli_4proc / w$stanli_1
w$sl_vs_cm_1 <- w$stanli_1 / w$cmdstan_1
w$sl_vs_cm_4 <- w$stanli_4 / w$cmdstan_4
w$sl4proc_vs_cm_4 <- w$stanli_4proc / w$cmdstan_4
mmm <- function(x) sprintf("%.2f [%.2f, %.2f]", stats::median(x), min(x), max(x))
cols <- c("sl_vs_cm_1", "sl_vs_cm_4", "sl4proc_vs_cm_4", "cm_4_alone", "sl_4_alone", "sl_4proc_alone")
tab <- do.call(rbind, lapply(split(w, list(w$model, w$cpu), drop = TRUE), function(x)
  data.frame(model = x$model[1], cpu = x$cpu[1], seeds = nrow(x), lapply(x[cols], mmm))))
cat("Sampling-phase ms per leapfrog, ratios (median [min, max] over seeds)\n",
    "  sl_vs_cm_1: stanli / CmdStan, one chain alone; sl_vs_cm_4: four chains each;\n",
    "  sl4proc_vs_cm_4: stanli as four processes / CmdStan four chains;\n",
    "  *_alone: four chains / the same engine's one chain\n\n", sep = "")
print(tab[order(tab$model, tab$cpu), ], row.names = FALSE, right = FALSE)
