# Read the result rows and compare the two schemes cell by cell: for every
# (experiment, data, family, warmup) the mean over seeds of each metric under
# each scheme, the paired difference, and a paired t-test across seeds.
# Usage: Rscript summarise.R [results_dir]
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results"
files <- list.files(dir, pattern = "^[AB]_[0-9]+[.]csv$", full.names = TRUE)
stopifnot(length(files) > 0)
res <- do.call(rbind, lapply(files, function(f) {
  r <- read.csv(f, stringsAsFactors = FALSE)
  r[, setdiff(names(r), "node"), drop = FALSE]
}))
res <- res[order(res$exp, res$id), ]
cat(nrow(res), "rows;", sum(res$status != "ok"), "not ok\n")
if (any(res$status != "ok")) print(res[res$status != "ok", c("exp", "id", "data", "family", "warmup", "scheme", "seed", "status")])
ok <- res[res$status == "ok", ]

metrics <- c("warmup_s", "sample_s", "warmup_leapfrog", "warmup_maxtree_frac",
             "transient_max", "transient_leapfrog", "sample_leapfrog_mean",
             "divergences", "max_rhat", "median_rhat", "min_ess_bulk",
             "median_ess_bulk", "min_ess_per_s", "median_ess_per_s",
             "min_ess_per_total_s")
metrics <- intersect(metrics, names(ok))

key <- c("exp", "data", "family", "warmup")
groups <- unique(ok[, key])
out <- list()
for (i in seq_len(nrow(groups))) {
  g <- merge(ok, groups[i, , drop = FALSE])
  a <- g[g$scheme == "constant", ]
  b <- g[g$scheme == "data", ]
  seeds <- intersect(a$seed, b$seed)
  a <- a[match(seeds, a$seed), ]
  b <- b[match(seeds, b$seed), ]
  for (m in metrics) {
    x <- a[[m]]; y <- b[[m]]
    d <- y - x
    p <- if (length(seeds) >= 3 && sd(d, na.rm = TRUE) > 0) tryCatch(t.test(y, x, paired = TRUE)$p.value, error = function(e) NA) else NA
    out[[length(out) + 1]] <- data.frame(
      groups[i, , drop = FALSE], metric = m, n_seeds = length(seeds),
      constant = mean(x, na.rm = TRUE), data = mean(y, na.rm = TRUE),
      ratio = mean(y, na.rm = TRUE) / mean(x, na.rm = TRUE),
      diff_mean = mean(d, na.rm = TRUE), diff_sd = sd(d, na.rm = TRUE),
      p_paired = p, stringsAsFactors = FALSE
    )
  }
}
out <- do.call(rbind, out)
write.csv(out, file.path(dir, "summary_long.csv"), row.names = FALSE)

# The headline table: one line per group, the ratios data / constant of the
# metrics that would show a real gain.
show <- c("warmup_leapfrog", "transient_max", "min_ess_per_total_s", "median_ess_per_s", "max_rhat")
wide <- do.call(rbind, lapply(split(out, out[, key], drop = TRUE), function(s) {
  r <- s[1, key]
  for (m in show) {
    z <- s[s$metric == m, ]
    r[[paste0(m, "_ratio")]] <- if (nrow(z)) round(z$ratio, 3) else NA
    r[[paste0(m, "_p")]] <- if (nrow(z)) round(z$p_paired, 3) else NA
  }
  r$n_seeds <- s$n_seeds[1]
  r
}))
wide <- wide[order(wide$exp, wide$data, wide$family, wide$warmup), ]
options(width = 200)
print(wide, row.names = FALSE)
write.csv(wide, file.path(dir, "summary.csv"), row.names = FALSE)

# Pooled over everything: how often the data start won on each metric.
cat("\nShare of groups where data start beat constant (ratio < 1 for costs, > 1 for ESS):\n")
for (m in metrics) {
  z <- out[out$metric == m, ]
  better <- if (grepl("ess", m)) z$ratio > 1 else z$ratio < 1
  cat(sprintf("  %-22s %3d / %3d  median ratio %.3f\n", m, sum(better, na.rm = TRUE), sum(!is.na(better)), median(z$ratio, na.rm = TRUE)))
}
