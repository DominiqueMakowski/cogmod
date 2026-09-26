# Compact view of summary_long.csv for one experiment: one line per group
# with data / constant ratios (and paired p) of the metrics that matter, then
# the pooled win shares by data set. Usage: Rscript view.R results A
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results"
exp_id <- if (length(args) > 1) args[2] else "A"
long <- read.csv(file.path(dir, "summary_long.csv"), stringsAsFactors = FALSE)
long <- long[long$exp == exp_id, ]
show <- c("warmup_leapfrog", "transient_max", "warmup_s", "min_ess_per_total_s", "median_ess_per_s", "max_rhat")
options(width = 220)
key <- c("data", "family", "warmup")
groups <- unique(long[, key])
rows <- lapply(seq_len(nrow(groups)), function(i) {
  s <- merge(long, groups[i, , drop = FALSE])
  r <- groups[i, , drop = FALSE]
  r$n <- s$n_seeds[1]
  for (m in show) {
    z <- s[s$metric == m, ]
    r[[m]] <- if (nrow(z)) sprintf("%.2f%s", z$ratio, ifelse(is.na(z$p_paired), "", ifelse(z$p_paired < 0.01, "**", ifelse(z$p_paired < 0.05, "*", "")))) else NA
  }
  r
})
tab <- do.call(rbind, rows)
tab <- tab[order(tab$data, tab$family, tab$warmup), ]
cat(sprintf("Experiment %s: ratio data / constant, * p<.05, ** p<.01 (paired over seeds)\n", exp_id))
print(tab, row.names = FALSE)

cat("\nWin share of the data start, by data set (ratio < 1 for costs, > 1 for ESS):\n")
for (d in unique(long$data)) {
  cat(" ", d, "\n")
  for (m in c("warmup_leapfrog", "transient_max", "warmup_s", "min_ess_per_total_s", "median_ess_per_s")) {
    z <- long[long$data == d & long$metric == m, ]
    better <- if (grepl("ess", m)) z$ratio > 1 else z$ratio < 1
    cat(sprintf("    %-22s %2d / %2d   median ratio %.3f   min %.2f  max %.2f\n", m, sum(better, na.rm = TRUE), sum(!is.na(better)), median(z$ratio, na.rm = TRUE), min(z$ratio, na.rm = TRUE), max(z$ratio, na.rm = TRUE)))
  }
}
