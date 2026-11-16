# Aggregate bench.R's runs over seeds: one row per (model, arm).
#
#   Rscript benchmarks/nutpieR/summarise.R [--dir benchmarks/nutpieR/results/hpc]
#                                          [--tag NAME] [--out FILE.csv]
#
# Reads every fit_*.csv under --dir: the cluster's fit[_<tag>]_<seed>/
# directories that `hpc/run.sh pull` brings back, or, with --dir
# benchmarks/nutpieR/results, a local run's fit_<model>_<seed>[_<tag>].csv.
# Ratios are against cmdstanr within each seed, where the arms ran one after
# the other on the same node, so they are paired; wall times across seeds are
# not comparable (shared nodes, several CPU types). Leapfrog steps per
# effective draw need no clock. Cells are median [min, max] over seeds.
source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(dir = "benchmarks/nutpieR/results/hpc", tag = "", out = ""))

f <- c(Sys.glob(file.path(opt$dir, "fit_*", "fit_*.csv")), Sys.glob(file.path(opt$dir, "fit_*.csv")))
# fit_<model>_<seed>[_<tag>].csv: the runs of one tag (none by default).
file_tag <- sub("^fit_[a-z0-9]+_[0-9]+_?", "", sub("\\.csv$", "", basename(f)))
f <- f[file_tag == opt$tag]
if (!length(f)) stop("no fit_*.csv under ", opt$dir, if (nzchar(opt$tag)) paste0(" with tag ", opt$tag), call. = FALSE)
fit <- do.call(rbind, lapply(f, function(x) cbind(dir = basename(dirname(x)), utils::read.csv(x))))
# The chains and pars files sit next to the fit files, with the same stem.
ch <- do.call(rbind, lapply(sub("/fit_([^/]*)$", "/chains_\\1", f), utils::read.csv))
pars <- do.call(rbind, lapply(sub("/fit_([^/]*)$", "/pars_\\1", f), utils::read.csv))
# r_1[j, k] duplicates r_Participant[j, ...] in runs before bench.R dropped it.
pars <- pars[!grepl("^r_[0-9]+\\[", pars$variable), ]

nodes <- unique(fit$dir[fit$dir != basename(opt$dir)])
if (length(nodes)) {
  cpu <- vapply(nodes, function(d) {
    p <- file.path(opt$dir, d, "node.txt")
    if (!file.exists(p)) return(NA_character_)
    x <- readLines(p)
    sub(".*: ", "", grep("model name", x, value = TRUE)[1])
  }, character(1))
  cat("seeds per CPU:\n")
  print(table(cpu))
}
v <- unique(fit[grep("^v_", names(fit))])
cat("versions:", paste(apply(v, 1, function(r) paste(sub("^v_", "", names(v)), r, collapse = " ")), collapse = " | "), "\n")

mmm <- function(x, d = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  if (length(x) == 1) return(formatC(x, format = "f", digits = d))
  sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), stats::median(x), min(x), max(x))
}
arm_levels <- c("cmdstanr", "stanli", "nutpie", "nutpie_lowrank")
cells <- unique(fit[c("model", "arm")])
cells <- cells[order(cells$model, match(cells$arm, arm_levels)), ]
tab <- do.call(rbind, lapply(seq_len(nrow(cells)), function(i) {
  x <- fit[fit$model == cells$model[i] & fit$arm == cells$arm[i], ]
  data.frame(
    model = x$model[1], arm = x$arm[1], seeds = nrow(x), N = x$N[1], n_par = x$n_par[1],
    warmup = paste(unique(x$warmup_iter), collapse = "/"),
    build_s = mmm(x$build_s, 1), wall_s = mmm(x$wall_s, 1),
    ess_bulk_min = mmm(x$ess_bulk_min, 0), ess_tail_min = mmm(x$ess_tail_min, 0),
    ess_per_wall_s = mmm(x$ess_bulk_per_wall_s, 1),
    wall_vs_cmdstanr = mmm(x$wall_vs_cmdstanr),
    ess_per_wall_s_vs_cmdstanr = mmm(x$ess_per_wall_s_vs_cmdstanr),
    ess_tail_per_wall_s_vs_cmdstanr = mmm(x$ess_tail_per_wall_s_vs_cmdstanr),
    ess_per_total_s_vs_cmdstanr = mmm(x$ess_per_total_s_vs_cmdstanr),
    leapfrog_warmup_k = mmm(x$leapfrog_warmup / 1000, 1), leapfrog_k = mmm(x$leapfrog / 1000, 1),
    leapfrog_per_ess = mmm(x$leapfrog_per_ess, 0),
    ess_1k_lf = mmm(x$ess_bulk_per_1k_leapfrog, 1), ms_per_leapfrog = mmm(x$ms_per_leapfrog, 3),
    stepsize = mmm(x$stepsize, 3), rhat_max = round(max(x$rhat_max), 3),
    divergent = sum(x$divergent), seeds_divergent = sum(x$divergent > 0),
    mean_z_max = round(max(x$mean_z_max), 1))
}))

show <- function(cols, title) {
  cat("\n", title, "\n", sep = "")
  print(tab[c("model", "arm", cols)], row.names = FALSE, right = FALSE)
}
show(c("seeds", "N", "n_par", "warmup", "build_s", "wall_s", "ess_bulk_min", "ess_per_wall_s"),
     "Fits (wall = sampling to brmsfit; build = compile or stanli build; ESS = min bulk over population-level parameters)")
show(c("wall_vs_cmdstanr", "ess_per_wall_s_vs_cmdstanr", "ess_tail_per_wall_s_vs_cmdstanr",
       "ess_per_total_s_vs_cmdstanr"),
     "Against cmdstanr within seed (> 1 = more ESS per second; total = build + wall)")
show(c("leapfrog_warmup_k", "leapfrog_k", "leapfrog_per_ess", "ess_1k_lf", "ms_per_leapfrog", "stepsize"),
     "Sampler efficiency without a clock (leapfrog steps over all chains, thousands; per ESS = warmup + sampling over min bulk ESS)")
show(c("rhat_max", "divergent", "seeds_divergent", "mean_z_max"),
     "Checks (mean_z_max = largest |mean difference| from cmdstanr / joint MCSE, over seeds)")
if (nzchar(opt$out)) utils::write.csv(tab, opt$out, row.names = FALSE)

# Stragglers. On a shared node a chain can run at several times its arm's
# usual cost per leapfrog step for a whole phase, with the same number of
# steps (README, "Stragglers"): that is the node, not the sampler, and it
# lengthens the arm's wall time. A seed is flagged when any chain of an arm
# that reports phase times (cmdstanr, stanli; nutpie times the run as one)
# spends more than 1.5x the median of that (model, arm) per step in warmup or
# sampling. Everything below is per model: a step's cost differs by model.
ch$ms_w <- 1000 * ch$warmup_s / ch$leapfrog_warmup
ch$ms_s <- 1000 * ch$sample_s / ch$leapfrog
key <- function(x) paste(x$model, x$seed)
lf <- aggregate(cbind(steps = leapfrog_warmup + leapfrog) ~ model + seed + arm, ch, max)
lf <- merge(lf, fit[c("model", "seed", "arm", "wall_s")])
for (m in sort(unique(fit$model))) {
  timed <- ch[ch$model == m & is.finite(ch$ms_w), ]
  ma <- paste(timed$arm)
  med <- tapply(c(timed$ms_w, timed$ms_s), c(ma, ma), median)
  slow <- timed$ms_w > 1.5 * med[ma] | timed$ms_s > 1.5 * med[ma]
  flag <- sort(unique(timed$seed[slow]))
  fm <- fit[fit$model == m, ]
  clean <- fm[!fm$seed %in% flag, ]
  cat(sprintf("\n== %s. Stragglers in %d of %d seeds%s\n", m, length(flag), length(unique(fm$seed)),
              if (length(flag)) paste0(" (", paste(flag, collapse = " "), ")") else ""))
  if (length(flag) && nrow(clean)) {
    cat("Without them (", length(unique(clean$seed)), " seeds), against cmdstanr:\n", sep = "")
    for (a in setdiff(unique(clean$arm), "cmdstanr")) {
      x <- clean[clean$arm == a, ]
      cat(sprintf("  %-15s wall %s | ESS/s %s | tail ESS/s %s\n", a, mmm(x$wall_vs_cmdstanr),
                  mmm(x$ess_per_wall_s_vs_cmdstanr), mmm(x$ess_tail_per_wall_s_vs_cmdstanr)))
    }
  }
  # nutpie against stanli in the same seed: the same gradient, so only the
  # sampler differs (all seeds: stragglers hit CmdStan's processes). Wall
  # time per step of the longest chain (warmup and sampling, fixed costs
  # included) says whether nutpie's steps cost more.
  for (a in intersect(c("nutpie", "nutpie_lowrank"), unique(fm$arm))) {
    w <- merge(fm[fm$arm == "stanli", ], fm[fm$arm == a, ], by = "seed")
    if (!nrow(w)) next
    cat(sprintf("%s / stanli within seed (%d seeds): wall %s | ESS/s %s | leapfrog per ESS %s | min bulk ESS %s\n",
                a, nrow(w), mmm(w$wall_s.y / w$wall_s.x), mmm(w$ess_bulk_per_wall_s.y / w$ess_bulk_per_wall_s.x),
                mmm(w$leapfrog_per_ess.y / w$leapfrog_per_ess.x), mmm(w$ess_bulk_min.y / w$ess_bulk_min.x)))
  }
  l <- lf[lf$model == m & !lf$seed %in% flag, ]
  if (nrow(l)) {
    s <- tapply(1000 * l$wall_s / l$steps, paste(l$arm), function(v) mmm(v, 3))
    cat("Wall ms per step of the longest chain, seeds without stragglers:",
        paste(names(s), s, collapse = " | "), "\n")
  }
  # Posterior agreement with cmdstanr, every parameter in every seed:
  # z = mean difference / joint MCSE should be N(0, 1). sd(z) above 1 means
  # the MCSEs are too small (ESS optimistic); a parameter whose z keeps one
  # sign across seeds (mean z x sqrt(seeds) far from 0) would be a bias.
  pm <- pars[pars$model == m, ]
  for (a in setdiff(unique(pm$arm), "cmdstanr")) {
    w <- merge(pm[pm$arm == "cmdstanr", ], pm[pm$arm == a, ], by = c("seed", "variable"))
    if (!nrow(w)) next
    w$z <- (w$mean.y - w$mean.x) / sqrt(w$mcse_mean.x^2 + w$mcse_mean.y^2)
    v <- aggregate(cbind(z, sd.x, sd.y) ~ variable, w, mean)
    cat(sprintf("Agreement, %s: |z| > 3 in %d of %d (%.1f by chance) | sd(z) %.2f | max |mean z| x sqrt(seeds) %.2f over %d parameters | posterior sd ratio %.3f-%.3f\n",
                a, sum(abs(w$z) > 3), nrow(w), 0.0027 * nrow(w), stats::sd(w$z),
                max(abs(v$z)) * sqrt(length(unique(w$seed))), nrow(v), min(v$sd.y / v$sd.x), max(v$sd.y / v$sd.x)))
  }
}
