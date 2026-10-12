# Aggregate the multilevel runs (multilevel.R, `run.sh submit ml`) that
# `run.sh pull` brought back.
#
#   Rscript benchmarks/stanli/hpc/summarise_ml.R [--dir benchmarks/stanli/results/hpc]
#                                                [--tag 0.20.0] [--out FILE.csv]
#
# Reads every ml_[<tag>_]<model>_<seed>/ under --dir (or, with --dir pointing
# at benchmarks/stanli/results, the local ml_*_<model>_<seed>[_<tag>].csv).
# One row per model over its seeds: median [min, max] of the within-seed
# ratios stanli / cmdstanr. Within a seed both arms ran one after the other
# on the same node, so the ratios are paired; wall times across seeds are
# not comparable (shared nodes, three CPU types). ESS per 1000 leapfrog
# steps needs no clock. --out writes the table as a CSV as well.
source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(dir = "benchmarks/stanli/results/hpc", tag = "", out = ""))
glob <- function(stem) {
  pre <- if (nzchar(opt$tag)) paste0("ml_", opt$tag, "_*") else "ml_*"
  f <- c(Sys.glob(file.path(opt$dir, pre, paste0("ml_", stem, "_*.csv"))),
         Sys.glob(file.path(opt$dir, paste0("ml_", stem, "_*.csv"))))
  if (nzchar(opt$tag)) f <- f[grepl(paste0("_", opt$tag, "\\.csv$"), f)]
  if (!length(f)) stop("no ml_", stem, "_*.csv under ", opt$dir, call. = FALSE)
  do.call(rbind, lapply(f, function(x) cbind(dir = basename(dirname(x)), utils::read.csv(x))))
}
fit <- glob("fit")
tm <- glob("time")
ck <- glob("check")
# The cost of a leapfrog step in the sampling phase, from the chains: sampling
# time over sampling leapfrogs (the leapfrogs nuts_params() counts). Runs
# before 2026-10-12 evening divided warmup + sampling time by them.
ch <- glob("chains")
lf <- aggregate(cbind(sample_s, leapfrog) ~ model + seed + arm, ch, sum)
fit$ms_per_leapfrog <- NULL
fit <- merge(fit, transform(lf, ms_per_leapfrog = 1000 * sample_s / leapfrog)[
  c("model", "seed", "arm", "ms_per_leapfrog")])
order_models <- c("lnr", "lnr6", "ddm4", "ddm5", "ddm7", "lba", "rdm")
models <- intersect(order_models, unique(fit$model))

# Nodes, where the run came from the cluster.
nodes <- unique(fit$dir[fit$dir != basename(opt$dir)])
if (length(nodes)) {
  cpu <- vapply(nodes, function(d) {
    p <- file.path(opt$dir, d, "node.txt")
    if (!file.exists(p)) return(NA_character_)
    x <- readLines(p)
    paste(sub("^node: ", "", x[1]), sub(".*: ", "", grep("model name", x, value = TRUE)[1]))
  }, character(1))
  cat("tasks per CPU:\n")
  print(table(sub("^\\S+ ", "", cpu)))
}

mmm <- function(x, d = 2) {
  x <- x[is.finite(x)]
  if (!length(x)) return(NA_character_)
  if (length(x) == 1) return(formatC(x, format = "f", digits = d))
  sprintf(paste0("%.", d, "f [%.", d, "f, %.", d, "f]"), stats::median(x), min(x), max(x))
}
# One value per seed: the stanli row's column over the cmdstanr row's.
paired <- function(x, col) {
  w <- merge(x[x$arm == "cmdstanr", c("seed", col)], x[x$arm == "stanli", c("seed", col)], by = "seed")
  w[[3]] / w[[2]]
}

tab <- do.call(rbind, lapply(models, function(m) {
  x <- fit[fit$model == m, ]
  t <- tm[tm$model == m, ]
  # Per gradient: within each seed (one node), the ratio within each block,
  # its median; then over seeds.
  g <- do.call(rbind, lapply(split(t, t$seed), function(s) {
    w <- stats::reshape(s[c("rep", "engine", "us_per_grad")], idvar = "rep", timevar = "engine",
                        direction = "wide")
    data.frame(cm = stats::median(w$us_per_grad.cmdstan), sl = stats::median(w$us_per_grad.stanli),
               r = stats::median(w$us_per_grad.stanli / w$us_per_grad.cmdstan))
  }))
  c_ <- ck[ck$model == m & ck$point != "init", ]
  cm <- x[x$arm == "cmdstanr", ]
  sl <- x[x$arm == "stanli", ]
  data.frame(
    model = m, seeds = length(unique(x$seed)), N = x$N[1], n_par = x$n_par[1],
    grad_us_cmdstan = round(stats::median(g$cm)), grad_us_stanli = round(stats::median(g$sl)),
    grad_ratio = mmm(g$r),
    build_s_cmdstanr = round(stats::median(cm$build_s)), build_s_stanli = round(stats::median(sl$build_s), 1),
    wall_s_cmdstanr = round(stats::median(cm$wall_s)), wall_s_stanli = round(stats::median(sl$wall_s)),
    wall_ratio = mmm(paired(x, "wall_s")),
    ms_per_leapfrog_ratio = mmm(paired(x, "ms_per_leapfrog")),
    ess_per_wall_s_ratio = mmm(paired(x, "ess_bulk_per_wall_s")),
    ess_1k_lf_cmdstanr = mmm(cm$ess_bulk_per_1k_leapfrog, 1),
    ess_1k_lf_stanli = mmm(sl$ess_bulk_per_1k_leapfrog, 1),
    leapfrog_ratio = mmm(paired(x, "leapfrog")),
    rhat_max = round(max(x$rhat_max), 3), divergent = paste(sum(cm$divergent), sum(sl$divergent), sep = "/"),
    mean_z_max = round(max(x$mean_z_max), 1),
    lp_identical = sprintf("%d/%d", sum(c_$lp_ulp == 0, na.rm = TRUE), nrow(c_)),
    lp_ulp_max = signif(max(c_$lp_ulp, na.rm = TRUE), 2),
    grad_relerr_max = signif(max(c_$grad_relerr, na.rm = TRUE), 2),
    nonfinite = sum(!c_$grad_finite_cmdstan | !c_$grad_finite_stanli))
}))

show <- function(cols, title) {
  cat("\n", title, "\n", sep = "")
  print(tab[c("model", cols)], row.names = FALSE, right = FALSE)
}
show(c("seeds", "N", "n_par", "grad_us_cmdstan", "grad_us_stanli", "grad_ratio"),
     "Per gradient at the posterior (single-threaded; ratio stanli / CmdStan, median [min, max] over seeds)")
show(c("build_s_cmdstanr", "build_s_stanli", "wall_s_cmdstanr", "wall_s_stanli", "wall_ratio",
       "ms_per_leapfrog_ratio", "ess_per_wall_s_ratio"),
     "Fit, 4 chains in parallel (wall = sampling to brmsfit; ms per leapfrog in the sampling phase; ratios stanli / cmdstanr within seed)")
show(c("ess_1k_lf_cmdstanr", "ess_1k_lf_stanli", "leapfrog_ratio", "rhat_max", "divergent", "mean_z_max"),
     "Sampling: min bulk ESS (population level) per 1000 leapfrog steps; divergent cmdstanr/stanli; mean_z_max = largest |mean difference| / joint MCSE")
show(c("lp_identical", "lp_ulp_max", "grad_relerr_max", "nonfinite"),
     "Exactness at the posterior draws (stanli against CmdStan)")
if (nzchar(opt$out)) utils::write.csv(tab, opt$out, row.names = FALSE)
