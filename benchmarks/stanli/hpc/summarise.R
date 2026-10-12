# Aggregate what `run.sh pull` brought back (benchmarks/stanli/results/hpc/).
#
#   Rscript benchmarks/stanli/hpc/summarise.R [--dir benchmarks/stanli/results/hpc]
#                                             [--kind fit|grad] [--tag 0.20.0]
#
# With no option: the untagged 0.19.0 runs, fit_<k>/ and grad_<k>/. --tag
# reads a tagged rerun instead (STANLI_TAG on run.sh submit: <kind>_<tag>_<k>/,
# files *_<tag>.csv in bench.R's current, long format), the grad rerun
# unless --kind fit.
#
# Fit: one row per arm over seeds. Wall time on a shared node is noise
# (3x between identical fits in the September ablation), so the timings are
# compared within a seed - the four arms ran one after another on the same
# node - and the hardware-free figure is ESS per 1000 leapfrog steps.
# Grad: per task (one node each), the median per-gradient cost of each
# program and its ratio to CmdStan `orig`, then the exactness checks.
source("benchmarks/gradient_programs.R")  # gp_args()
opt <- gp_args(list(dir = "benchmarks/stanli/results/hpc", kind = "", tag = ""))
rd <- function(pattern) {
  f <- Sys.glob(file.path(opt$dir, pattern))
  do.call(rbind, lapply(f, function(x) cbind(task = basename(dirname(x)), utils::read.csv(x))))
}
# fit_<k>/fit_summary.csv untagged; fit_<tag>_<k>/fit_summary_<tag>.csv tagged.
# The untagged glob fit_*/ would also catch the tagged directories.
fit_tag <- if (nzchar(opt$tag) && opt$kind == "fit") opt$tag else ""
rd_fit <- function(stem) {
  if (nzchar(fit_tag)) return(rd(sprintf("fit_%s_*/%s_%s.csv", fit_tag, stem, fit_tag)))
  x <- rd(sprintf("fit_*/%s.csv", stem))
  if (!is.null(x)) x[grepl("^fit_[0-9]+$", x$task), ]
}
node <- function(task) vapply(task, function(t) {
  x <- readLines(file.path(opt$dir, t, "node.txt"))
  paste(sub("^node: ", "", x[1]), sub(".*: ", "", grep("model name", x, value = TRUE)))
}, character(1))
q <- function(x) sprintf("%.2f [%.2f, %.2f]", stats::median(x), stats::quantile(x, 0.1), stats::quantile(x, 0.9))

# ---- tagged grad rerun -------------------------------------------------------
# Per task, the median us per gradient and the ratio to CmdStan `orig` taken
# within each block (median [10th, 90th percentile] over blocks); then the
# pairs that split the variants into sigmas and likelihood; then exactness.
if (nzchar(opt$tag) && opt$kind != "fit") {
  tg <- opt$tag
  tm <- rd(sprintf("grad_%s_*/time_%s.csv", tg, tg))
  if (is.null(tm)) stop("no grad_", tg, "_*/time_", tg, ".csv under ", opt$dir, call. = FALSE)
  cat(sprintf("grad %s: %d tasks\n", tg, length(unique(tm$task))))
  print(data.frame(task = unique(tm$task), node = node(unique(tm$task))), row.names = FALSE)
  pair <- function(x, a, b) {
    w <- stats::reshape(x[c("rep", "program", "us_per_grad")], idvar = "rep",
                        timevar = "program", direction = "wide")
    names(w) <- sub("us_per_grad.", "", names(w), fixed = TRUE)
    if (!all(c(a, b) %in% names(w))) return(NA_character_)
    q(w[[a]] / w[[b]])
  }
  progs <- unique(tm$program)
  for (t in unique(tm$task)) {
    x <- tm[tm$task == t, ]
    med <- tapply(x$us_per_grad, x$program, stats::median)
    cat("\n", t, "\n", sep = "")
    tab <- data.frame(program = names(med), us = round(med),
                      vs_cmdstan_orig = vapply(names(med), function(p) pair(x, p, "cmdstan_orig"), ""))
    print(tab[order(tab$us), ], row.names = FALSE, right = FALSE)
  }
  pairs <- list(c("_orig_s", "_orig", "sigmas, brms likelihood"),
                c("_rw", "_rw_v", "sigmas, rewritten likelihood"),
                c("_rw_v", "_orig", "likelihood, vector sigmas"),
                c("_rw", "_orig_s", "likelihood, scalar sigmas"),
                c("_rw", "_orig", "both"))
  cat("\nwithin-block pairs, median [10th, 90th] per task:\n")
  for (e in c("cmdstan", "stanli")) for (p in pairs) {
    v <- vapply(unique(tm$task), function(t) pair(tm[tm$task == t, ], paste0(e, p[1]), paste0(e, p[2])), "")
    cat(sprintf("%-8s %-30s %s\n", e, p[3], paste(v, collapse = "  ")))
  }
  for (p in list(c("stanli_orig", "cmdstan_orig"), c("stanli_orig_fast", "cmdstan_orig"),
                 c("stanli_orig_fast", "stanli_orig"), c("stanli_orig_s_fast", "stanli_orig_s")))
    if (all(p %in% progs))
      cat(sprintf("%-8s %-30s %s\n", "", paste(p, collapse = " / "),
                  paste(vapply(unique(tm$task), function(t) pair(tm[tm$task == t, ], p[1], p[2]), ""),
                        collapse = "  ")))
  ck <- rd(sprintf("grad_%s_*/check_%s.csv", tg, tg))
  edges <- !grepl("^(init|near|wide)", ck$point)
  cat("\nexactness against CmdStan `orig`, all tasks: interior points\n")
  print(do.call(rbind, lapply(split(ck[!edges, ], ck$engine[!edges]), function(x) data.frame(
    lp_identical = sprintf("%d/%d", sum(x$lp_ulp == 0, na.rm = TRUE), nrow(x)),
    lp_ulp_max = max(x$lp_ulp, na.rm = TRUE),
    grad_relerr_max = signif(max(x$grad_relerr, na.rm = TRUE), 2),
    grad_nonfinite = sum(!x$grad_finite)))))
  cat("\nedge points where an engine's lp differs from CmdStan `orig`, or its gradient's finiteness does:\n")
  ref <- ck[ck$engine == "cmdstan_orig", c("task", "point", "grad_finite")]
  m <- merge(ck[edges, ], ref, by = c("task", "point"), suffixes = c("", "_orig"))
  print(unique(m[m$lp_ulp != 0 | m$grad_finite != m$grad_finite_orig,
                 c("point", "engine", "lp_ulp", "grad_finite", "grad_finite_orig")]), row.names = FALSE)
  tl <- rd(sprintf("grad_%s_*/tail_%s.csv", tg, tg))
  cat("\ntail (12 s trial, loser's x ~ -65):\n")
  print(stats::aggregate(cbind(lp_ulp, grad_finite, grad_relerr) ~ engine, tl, max, na.action = stats::na.pass))
  # bisect.R took --tag from 0.20.0 on; the 0.19.1 tasks wrote bisect.csv.
  bs <- rd(sprintf("grad_%s_*/bisect_%s.csv", tg, tg))
  if (is.null(bs)) bs <- rd(sprintf("grad_%s_*/bisect.csv", tg))
  pin <- Sys.glob(file.path(opt$dir, sprintf("grad_%s_*/repro_parallel_init.txt", tg)))
  if (length(pin)) {
    cat("\nrepro_parallel_init.R (wall/sum near 0.25-0.33: the chains overlapped; 1.0: serial):\n")
    for (f in pin) cat(basename(dirname(f)), grep("^init", readLines(f), value = TRUE), sep = "\n  ")
  }
  if (!is.null(bs)) {
    cat("\nbisect: median us per gradient, per task\n")
    print(stats::reshape(bs[c("task", "variant", "us_per_grad")], idvar = "variant",
                         timevar = "task", direction = "wide"), row.names = FALSE)
  }
  quit(save = "no")
}

# ---- fit ---------------------------------------------------------------------
fs <- rd_fit("fit_summary")
if (!is.null(fs)) {
  if (is.null(fs$ess_bulk_per_wall_s)) fs$ess_bulk_per_wall_s <- fs$ess_bulk_min / fs$wall_s
  ref <- fs[fs$arm == "cmdstanr_orig", c("seed", "ms_per_leapfrog", "ess_bulk_per_s", "ess_bulk_per_wall_s", "wall_s")]
  names(ref)[-1] <- paste0(names(ref)[-1], "_ref")
  fs <- merge(fs, ref, by = "seed")
  # The 2x2 first, in that order, then any other arm.
  known <- c("cmdstanr_orig", "cmdstanr_sel", "stanli_orig", "stanli_sel")
  arms <- c(intersect(known, fs$arm), sort(setdiff(unique(fs$arm), known)))
  cat(sprintf("fit: %d seeds; nodes: %s\n\n", length(unique(fs$seed)),
              paste(names(table(node(unique(fs$task)))), collapse = "; ")))
  tab <- do.call(rbind, lapply(arms, function(a) {
    x <- fs[fs$arm == a, ]
    data.frame(arm = a,
               ms_per_leapfrog = q(x$ms_per_leapfrog),
               ms_per_lf_vs_cmdstan = q(x$ms_per_leapfrog / x$ms_per_leapfrog_ref),
               ess_per_s_vs_cmdstan = q(x$ess_bulk_per_s / x$ess_bulk_per_s_ref),
               wall_s = q(x$wall_s),
               ess_per_wall_s_vs_cmdstan = q(x$ess_bulk_per_wall_s / x$ess_bulk_per_wall_s_ref),
               ess_per_1k_leapfrog = q(x$ess_bulk_per_1k_leapfrog),
               leapfrog = q(x$leapfrog),
               rhat_max = sprintf("%.3f", max(x$rhat_max)),
               divergent = sum(x$divergent))
  }))
  print(tab, row.names = FALSE, right = FALSE)
  cat("\nmedian [10th, 90th percentile] over seeds; ratios within seed\n")

  # Same posterior? Each arm's means against cmdstanr_orig's, in posterior SDs.
  fp <- rd_fit("fit_pars")
  r0 <- fp[fp$arm == "cmdstanr_orig", c("seed", "variable", "mean", "sd")]
  m <- merge(fp, r0, by = c("seed", "variable"), suffixes = c("", "_ref"))
  m$z <- abs(m$mean - m$mean_ref) / m$sd_ref
  cat("\nlargest |mean - cmdstanr_orig mean| / sd, per arm:\n")
  print(round(tapply(m$z, m$arm, max), 3))
}

# ---- grad --------------------------------------------------------------------
# The untagged 0.19.0 tasks only (grad_<k>/); tagged reruns are read above.
rd0 <- function(pattern) { x <- rd(pattern); if (!is.null(x)) x[grepl("^grad_[0-9]+$", x$task), ] }
tm <- if (!nzchar(opt$tag)) rd0("grad_*/time.csv")
if (!is.null(tm)) {
  med <- stats::aggregate(us_per_grad ~ task + program, tm, stats::median)
  w <- stats::reshape(med, idvar = "task", timevar = "program", direction = "wide")
  names(w) <- sub("us_per_grad.", "", names(w), fixed = TRUE)
  progs <- setdiff(names(w), "task")
  cat("\ngrad: median us per gradient, then / cmdstan_orig, per task\n")
  print(cbind(w["task"], node = node(w$task), round(w[progs])), row.names = FALSE)
  print(cbind(w["task"], round(w[progs] / w$cmdstan_orig, 2)), row.names = FALSE)
  ck <- rd0("grad_*/check.csv")
  g <- grep("^grad_", names(ck), value = TRUE)
  cat("\nlargest gradient relative error against CmdStan orig, over tasks and points:\n")
  print(signif(sapply(ck[g], max), 2))
  tl <- rd0("grad_*/tail.csv")
  cat("\ntail (12 s trial, loser's x ~ -60):\n")
  print(stats::aggregate(cbind(grad_finite, grad_relerr) ~ variant, tl, max, na.action = stats::na.pass))
  bs <- rd0("grad_*/bisect.csv")
  if (!is.null(bs)) {
    cat("\nbisect: median us per gradient, per task\n")
    print(stats::reshape(bs[c("task", "variant", "us_per_grad")], idvar = "variant",
                         timevar = "task", direction = "wide"), row.names = FALSE)
  }
}
