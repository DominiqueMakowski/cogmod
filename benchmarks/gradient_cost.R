# Gradient cost: did a change make a family's gradient dearer?
#
# Sampling time is gradient evaluations times the cost of one, and only the
# second is a property of the Stan code. So that is what is measured here: the
# wall time of grad_log_prob() on a fixed, simulated data set, for the same
# family's program generated from two working trees - the base of a pull
# request and its head - compiled side by side and timed in alternating
# blocks on the same machine, so that whatever else the machine is doing hits
# both alike. The headline is the ratio of medians; absolute times on a shared
# runner mean little, and only the ratio is worth reading.
#
# How much to trust that ratio: the alternating blocks handle whatever the
# machine is doing *within* a run, but nothing here controls for what it is
# doing between two runs. At the default --reps 7, measured 2026-09-18 on a
# Windows laptop, the same ratio came back 0.15 apart across runs, and
# byte-identical code read 1.02 once and 1.19 another time. Treat 7 blocks as
# a screen: it will catch a real regression, which reproduced at 1.47 and 1.48
# there, but it will not tell 1.1 from 1.25. Before quoting a number, or
# before concluding anything about a ratio near --fail, re-run it with
# --reps 21; the same comparison settled to 1.11 and 1.12 at that setting.
#
# A family whose generated Stan program is byte-for-byte the same in both
# trees is reported as unchanged and neither compiled nor timed: that is
# where nearly every PR lands, and it keeps the job short.
#
# Two steps, because two versions of a package cannot be loaded into one R
# session. Each tree emits its programs in a process of its own, then a third
# process - which needs cmdstanr and nothing of cogmod - times them:
#
#   Rscript benchmarks/gradient_cost.R emit --pkg ../base --out /tmp/base
#   Rscript benchmarks/gradient_cost.R emit --pkg .       --out /tmp/pr
#   Rscript benchmarks/gradient_cost.R time --base /tmp/base --pr /tmp/pr
#
#   emit:  --pkg (tree to load), --out (directory), --families (subset),
#          --n (trials, default 5000)
#   time:  --base, --pr (the two emitted directories), --out (report directory,
#          default benchmarks/results/gradient_cost), --reps (alternating
#          blocks, 7), --block (gradients per block, 25), --fail (ratio above
#          which the script exits 1; 1.3), --warn (ratio to flag; 1.1)
#
# The `time` step writes <out>/gradient_cost.csv and <out>/summary.md, the
# latter shaped for a pull-request comment.

source("benchmarks/gradient_programs.R")

mode <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(mode) || !mode %in% c("emit", "time")) {
  stop("usage: Rscript benchmarks/gradient_cost.R emit|time [--options]", call. = FALSE)
}
# gp_args() reads all trailing args; drop the mode first.
commandArgs <- local({
  orig <- base::commandArgs
  function(trailingOnly = FALSE) { a <- orig(trailingOnly); if (trailingOnly) a[-1] else a }
})

# ---- emit --------------------------------------------------------------------
if (mode == "emit") {
  opt <- gp_args(list(pkg = ".", out = "", families = "", n = 5000L))
  if (!nzchar(opt$out)) stop("--out is required", call. = FALSE)
  gp_load(opt$pkg)
  fams <- if (nzchar(opt$families)) strsplit(opt$families, ",")[[1]] else gp_families()
  ver <- as.character(utils::packageVersion("cogmod"))
  sha <- tryCatch(system2("git", c("-C", opt$pkg, "rev-parse", "--short", "HEAD"), stdout = TRUE),
                  error = function(e) NA_character_)
  cat(sprintf("Emitting %d programs from %s (cogmod %s, %s)\n", length(fams), opt$pkg, ver, sha))
  gp_emit(opt$out, fams, n = opt$n, extremes = FALSE)
  writeLines(c(ver, sha), file.path(opt$out, "VERSION"))
  quit(status = 0)
}

# ---- time --------------------------------------------------------------------
opt <- gp_args(list(base = "", pr = "", out = "benchmarks/results/gradient_cost",
                    reps = 7L, block = 25L, fail = 1.3, warn = 1.1))
if (!nzchar(opt$base) || !nzchar(opt$pr)) stop("--base and --pr are required", call. = FALSE)
suppressPackageStartupMessages(library(cmdstanr))
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)

label <- function(dir) {
  v <- file.path(dir, "VERSION")
  if (file.exists(v)) paste(readLines(v), collapse = " @ ") else basename(dir)
}
pr_stan <- list.files(opt$pr, pattern = "\\.stan$", full.names = TRUE)
fams <- sub("\\.stan$", "", basename(pr_stan))

# One timing: `block` gradients at the start point, in seconds.
time_block <- function(fit, up, block) {
  t0 <- Sys.time()
  for (i in seq_len(block)) fit$grad_log_prob(up, jacobian = TRUE)
  as.numeric(Sys.time() - t0, units = "secs")
}

# Compile and stand up model methods for one emitted program, or NULL with
# the reason on failure.
stand_up <- function(dir, fam) {
  stan <- file.path(dir, paste0(fam, ".stan"))
  data <- file.path(dir, paste0(fam, ".data.json"))
  init <- file.path(dir, paste0(fam, ".init.json"))
  if (!file.exists(init)) init <- NA_character_
  tryCatch({
    mod <- gp_compile(stan, dir = dir)
    gp_methods(mod, data, init)
  }, error = function(e) structure(list(), reason = conditionMessage(e)))
}

rows <- list()
cat(sprintf("base: %s\npr:   %s\n", label(opt$base), label(opt$pr)))
for (fam in fams) {
  base_stan <- file.path(opt$base, paste0(fam, ".stan"))
  pr_code <- readLines(file.path(opt$pr, paste0(fam, ".stan")))
  if (!file.exists(base_stan)) {
    rows[[fam]] <- data.frame(family = fam, status = "new in PR", base_ms = NA, pr_ms = NA, ratio = NA)
    cat(sprintf("  %-22s new in PR\n", fam)); next
  }
  # Compared with comments and blank lines stripped: the preludes carry their
  # reasoning as Stan comments, and rewording one should not cost two
  # compilations and a timing.
  strip <- function(x) { x <- sub("//.*$", "", x); x <- trimws(x); x[nzchar(x)] }
  if (identical(strip(readLines(base_stan)), strip(pr_code))) {
    rows[[fam]] <- data.frame(family = fam, status = "unchanged", base_ms = NA, pr_ms = NA, ratio = NA)
    cat(sprintf("  %-22s unchanged\n", fam)); next
  }
  cat(sprintf("  %-22s changed; compiling both ...", fam))
  m <- list(base = stand_up(opt$base, fam), pr = stand_up(opt$pr, fam))
  bad <- vapply(m, function(x) length(x) == 0, logical(1))
  if (any(bad)) {
    why <- paste(sprintf("%s: %s", names(m)[bad], vapply(m[bad], attr, "", "reason")), collapse = "; ")
    rows[[fam]] <- data.frame(family = fam, status = paste("failed -", why), base_ms = NA, pr_ms = NA, ratio = NA)
    cat(" FAILED\n"); next
  }
  # Warm both up once, then alternate, base first on odd blocks and PR first on
  # even ones, so neither always follows the other.
  for (k in names(m)) time_block(m[[k]]$fit, m[[k]]$up0, 3)
  secs <- list(base = numeric(0), pr = numeric(0))
  for (r in seq_len(opt$reps)) {
    order <- if (r %% 2 == 1) c("base", "pr") else c("pr", "base")
    for (k in order) secs[[k]] <- c(secs[[k]], time_block(m[[k]]$fit, m[[k]]$up0, opt$block))
  }
  ms <- vapply(secs, function(s) 1000 * stats::median(s) / opt$block, numeric(1))
  ratio <- unname(ms["pr"] / ms["base"])
  status <- if (ratio > opt$fail) "SLOWER" else if (ratio > opt$warn) "slower" else if (ratio < 1 / opt$warn) "faster" else "same"
  rows[[fam]] <- data.frame(family = fam, status = status, base_ms = unname(ms["base"]),
                            pr_ms = unname(ms["pr"]), ratio = ratio)
  cat(sprintf(" base %.3f ms, PR %.3f ms per gradient, ratio %.2f (%s)\n", ms["base"], ms["pr"], ratio, status))
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(opt$out, "gradient_cost.csv"), row.names = FALSE)

timed <- res[!is.na(res$ratio), ]
worst <- if (nrow(timed)) max(timed$ratio) else NA
verdict <- if (any(grepl("^failed", res$status))) "FAILED" else
  if (isTRUE(worst > opt$fail)) sprintf("regression (worst ratio %.2f)", worst) else
  if (nrow(timed) == 0) "no Stan program changed" else sprintf("no regression (worst ratio %.2f)", worst)
# The number of trials as emitted (--n at the emit step, plus nothing: the cost
# programs carry no tail responses), read back from the data rather than
# assumed.
n_trials <- vapply(fams, function(fam) {
  d <- file.path(opt$pr, paste0(fam, ".data.json"))
  if (file.exists(d)) as.integer(jsonlite::fromJSON(d)$N) else NA_integer_
}, integer(1))
n_label <- if (length(unique(stats::na.omit(n_trials))) == 1) {
  sprintf("%d simulated trials", unique(stats::na.omit(n_trials)))
} else {
  sprintf("%d to %d simulated trials", min(n_trials, na.rm = TRUE), max(n_trials, na.rm = TRUE))
}
md <- c(
  sprintf("## Gradient cost: %s", verdict),
  "",
  sprintf("Base %s vs PR %s. Milliseconds per `grad_log_prob()` on %s, median of %d alternating blocks of %d; ratio = PR / base, flagged above %.2f, failing above %.2f.",
          label(opt$base), label(opt$pr), n_label, opt$reps, opt$block, opt$warn, opt$fail),
  ""
)
if (nrow(timed)) {
  md <- c(md,
    "| family | base ms | PR ms | ratio | |",
    "|---|---:|---:|---:|---|",
    sprintf("| %s | %.3f | %.3f | %.2f | %s |", timed$family, timed$base_ms, timed$pr_ms, timed$ratio, timed$status),
    "")
}
other <- res[is.na(res$ratio), ]
if (nrow(other)) {
  md <- c(md, sprintf("%s: %s.", sub(" -.*$", "", other$status), other$family))
  # collapse the unchanged ones into one line
  unch <- other$family[other$status == "unchanged"]
  md <- c(md[!grepl("^unchanged:", md)],
          if (length(unch)) sprintf("Unchanged Stan program (not timed): %s.", paste(unch, collapse = ", ")))
}
writeLines(md, file.path(opt$out, "summary.md"))
cat("\n", paste(md, collapse = "\n"), "\n", sep = "")

if (any(grepl("^failed", res$status)) || isTRUE(worst > opt$fail)) {
  cat("::error::gradient cost:", verdict, "\n")
  quit(status = 1)
}
