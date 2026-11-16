# What each between-trial variability costs per gradient, relative to the
# model with sigmadrift alone (the ratios ?cogmod_ddm quotes). One intercept-
# only program per configuration from the tree at --pkg, the same simulated
# data for all, timed in alternating blocks through CmdStan's model methods.
#
#   Rscript benchmarks/wiener_stan/ladder.R [--pkg .] [--n 300] [--reps 11] [--out DIR]
source("benchmarks/gradient_programs.R")
opt <- gp_args(list(pkg = ".", n = 300L, reps = 11L,
                    out = file.path(Sys.getenv("TEMP"), "wiener_stan_ladder")))
gp_load(opt$pkg)
dir.create(opt$out, recursive = TRUE, showWarnings = FALSE)
.GP_VARIANTS <- c(.GP_VARIANTS, list(
  ddm6_st0 = list(family = "cogmod_ddm", fixed = c(sigmabias = 0)),
  ddm6_sw = list(family = "cogmod_ddm", fixed = c(sigmandt = 0))))
keys <- c(ddm4 = "cogmod_ddm4", ddm5 = "cogmod_ddm5", ddm6_st0 = "ddm6_st0",
          ddm6_sw = "ddm6_sw", ddm7 = "cogmod_ddm")
em <- gp_emit(file.path(opt$out, "programs"), unname(keys), n = opt$n, extremes = FALSE)
fits <- lapply(seq_len(nrow(em)), function(i) {
  cat("compiling", em$family[i], "\n")
  gp_methods(gp_compile(em$stan[i]), em$data[i], em$init[i])
})
names(fits) <- names(keys)
for (k in names(fits)) {
  g <- fits[[k]]$fit$grad_log_prob(fits[[k]]$up0)
  cat(sprintf("%-9s lp %.6g, gradient finite %s\n", k, attr(g, "log_prob"), all(is.finite(g))))
}
blk <- function(k, m) { f <- fits[[k]]; for (i in seq_len(m)) f$fit$grad_log_prob(f$up0) }
# Windows' clock ticks at about 10 ms, so a block is sized to about 0.5 s by
# doubling until one takes at least 0.1 s.
passes <- vapply(names(fits), function(k) {
  m <- 1L
  repeat {
    el <- system.time(blk(k, m))[["elapsed"]]
    if (el >= 0.1) break
    m <- m * 2L
  }
  max(1L, as.integer(ceiling(m * 0.5 / el)))
}, integer(1))
tm <- list()
for (r in seq_len(opt$reps)) for (k in sample(names(fits))) {
  st <- system.time(blk(k, passes[[k]]))[["elapsed"]]
  tm[[length(tm) + 1]] <- data.frame(rep = r, program = k, ms = st / passes[[k]] * 1e3)
}
tm <- do.call(rbind, tm)
utils::write.csv(tm, file.path(opt$out, "ladder.csv"), row.names = FALSE)
b <- tm[tm$program == "ddm5", ]
cat(sprintf("\nms per gradient on %d trials, median of %d rounds; ratio to ddm5 within round, median [min, max]\n",
            opt$n, opt$reps))
for (k in names(fits)) {
  x <- tm[tm$program == k, ]
  rr <- x$ms[order(x$rep)] / b$ms[order(b$rep)]
  cat(sprintf("  %-9s %9.3f ms  %7.2f [%.2f, %.2f]\n", k, stats::median(x$ms), stats::median(rr), min(rr), max(rr)))
}
