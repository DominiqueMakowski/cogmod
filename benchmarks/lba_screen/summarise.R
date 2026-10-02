# Summarise the LBA screen: per arm, do the eight chains agree?
#
#   Rscript summarise.R [results dir]
#
# Per arm: lp__ of each chain relative to the arm's median chain (a chain
# more than 100 units away is in another place), frozen chains (step size
# below a hundredth of the arm's median), the pooled Rhat over the
# population-level parameters with the five worst, warmup leapfrogs and
# divergences. Then the decision rules of the README, applied mechanically,
# as a starting point for reading the tables rather than a verdict.
suppressPackageStartupMessages(library(posterior))
options(width = 160)
args <- commandArgs(trailingOnly = TRUE)
dir <- if (length(args)) args[1] else "results"
source("cells.R")

files <- list.files(dir, pattern = "\\.rds$", full.names = TRUE)
if (!length(files)) stop("no .rds results in ", dir)
res <- lapply(files, readRDS)
cells <- do.call(rbind, lapply(res, `[[`, "cell"))
arms <- screen_arms()
arms <- arms[arms$arm %in% cells$arm, ]

verdict <- list()
for (a in arms$arm) {
  k <- which(cells$arm == a)
  k <- k[order(cells$chain[k])]
  cc <- cells[k, ]
  cat(sprintf("\n== %s (%s, %s defaults): %d chain(s) done\n", a,
              arms$model[arms$arm == a], arms$defaults[arms$arm == a], length(k)))
  dlp <- cc$lp_mean - stats::median(cc$lp_mean)
  frozen <- cc$step_size < stats::median(cc$step_size) / 100
  tab <- data.frame(chain = cc$chain, d_lp = round(dlp), lp_sd = cc$lp_sd,
                    step = cc$step_size, frozen = frozen,
                    warmup_lf = cc$warmup_leapfrog, maxtree = cc$sample_maxtree_frac,
                    div = cc$divergences, hours = cc$wall_h)
  print(tab, row.names = FALSE)

  rh <- NA_real_
  if (length(k) >= 2) {
    dr <- do.call(bind_draws, c(lapply(res[k], `[[`, "draws"), along = "chain"))
    sm <- summarise_draws(subset_draws(dr, variable = "lp__", exclude = TRUE), "rhat")
    sm <- sm[order(-sm$rhat), ]
    rh <- max(sm$rhat, na.rm = TRUE)
    cat(sprintf("pooled max Rhat %.2f over %d population-level parameters; worst:\n",
                rh, nrow(sm)))
    print(as.data.frame(utils::head(sm, 5)), row.names = FALSE, digits = 3)
  }
  verdict[[a]] <- data.frame(
    arm = a, chains = length(k), off = sum(abs(dlp) > 100), frozen = sum(frozen),
    max_rhat = round(rh, 2),
    agrees = length(k) >= 8 && all(abs(dlp) <= 100) && !any(frozen) && rh < 1.05
  )
}

v <- do.call(rbind, verdict)
cat("\n== Summary (agrees = 8 chains, all within 100 lp__ of the median, none",
    "frozen, pooled Rhat < 1.05)\n")
print(v, row.names = FALSE)

ok <- stats::setNames(v$agrees, v$arm)
say <- function(...) cat("->", ..., "\n")
cat("\n== Decision rules (README)\n")
if (all(c("lba_old", "lba_new", "lba_a", "lba_b") %in% names(ok))) {
  if (ok[["lba_old"]]) {
    say("lba_old agrees: 200 participants cannot see the split. Step up to 480 before",
        "spending full-data time.")
  } else if (ok[["lba_new"]]) {
    say("lba_new agrees where lba_old splits: the defaults mattered. Re-run full data",
        "unchanged with 0.3.4.")
  } else if (ok[["lba_b"]] && !ok[["lba_a"]]) {
    say("Start-point range: refit full data as (b), and screen gam_rdm5 / gam_lnr6 the same way.")
  } else if (ok[["lba_a"]] && !ok[["lba_b"]]) {
    say("The ray: option (a), sigmaone ~ 1 + (1 | Participant).")
  } else if (ok[["lba_a"]] && ok[["lba_b"]]) {
    say("Both variants agree: choose on elpd at 200 and on what needs interpreting.")
  } else {
    say("Neither variant agrees: try (a) and (b) together, then sigmaone = 1, then drop the LBA.")
  }
} else {
  say("LBA arms incomplete; no rule applied.")
}
if (all(c("rdm5_new", "rdm5_b") %in% names(ok))) {
  say(sprintf("RDM: rdm5_new %s, rdm5_b %s.",
              if (ok[["rdm5_new"]]) "agrees" else "splits",
              if (ok[["rdm5_b"]]) "agrees" else "splits"))
}
