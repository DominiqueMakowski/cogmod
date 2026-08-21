# Regenerates the two LOO tables of the paper (@tbl-eam and @tbl-loo) and the
# handful of parameter estimates quoted in Examples 1 and 2, from the model
# files fitted by the package vignettes.
#
#   Rscript make_tables.R      (run from the `paper/` directory)
#
# Nothing is refitted here; the vignettes own those fits.

library(cogmod)
library(brms)

VIG <- file.path("..", "vignettes", "models")
rds <- function(f) readRDS(file.path(VIG, paste0(f, ".rds")))

# ---------------------------------------------------------------------------
# Shared table builder: LOOIC, ENP, ELPD, difference against the best model,
# SE of that difference, and median per-chain sampling time.
# ---------------------------------------------------------------------------

loo_table <- function(models) {
  elpd <- vapply(models, function(m) m$criteria$loo$estimates["elpd_loo", "Estimate"], numeric(1))
  p_loo <- vapply(models, function(m) m$criteria$loo$estimates["p_loo", "Estimate"], numeric(1))
  looic <- vapply(models, function(m) m$criteria$loo$estimates["looic", "Estimate"], numeric(1))
  best <- names(which.max(elpd))
  a <- models[[best]]$criteria$loo$pointwise[, "elpd_loo"]

  se <- vapply(names(models), function(nm) {
    if (nm == best) return(0)
    b <- models[[nm]]$criteria$loo$pointwise[, "elpd_loo"]
    sqrt(length(a)) * stats::sd(b - a)
  }, numeric(1))

  time <- vapply(models, function(m) {
    t <- try(attributes(m$fit)$metadata$time$chain$total, silent = TRUE)
    if (inherits(t, "try-error") || is.null(t)) NA_real_ else stats::median(t)
  }, numeric(1))

  # Worst Pareto k, so the caption can say whether the differences are safe
  khat <- vapply(models, function(m) {
    k <- try(m$criteria$loo$diagnostics$pareto_k, silent = TRUE)
    if (inherits(k, "try-error") || is.null(k)) NA_real_ else max(k)
  }, numeric(1))

  out <- data.frame(
    Model = names(models), LOOIC = looic, ENP = p_loo, ELPD = elpd,
    Difference = elpd - max(elpd), SE = se, Time = time, max_k = khat,
    row.names = NULL
  )
  out[order(-out$ELPD), ]
}

md <- function(tab) {
  cat("\n| Model | LOOIC | ENP | ELPD | Difference | SE | Time (s) |\n")
  cat("|:------|------:|----:|-----:|-----------:|---:|---------:|\n")
  for (i in seq_len(nrow(tab))) {
    r <- tab[i, ]
    cat(sprintf("| %s | %s | %.2f | %.1f | %s | %.1f | %.1f |\n",
                r$Model,
                formatC(r$LOOIC, format = "f", digits = 1, big.mark = ","),
                r$ENP, r$ELPD,
                if (r$Difference == 0) "0.0" else sprintf("%.1f", r$Difference),
                r$SE, r$Time))
  }
  cat(sprintf("\nworst Pareto k by model:\n"))
  print(round(setNames(tab$max_k, tab$Model), 2))
}

# ---------------------------------------------------------------------------
# Table 1 (@tbl-eam): five choice+RT architectures
# ---------------------------------------------------------------------------

eam <- list(
  DDM = rds("m_ddm"), `DDM-5` = rds("m_ddm5"), LNR = rds("m_lnr"),
  LBA = rds("m_lba2"), RDM = rds("m_rdm")
)
cat("=== @tbl-eam: n =", nrow(eam$DDM$data), "===\n")
t_eam <- loo_table(eam)
md(t_eam)

cat("\n--- DDM condition effects (linear predictor) ---\n")
print(round(brms::fixef(eam$DDM)[, c("Estimate", "Q2.5", "Q97.5")], 3))
cat("\n--- DDM-5 sigmadrift ---\n")
print(round(brms::fixef(eam$`DDM-5`)[, c("Estimate", "Q2.5", "Q97.5")], 3))

# ---------------------------------------------------------------------------
# Table 2 (@tbl-loo): RT-only families
# ---------------------------------------------------------------------------

rt <- list(
  Gaussian = rds("m_normal"),
  ExGaussian = rds("m_exgauss"),
  `Shifted LogNormal` = rds("m_lognormal"),
  `Shifted Wald` = rds("m_wald"),
  `Wald + drift variability` = rds("m_wald4"),
  `Single-accumulator LBA` = rds("m_lba1"),
  `Recinormal (LATER)` = rds("m_recinormal"),
  `Ex-Wald` = rds("m_exwald"),
  `Birnbaum-Saunders` = rds("m_bisa"),
  `Log-Student` = rds("m_logstudent"),
  LogGamma = rds("m_loggamma"),
  Weibull = rds("m_weibull"),
  LogWeibull = rds("m_logweibull"),
  `Inverse Weibull` = rds("m_invweibull"),
  Gamma = rds("m_gamma"),
  `Inverse Gamma` = rds("m_invgamma")
)
cat("\n\n=== @tbl-loo: n =", nrow(rt$Gaussian$data), "===\n")
t_rt <- loo_table(rt)
md(t_rt)

# ---------------------------------------------------------------------------
# Estimates quoted in the Example 2 prose
# ---------------------------------------------------------------------------

cat("\n--- Gaussian: condition effect (s) ---\n")
print(round(brms::fixef(rt$Gaussian)[, c("Estimate", "Q2.5", "Q97.5")], 4))

cat("\n--- ExGaussian: mu / sigma / tau by condition (s) ---\n")
d <- as_draws_df(rt$ExGaussian)
sp <- function(x) log1p(exp(x))
eg <- list(
  mu    = list(A = d$b_Intercept,           B = d$b_Intercept + d$b_ConditionSpeed),
  sigma = list(A = sp(d$b_sigma_Intercept), B = sp(d$b_sigma_Intercept + d$b_sigma_ConditionSpeed)),
  tau   = list(A = sp(d$b_tau_Intercept),   B = sp(d$b_tau_Intercept + d$b_tau_ConditionSpeed))
)
for (nm in names(eg)) {
  cat(sprintf("%-6s Accuracy = %.3f [%.3f, %.3f]   Speed = %.3f [%.3f, %.3f]\n", nm,
              mean(eg[[nm]]$A), quantile(eg[[nm]]$A, .025), quantile(eg[[nm]]$A, .975),
              mean(eg[[nm]]$B), quantile(eg[[nm]]$B, .025), quantile(eg[[nm]]$B, .975)))
}
mean_A <- eg$mu$A + eg$tau$A
mean_B <- eg$mu$B + eg$tau$B
cat(sprintf("implied mean: Accuracy = %.3f s, Speed = %.3f s, diff = %.3f s\n",
            mean(mean_A), mean(mean_B), mean(mean_B - mean_A)))
cat(sprintf("share of the speed-up carried by tau: %.0f%%   by mu: %.0f%%\n",
            100 * mean(eg$tau$A - eg$tau$B) / mean(mean_A - mean_B),
            100 * mean(eg$mu$A - eg$mu$B) / mean(mean_A - mean_B)))

cat("\n--- Shifted Wald: drift and boundary by condition ---\n")
print(round(brms::fixef(rt$`Shifted Wald`)[, c("Estimate", "Q2.5", "Q97.5")], 3))
dw <- as_draws_df(rt$`Shifted Wald`)
cat(sprintf("drift: Accuracy = %.2f, Speed = %.2f\n",
            mean(dw$b_Intercept), mean(dw$b_Intercept + dw$b_ConditionSpeed)))
