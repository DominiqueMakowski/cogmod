# Reports every number Example 2's clinical section quotes, from the objects
# example_cnp.R saved. Separate from the fitting script so the numbers can be
# re-read without a 28-minute refit.
#
#   Rscript report_cnp.R      (run from the `paper/` directory)

suppressMessages({
  library(dplyr); library(brms); library(cogmod); library(bayestestR)
})

DIR <- "cnp_data"
d  <- readRDS(file.path(DIR, "d.rds"))
sc <- readRDS(file.path(DIR, "sc.rds"))
m  <- readRDS(file.path(DIR, "m_ln.rds"))

cat(sprintf("%s trials, %d participants (%s)\n",
            format(nrow(d), big.mark = ","), n_distinct(d$Participant),
            paste(sprintf("%s = %d", levels(d$Diagnosis),
                          tapply(d$Participant, d$Diagnosis,
                                 function(x) length(unique(x)))),
                  collapse = ", ")))
cat("max trials per participant:", max(table(d$Participant)), "\n\n")

# ---- The conventional analysis ---------------------------------------------

cat("--- switch cost (the standard index) ---\n")
tt <- t.test(cost ~ Diagnosis, data = sc)
cat(sprintf("CONTROL %.1f ms, ADHD %.1f ms, t(%.0f) = %.2f, p = %.3f\n",
            1000 * mean(sc$cost[sc$Diagnosis == "CONTROL"]),
            1000 * mean(sc$cost[sc$Diagnosis == "ADHD"]),
            tt$parameter, tt$statistic, tt$p.value))

cat("\n--- every other summary statistic ---\n")
stats <- d |>
  group_by(Participant, Diagnosis) |>
  summarise(mean = mean(RT), sd = sd(RT), cv = sd(RT) / mean(RT),
            skew = mean((RT - mean(RT))^3) / sd(RT)^3,
            q10 = quantile(RT, .1), q90 = quantile(RT, .9), .groups = "drop")
for (v in c("mean", "sd", "cv", "skew", "q10", "q90")) {
  tt <- t.test(stats[[v]] ~ stats$Diagnosis)
  cat(sprintf("%-5s t = %6.2f, p = %.3f\n", v, tt$statistic, tt$p.value))
}

cat("\n--- linear mixed model ---\n")
m_lm <- lme4::lmer(RT ~ Diagnosis * Switching + (Switching | Participant), data = d)
print(parameters::parameters(m_lm))

# ---- The model -------------------------------------------------------------

p <- posterior::as_draws_df(m)
sp <- function(x) log1p(exp(x))
pd <- function(x) as.numeric(p_direction(x))

ndt_c <- exp(p$b_ndt_Intercept)
ndt_a <- exp(p$b_ndt_Intercept + p$b_ndt_DiagnosisADHD)
sig_c <- sp(p$b_sigma_Intercept)
sig_a <- sp(p$b_sigma_Intercept + p$b_sigma_DiagnosisADHD)
dec_c <- exp(p$b_Intercept + sig_c^2 / 2)
dec_a <- exp(p$b_Intercept + p$b_DiagnosisADHD + sig_a^2 / 2)

cat("\n--- response-scale group differences (ADHD - control) ---\n")
report <- function(x, label) {
  cat(sprintf("%-22s %7.1f ms [%6.1f, %6.1f]  pd = %.3f\n", label,
              1000 * median(x), 1000 * quantile(x, .025),
              1000 * quantile(x, .975), pd(x)))
}
report(dec_a - dec_c, "decision component")
report(ndt_a - ndt_c, "non-decision time")
report((ndt_a + dec_a) - (ndt_c + dec_c), "total mean RT")

cat("\n--- link-scale effects the text quotes ---\n")
cat(sprintf("sigma  Diagnosis            = %+.3f, pd = %.3f\n",
            median(p$b_sigma_DiagnosisADHD), pd(p$b_sigma_DiagnosisADHD)))
cat(sprintf("mu     Diagnosis x Switching = %+.3f, pd = %.3f\n",
            median(p$`b_DiagnosisADHD:SwitchingSWITCH`),
            pd(p$`b_DiagnosisADHD:SwitchingSWITCH`)))
cat(sprintf("sigma  Diagnosis x Switching = %+.3f, pd = %.3f\n",
            median(p$`b_sigma_DiagnosisADHD:SwitchingSWITCH`),
            pd(p$`b_sigma_DiagnosisADHD:SwitchingSWITCH`)))
cat(sprintf("posterior cor(mu_ADHD, ndt_ADHD) = %.3f\n",
            cor(p$b_DiagnosisADHD, p$b_ndt_DiagnosisADHD)))
