# =============================================================================
# Example 2 of the paper, clinical reanalysis section.
#
# Task-switching data from the UCLA Consortium for Neuropsychiatric Phenomics
# (Poldrack et al., 2016), OpenNeuro ds000030. The data are NOT redistributed
# with this repository; run `download_cnp.R` first to fetch the trial-level
# event files from the OpenNeuro S3 bucket into cnp_data/.
#
# Runtime: the brms fit took ~27 min on 16 cores.
# =============================================================================

library(dplyr)
library(brms)
library(cogmod)
library(bayestestR)

DIR <- "cnp_data"

# ---- 1. Load the downloaded files --------------------------------------------

pfile <- file.path(DIR, "participants.tsv")
if (!file.exists(pfile)) stop("Run download_cnp.R first.")
participants <- read.delim(pfile)

subs <- participants |>
  filter(diagnosis %in% c("CONTROL", "ADHD"), taskswitch == "1") |>
  pull(participant_id)

# ---- 2. Assemble trial-level data --------------------------------------------

# Some event files code missing values as "n/a", which makes read.delim() give
# the column back as character for those participants and only those; declare
# the string so every file yields the same types.
d <- lapply(subs, function(s) {
  f <- file.path(DIR, "events", paste0(s, ".tsv"))
  if (!file.exists(f)) return(NULL)
  e <- read.delim(f, na.strings = c("n/a", "NA", ""))
  data.frame(
    Participant = s,
    Diagnosis   = participants$diagnosis[participants$participant_id == s],
    RT          = as.numeric(e$ReactionTime),
    Correct     = as.numeric(e$CorrectResp),
    Switching   = as.character(e$Switching),
    Congruency  = as.character(e$Congruency)
  )
}) |> bind_rows()

# Correct responses only; standard RT trimming. RT == 0 codes a non-response.
d <- d |>
  filter(!is.na(RT), !is.na(Correct), Correct == 1, RT > 0.2, RT < 2) |>
  mutate(
    Diagnosis = factor(Diagnosis, levels = c("CONTROL", "ADHD")),
    Switching = factor(Switching, levels = c("NOSWITCH", "SWITCH"))
  )

# keep participants contributing trials in both switching conditions
keep <- d |> count(Participant, Switching) |> count(Participant) |>
  filter(n == 2) |> pull(Participant)
d <- filter(d, Participant %in% keep)

cat(sprintf("%d trials, %d participants\n", nrow(d), n_distinct(d$Participant)))

# ---- 3. The conventional analysis --------------------------------------------

# Switch cost: per-participant difference in mean RT
sc <- d |>
  group_by(Participant, Diagnosis, Switching) |>
  summarise(m = mean(RT), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Switching, values_from = m) |>
  mutate(cost = SWITCH - NOSWITCH)

print(t.test(cost ~ Diagnosis, data = sc))

# Every other summary statistic a neuropsychologist might compute
stats <- d |>
  group_by(Participant, Diagnosis) |>
  summarise(mean = mean(RT), sd = sd(RT), cv = sd(RT) / mean(RT),
            skew = mean((RT - mean(RT))^3) / sd(RT)^3,
            q10 = quantile(RT, .1), q90 = quantile(RT, .9), .groups = "drop")

for (v in c("mean", "sd", "cv", "skew", "q10", "q90")) {
  tt <- t.test(stats[[v]] ~ stats$Diagnosis)
  cat(sprintf("%-5s t = %6.2f, p = %.4f\n", v, tt$statistic, tt$p.value))
}

# Linear mixed model
m_lm <- lme4::lmer(RT ~ Diagnosis * Switching + (Switching | Participant), data = d)
print(parameters::parameters(m_lm))

# ---- 4. The cogmod model -----------------------------------------------------

f <- bf(
  RT ~ Diagnosis * Switching + (Switching | Participant),
  sigma ~ Diagnosis * Switching + (1 | Participant),
  ndt ~ Diagnosis + (1 | Participant),
  poutlier ~ 1,
  family = cogmod_lognormal()
)

m_ln <- brm(
  f,
  data = d,
  prior = cogmod_priors(f, d),
  init = cogmod_inits(f, d),
  stanvars = cogmod_stanvars(f),
  chains = 4, cores = 4, iter = 1000, warmup = 500,
  threads = threading(3),
  backend = "cmdstanr", seed = 42
)

print(summary(m_ln))

# ---- 5. Response-scale decomposition -----------------------------------------

p        <- posterior::as_draws_df(m_ln)
softplus <- function(x) log1p(exp(x))

# `ndt` is on a `log` link and is in seconds outright. It used to be called
# `tau` and to be a logit-scaled fraction of a `minrt` argument; both are gone.
ndt_c <- exp(p$b_ndt_Intercept)
ndt_a <- exp(p$b_ndt_Intercept + p$b_ndt_DiagnosisADHD)
sig_c <- softplus(p$b_sigma_Intercept)
sig_a <- softplus(p$b_sigma_Intercept + p$b_sigma_DiagnosisADHD)
dec_c <- exp(p$b_Intercept + sig_c^2 / 2)
dec_a <- exp(p$b_Intercept + p$b_DiagnosisADHD + sig_a^2 / 2)

report <- function(x, label) {
  cat(sprintf("%-22s %7.1f ms [%6.1f, %6.1f]  pd = %.3f\n", label,
              1000 * median(x), 1000 * quantile(x, .025), 1000 * quantile(x, .975),
              as.numeric(p_direction(x))))
}
report(dec_a - dec_c, "decision component")
report(ndt_a - ndt_c, "non-decision time")
report((ndt_a + dec_a) - (ndt_c + dec_c), "total mean RT")

# The dispersion effect is reported on the link scale in the text, so give its
# pd directly, along with the two interaction terms the model agrees are null.
pd <- function(x) as.numeric(p_direction(x))
cat(sprintf("sigma: Diagnosis effect (link) = %+.3f, pd = %.3f\n",
            median(p$b_sigma_DiagnosisADHD), pd(p$b_sigma_DiagnosisADHD)))
cat(sprintf("mu:    Diagnosis x Switching   = %+.3f, pd = %.3f\n",
            median(p$`b_DiagnosisADHD:SwitchingSWITCH`),
            pd(p$`b_DiagnosisADHD:SwitchingSWITCH`)))
cat(sprintf("sigma: Diagnosis x Switching   = %+.3f, pd = %.3f\n",
            median(p$`b_sigma_DiagnosisADHD:SwitchingSWITCH`),
            pd(p$`b_sigma_DiagnosisADHD:SwitchingSWITCH`)))

# The decision and non-decision parameters trade off; report the correlation
cat("posterior cor(mu_ADHD, ndt_ADHD) =",
    round(cor(p$b_DiagnosisADHD, p$b_ndt_DiagnosisADHD), 3), "\n")

# `sc` is the per-participant switch cost computed in section 3; make_fig_adhd.R
# needs all three objects.
saveRDS(d, file.path(DIR, "d.rds"))
saveRDS(sc, file.path(DIR, "sc.rds"))
saveRDS(m_ln, file.path(DIR, "m_ln.rds"))
