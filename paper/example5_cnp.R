# =============================================================================
# Example 5: A Clinical Illustration
#
# Task-switching data from the UCLA Consortium for Neuropsychiatric Phenomics
# (Poldrack et al., 2016), OpenNeuro ds000030. The data are NOT redistributed
# with this repository; this script downloads the trial-level event files
# directly from the OpenNeuro S3 bucket.
#
# Runtime: download ~1 min; the brms fit took ~27 min on 16 cores.
# =============================================================================

library(dplyr)
library(brms)
library(cogmod)
library(bayestestR)

BASE <- "https://s3.amazonaws.com/openneuro.org/ds000030"
DIR  <- "cnp_data"
dir.create(file.path(DIR, "events"), recursive = TRUE, showWarnings = FALSE)

# ---- 1. Download -------------------------------------------------------------

pfile <- file.path(DIR, "participants.tsv")
if (!file.exists(pfile)) download.file(file.path(BASE, "participants.tsv"), pfile)
participants <- read.delim(pfile)

subs <- participants |>
  filter(diagnosis %in% c("CONTROL", "ADHD"), taskswitch == "1") |>
  pull(participant_id)

for (s in subs) {
  out <- file.path(DIR, "events", paste0(s, ".tsv"))
  if (file.exists(out) && file.size(out) > 500) next
  url <- sprintf("%s/%s/func/%s_task-taskswitch_events.tsv", BASE, s, s)
  try(download.file(url, out, quiet = TRUE), silent = TRUE)
}

# ---- 2. Assemble trial-level data --------------------------------------------

d <- lapply(subs, function(s) {
  f <- file.path(DIR, "events", paste0(s, ".tsv"))
  if (!file.exists(f)) return(NULL)
  e <- read.delim(f)
  data.frame(
    Participant = s,
    Diagnosis   = participants$diagnosis[participants$participant_id == s],
    RT          = e$ReactionTime,
    Correct     = e$CorrectResp,
    Switching   = e$Switching,
    Congruency  = e$Congruency
  )
}) |> bind_rows()

# Correct responses only; standard RT trimming. RT == 0 codes a non-response.
d <- d |>
  filter(Correct == 1, RT > 0.2, RT < 2) |>
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
  tau ~ Diagnosis + (1 | Participant),
  minrt = min(d$RT),
  family = rt_lognormal()
)

m_ln <- brm(
  f,
  data = d,
  stanvars = rt_lognormal_stanvars(),
  chains = 4, cores = 4, iter = 1000, warmup = 500,
  threads = threading(3),
  backend = "cmdstanr", seed = 42
)

print(summary(m_ln))

# ---- 5. Response-scale decomposition -----------------------------------------

p       <- posterior::as_draws_df(m_ln)
softplus <- function(x) log1p(exp(x))
minrt   <- min(d$RT)

ndt_c <- plogis(p$b_tau_Intercept) * minrt
ndt_a <- plogis(p$b_tau_Intercept + p$b_tau_DiagnosisADHD) * minrt
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

# The decision and non-decision parameters trade off; report the correlation
cat("posterior cor(mu_ADHD, tau_ADHD) =",
    round(cor(p$b_DiagnosisADHD, p$b_tau_DiagnosisADHD), 3), "\n")

saveRDS(d, file.path(DIR, "d.rds"))
saveRDS(m_ln, file.path(DIR, "m_ln.rds"))
