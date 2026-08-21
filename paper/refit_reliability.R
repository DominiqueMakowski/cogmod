# Refits the mixed shifted-LogNormal used for @tbl-dvour in Example 2, and
# prints the D-vour table. Mirrors the `Assessing Reliability` vignette exactly,
# so the paper and the vignette report the same fit.
#
#   Rscript refit_reliability.R      (run from the `paper/` directory)
#
# Writes ../vignettes/models/m_reliability.rds.

library(dplyr)
library(brms)
library(cogmod)

set.seed(14)

n_participants <- 30
n_trials <- 80 # Per condition
ndt <- 0.15 # Non-decision time (s): no response can occur before that

participants <- data.frame(
  Participant = sprintf("S%02d", 1:n_participants),
  Intercept = rnorm(n_participants, mean = log(0.6), sd = 0.20),
  Effect_B = rnorm(n_participants, mean = 0.10, sd = 0.01),
  Effect_C = rnorm(n_participants, mean = 0.00, sd = 0.30)
)

sim <- expand.grid(
  Trial = 1:n_trials,
  Condition = c("A", "B", "C"),
  Participant = participants$Participant,
  stringsAsFactors = FALSE
) |>
  left_join(participants, by = "Participant") |>
  mutate(
    meanlog = Intercept +
      if_else(Condition == "B", Effect_B, 0) +
      if_else(Condition == "C", Effect_C, 0),
    RT = ndt + rlnorm(n(), meanlog = meanlog, sdlog = 0.22),
    Participant = factor(Participant),
    Condition = factor(Condition)
  ) |>
  arrange(Participant, Condition, Trial) |>
  select(Participant, Trial, Condition, RT)

f <- bf(
  RT ~ Condition + (Condition | Participant),
  sigma ~ Condition + (Condition | Participant),
  ndt ~ Condition + (Condition | Participant),
  family = cogmod_lognormal()
)

m <- brm(f,
  data = sim,
  prior = cogmod_priors(f, sim),
  init = cogmod_inits(f, sim),
  stanvars = cogmod_stanvars(f),
  chains = 4, cores = 4, iter = 1000, seed = 14, backend = "cmdstanr"
)

saveRDS(m, file = "../vignettes/models/m_reliability.rds")
print(summary(m))

# ---- The numbers quoted in the paper ---------------------------------------

cat("\n==== fixed effects on mu (response formula) ====\n")
print(bayestestR::describe_posterior(m, effects = "fixed", parameters = "^b_Condition"))

cat("\n==== random-effect SDs ====\n")
print(as.data.frame(brms::VarCorr(m)$Participant$sd))

cat("\n==== D-vour ====\n")
random <- modelbased::estimate_grouplevel(m)
dv <- performance::performance_dvour(random)
print(as.data.frame(dv), digits = 3)

# Empirical difference scores, for the comparison in the text
emp <- sim |>
  group_by(Participant, Condition) |>
  summarise(m = mean(RT), .groups = "drop") |>
  tidyr::pivot_wider(names_from = Condition, values_from = m) |>
  mutate(B_A = B - A, C_A = C - A)
cat("\nempirical between-participant SD of B - A:", round(sd(emp$B_A), 4), "s\n")
cat("empirical between-participant SD of C - A:", round(sd(emp$C_A), 4), "s\n")

# Bootstrap the per-participant standard error of the B - A difference
set.seed(1)
ses <- sapply(levels(sim$Participant), function(p) {
  d <- sim[sim$Participant == p, ]
  a <- d$RT[d$Condition == "A"]
  b <- d$RT[d$Condition == "B"]
  sd(replicate(2000, mean(sample(b, replace = TRUE)) - mean(sample(a, replace = TRUE))))
})
cat("mean bootstrap SE of the per-participant B - A difference:",
    round(mean(ses), 4), "s\n")

# D-vour on the raw empirical difference scores
emp_dvour <- var(emp$B_A) / (var(emp$B_A) + mean(ses^2))
cat("D-vour of the empirical B - A difference scores:", round(emp_dvour, 3), "\n")
