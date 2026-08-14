# Script used to generate data/badlm.rda
#
# A simulated dataset in which two conditions have (by construction) the same
# mean reaction time, but radically different distributions. It exists to
# demonstrate what a linear model - even a correctly specified linear *mixed*
# model - cannot see. See ?badlm and the "RT models" vignette.
#
# Condition A: short non-decision time (50 ms), wide and heavily skewed.
# Condition B: long non-decision time (450 ms), narrow and nearly symmetric.
#
# Participants differ in overall speed through an *additive* offset applied
# identically to both conditions. This matters: a multiplicative offset would
# make the two conditions' means diverge participant by participant and so
# introduce a spurious condition effect. Because the offset is additive and
# condition-independent, the within-participant contrast is exactly zero for
# everybody, and a random intercept is the correctly specified model for it.

set.seed(5)

n_participants <- 20
n_trials <- 25 # Per condition

participants <- data.frame(
  Participant = sprintf("S%02d", 1:n_participants),
  Offset = rnorm(n_participants, mean = 0, sd = 0.03)
)

badlm <- data.frame(
  Participant = rep(participants$Participant, each = 2 * n_trials),
  Condition = rep(rep(c("A", "B"), each = n_trials), n_participants)
)
badlm$Offset <- participants$Offset[match(badlm$Participant, participants$Participant)]

n <- sum(badlm$Condition == "A")
badlm$RT <- NA_real_
# Wide distribution, very short time-to-first-response (shift = 0.05)
badlm$RT[badlm$Condition == "A"] <- 0.05 + badlm$Offset[badlm$Condition == "A"] +
  rlnorm(n, meanlog = log(0.65) - 0.5^2 / 2, sdlog = 0.5)
# Narrow distribution, long non-decision time (shift = 0.45)
badlm$RT[badlm$Condition == "B"] <- 0.45 + badlm$Offset[badlm$Condition == "B"] +
  rlnorm(n, meanlog = log(0.25) - 0.15^2 / 2, sdlog = 0.15)

badlm$Participant <- factor(badlm$Participant)
badlm$Condition <- factor(badlm$Condition)
badlm$Offset <- NULL # The participant offsets are latent: drop them
badlm <- badlm[c("Participant", "Condition", "RT")]
row.names(badlm) <- NULL

# Sanity checks: the conditions must remain mean-matched and the shapes distinct
stopifnot(
  nrow(badlm) == 1000,
  abs(diff(tapply(badlm$RT, badlm$Condition, mean))) < 0.005,
  min(badlm$RT) > 0
)

save(badlm, file = "data/badlm.rda", compress = "xz")
