# Assessing Reliability and Interindividual Variability

``` r

library(cogmod)
library(easystats)
library(dplyr)
library(ggplot2)
library(brms)
library(cmdstanr)

options(mc.cores = parallel::detectCores() - 2)
```

Imagine we have designed a brand new cognitive task. It has three
conditions - **A**, **B** and **C** - and our ambition is to use it in a
**clinical context**: to quantify something about *individuals*, so that
a patient’s score can be related to their symptoms, compared to a norm,
or tracked over time.

The first thing we do is run a validation study and test whether the
conditions affect reaction times. As we will see, only one of them turns
out to have a significant effect - and the tempting conclusion is that
this effect is the one the task actually captures, and therefore the one
we should be computing for each patient.

That conclusion would be premature, and most likely wrong. An average
effect tells us what the manipulation does *to the group*. It says
nothing about whether it does the same thing to everybody, and nothing
about whether we can measure it in a single person with any precision.
Yet those two properties - **interindividual variability** and
**reliability** - are exactly what determines whether a task can serve
as an individual-level instrument, and they have to be assessed on their
own terms.

## Data Description

Our validation study is a typical repeated-measures design: 30
participants each perform 80 trials in each of the three conditions, and
we record their reaction times.

Code

``` r

# Because we are simulating the data, we get to decide what is true, and we make
# the two effects *deliberately* different:
# 
# - **A vs. B**: a small but extremely **consistent** effect. Every participant
#   is slowed down by about the same amount (the between-participant SD of the
#   effect is tiny, 0.01). There is essentially no interindividual variability.
# - **A vs. C**: an effect with a population mean of **exactly zero**, but a very
#   large between-participant SD (0.30). Some participants are dramatically
#   slowed down by condition C, others are dramatically speeded up.

set.seed(14)

n_participants <- 30
n_trials <- 80  # Per condition

participants <- data.frame(
  Participant = sprintf("S%02d", 1:n_participants),
  # Overall speed: participants differ substantially from one another
  Intercept = rnorm(n_participants, mean = log(0.6), sd = 0.20),
  # A vs. B: small effect, (almost) identical for everybody
  Effect_B = rnorm(n_participants, mean = 0.10, sd = 0.01),
  # A vs. C: no effect on average, but huge interindividual variability
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
    # Trial-level mean on the log scale
    meanlog = Intercept +
      if_else(Condition == "B", Effect_B, 0) +
      if_else(Condition == "C", Effect_C, 0),
    # Generate RTs, adding within-participant (trial-to-trial) noise
    RT = rlnorm(n(), meanlog = meanlog, sdlog = 0.22),
    Participant = factor(Participant),
    Condition = factor(Condition)
  ) |>
  arrange(Participant, Condition, Trial) |>
  select(Participant, Trial, Condition, RT)
```

Looking at the raw distributions, the three conditions are nearly
indistinguishable, and nothing suggests that B and C differ in any
interesting way.

``` r

ggplot(sim, aes(x = RT, fill = Condition)) +
  geom_density(alpha = 0.5) +
  scale_fill_manual(values = c("A" = "grey40", "B" = "purple", "C" = "orange")) +
  labs(x = "Reaction Time (s)", y = "Density") +
  theme_minimal()
```

![](reliability_files/figure-html/unnamed-chunk-2-1.png)

## Empirical Estimates

Let us now compute the score we would hand to a clinician: for each
participant, average the RTs within each condition and subtract, giving
an “empirical” effect of B - A and of C - A (expressed here in
milliseconds).

``` r

empirical <- sim |>
  summarise(RT = mean(RT), .by = c(Participant, Condition)) |>
  reframe(
    Effect = c("B - A", "C - A"),
    Difference = c(
      RT[Condition == "B"] - RT[Condition == "A"],
      RT[Condition == "C"] - RT[Condition == "A"]
    ) * 1000,
    .by = Participant
  )

empirical |>
  summarise(
    Mean = mean(Difference),
    SD = sd(Difference),
    Min = min(Difference),
    Max = max(Difference),
    .by = Effect
  )
#>   Effect      Mean        SD       Min      Max
#> 1  B - A 73.956648  30.80423   22.2618 159.6227
#> 2  C - A  0.444845 170.68889 -240.4986 408.4365
```

The two effects could not be more different, and yet **the number we
usually report - the mean - is blind to it**. The B - A effect is about
+74 ms on average, with a between-participant SD of ~31 ms. The C - A
effect has a mean of essentially zero (+0.4 ms) - but a
between-participant SD of ~171 ms, and it ranges from -240 ms to +408 ms
depending on which participant we look at.

``` r

# Sort participants by the size of their C - A effect
order_C <- empirical |>
  filter(Effect == "C - A") |>
  arrange(Difference) |>
  pull(Participant)

empirical |>
  mutate(Participant = factor(Participant, levels = order_C)) |>
  ggplot(aes(x = Difference, y = Participant, color = Effect)) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_point(size = 2.5) +
  facet_wrap(~Effect) +
  scale_color_manual(values = c("B - A" = "purple", "C - A" = "orange")) +
  labs(x = "RT difference (ms)", y = "Participant") +
  theme_minimal() +
  theme(legend.position = "none")
```

![](reliability_files/figure-html/unnamed-chunk-4-1.png)

The picture is unambiguous (note that participants are ordered by the
size of their C - A effect, so that a given row refers to the same
person in both panels). On the left, every single participant shows the
same small slowing in condition B: the effect is *reliable* in the sense
of being reproducible across people, but there is almost nothing to
distinguish one participant from another. On the right, participants are
spread all over the place, with some being massively slower and others
massively faster in condition C - and these differences cancel out to
nothing at the group level.

This is bad news for our original plan. Condition B - the one that
reached significance, and on which we were about to build our clinical
index - gives essentially the same score to everybody, and a measure
that does not vary between people cannot possibly correlate with
symptoms, separate patients from controls, or track change. Condition
C - the “null” effect that we would have discarded - is the one carrying
the individual differences we are after.

> **Careful: empirical estimates are noisy!**

These per-participant differences are *not* clean measurements of each
participant’s true effect. Each one is computed from a finite number of
trials and therefore carries measurement error. Note for instance that
the B - A effects show an apparent spread of ~31 ms even though we
generated them with a between-participant SD of only ~6 ms (0.01 on the
log scale): almost all of that visible spread is noise, not real
interindividual variability. Disentangling the two is precisely what the
model below does.
