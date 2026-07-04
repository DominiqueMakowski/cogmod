# Linear Ballistic Accumulator (LBA) Model Simulation

Simulates random draws (reaction times and choices) from a two-choice
Linear Ballistic Accumulator (LBA) model.

In this parametrization, each accumulator has its own independent drift
rate distribution:

- Accumulator 0 has drift drawn from N(driftzero, sigmazero^2).

- Accumulator 1 has drift drawn from N(driftone, sigmaone^2).

For each trial, drift rates are sampled on an individual basis until at
least one of the two is positive. The starting point for each
accumulator is sampled uniformly from (0, sigmabias). The decision
threshold is defined as sigmabias + bs. The decision time for an
accumulator is calculated as (b - start)/drift, and if its drift is not
positive, its decision time is set to Inf. The winning accumulator (the
one whose decision time is minimal) determines the response, and the
final reaction time is the sum of its decision time and a fixed
non-decision time (ndt).

## Usage

``` r
rlba(
  n,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  bs = 0.5,
  ndt = 0.3,
  max_iter = 100
)

dlba(
  x,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  bs = 0.5,
  ndt = 0.3,
  response,
  log = FALSE
)

lba_lpdf_expose()

lba_stanvars()

lba(
  link_mu = "identity",
  link_driftone = "identity",
  link_sigmazero = "softplus",
  link_sigmaone = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

log_lik_lba(i, prep)

posterior_predict_lba(i, prep, ...)

posterior_epred_lba(prep)
```

## Arguments

- n:

  Number of simulated trials. Must be a positive integer.

- driftzero:

  Mean drift rate for the first accumulator (accumulator 0). Range:
  (-Inf, Inf).

- driftone:

  Mean drift rate for the second accumulator (accumulator 1). Range:
  (-Inf, Inf).

- sigmazero:

  Standard deviation of the drift rate for the first accumulator. Must
  be positive. Range: (0, Inf).

- sigmaone:

  Standard deviation of the drift rate for the second accumulator. Must
  be positive. Range: (0, Inf).

- sigmabias:

  Maximum starting point for the uniform distribution of starting
  evidence (0, sigmabias). Must be positive. Range: (0, Inf). Default:
  0.5.

- bs:

  Additional amount beyond `sigmabias` to set the decision threshold (b
  = sigmabias + bs). Must be positive. Range: (0, Inf). Default: 0.5.

- ndt:

  Non-decision time, representing processes such as encoding and motor
  response. Must be non-negative. Range: \[0, Inf). Default: 0.3.

- max_iter:

  Maximum iterations allowed (per trial) for resampling drift rates if
  both are non-positive. Default: 100.

- x:

  The observed reaction time (RT). Must be greater than `ndt`.

- response:

  The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.

- log:

  Logical; if TRUE, returns the log-density. Default: FALSE.

- link_mu, link_driftone, link_sigmazero, link_sigmaone, link_sigmabias,
  link_bs, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Details

**Psychological Interpretation:**

- **Drift Rate (`driftzero`, `driftone`)**: Reflects the rate at which
  evidence accumulates for each choice. Higher drift rates indicate
  faster evidence accumulation and a higher likelihood of selecting the
  corresponding choice. Differences in drift rates between accumulators
  can represent differences in preference, difficulty, or bias between
  the two options.

- **Drift Rate Variability (`sigmazero`, `sigmaone`)**: Captures
  trial-to-trial variability in the evidence accumulation process.
  Higher variability indicates less consistent evidence accumulation,
  leading to greater variability in reaction times and choices.

- **Start Point Variability (`sigmabias`)**: Represents the range of
  initial evidence levels for each accumulator. Larger values of
  `sigmabias` introduce more variability in reaction times, as the
  starting point can vary more widely between trials.

- **Threshold (`b = sigmabias + bs`)**: Boundary separation (`bs`).
  Represents the amount of evidence required to make a decision. Higher
  thresholds lead to longer reaction times but more accurate decisions,
  as more evidence is required before a choice is made.

- **Non-Decision Time (`ndt`)**: Accounts for processes unrelated to
  evidence accumulation, such as sensory encoding and motor response.
  This parameter shifts all reaction times by a constant amount.

## References

- Brown, S. D., & Heathcote, A. (2008). The simplest complete model of
  choice response time: Linear ballistic accumulation. *Cognitive
  Psychology*, *57*(3), 153-178.
  [doi:10.1016/j.cogpsych.2007.12.002](https://doi.org/10.1016/j.cogpsych.2007.12.002)

## Examples

``` r
df <- rlba(n = 1000, driftzero = 3, driftone = 2,
           sigmazero = 0.5, sigmaone = 0.5,
           sigmabias = 0.5, bs = 0.5, ndt = 0.3)
hist(df$rt[df$response == 0], breaks = 50, col = rgb(0,0,1,0.5))
hist(df$rt[df$response == 1], breaks = 50, col = rgb(1,0,0,0.5), add = TRUE)


# You can expose the lpdf function as follows:
# lba_lpdf <- lba_lpdf_expose()
# lba_lpdf(...)
```
