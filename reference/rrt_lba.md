# Single-Accumulator LBA model

This function simulates reaction times using a single-accumulator
version of the Linear Ballistic Accumulator (LBA) model.

## Usage

``` r
rrt_lba(n, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3)

drt_lba(
  x,
  drift = 3,
  sigma = 1,
  sigmabias = 0.5,
  bs = 0.5,
  ndt = 0.3,
  log = FALSE
)

rt_lba_lpdf_expose()

rt_lba_stanvars()

rt_lba(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

log_lik_rt_lba(i, prep)

posterior_predict_rt_lba(i, prep, ...)

posterior_epred_rt_lba(prep)
```

## Arguments

- n:

  A single positive integer giving the number of trials.

- drift:

  Mean drift rate.

- sigma:

  Standard deviation of the drift rate.

- sigmabias:

  The starting-point range (A); must be positive.

- bs:

  The additional offset such that b = sigmabias + bs; must be positive.

- ndt:

  Non-decision time.

- x:

  The observed reaction time (RT). Must be greater than `ndt`.

- log:

  Logical; if TRUE, returns the log-density. Default: FALSE.

- link_mu, link_sigma, link_sigmabias, link_bs, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Examples

``` r
# Simulate 1000 trials with specified parameters
rts <- rrt_lba(n = 1000, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3)
hist(rts, breaks = 100, xlab = "RT (s)", col = "lightblue")

```
