# Shifted Log-Weibull Model (Gumbel Distribution)

Provides the necessary functions to use a shifted Log-Weibull
distribution as a custom family in `brms`. This version is
reparameterized to use `tau` (proportion of non-decision time relative
to minimum RT) and `minrt` (minimum possible RT), where the non-decision
time `ndt = tau * minrt`. The distribution is parameterized by the shape
(`mu`), scale (`sigma`), as well as `tau`, and `minrt` (which is
typically fixed).

## Usage

``` r
rt_logweibull(
  link_mu = "identity",
  link_sigma = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

rt_logweibull_lpdf_expose()

rt_logweibull_stanvars()

log_lik_rt_logweibull(i, prep)

posterior_predict_rt_logweibull(i, prep, ...)

posterior_epred_rt_logweibull(prep)
```

## Arguments

- link_mu, link_sigma, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.
