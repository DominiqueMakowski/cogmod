# Shifted Inverse Weibull Model (Fréchet Distribution)

Provides the necessary functions to use a shifted Inverse Weibull
distribution as a custom family in `brms`. This version is
reparameterized to use `tau` (proportion of non-decision time relative
to minimum RT) and `minrt` (minimum possible RT), where the non-decision
time `ndt = tau * minrt`. The distribution is parameterized by the shape
(`mu`), scale (`sigma`), as well as `tau`, and `minrt` (which is
typically fixed).

## Usage

``` r
rt_invweibull(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

rt_invweibull_lpdf_expose()

rt_invweibull_stanvars()

log_lik_rt_invweibull(i, prep)

posterior_predict_rt_invweibull(i, prep, ...)

posterior_epred_rt_invweibull(prep)
```

## Arguments

- link_mu, link_sigma, link_tau, link_minrt:

  Link function for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.
