# Shifted Gamma Model for brms

Provides the necessary functions to use a shifted Gamma distribution as
a custom family in `brms`. This version is parameterized using `tau`
(proportion of non-decision time relative to minimum RT) and `minrt`
(minimum possible RT), where the non-decision time `ndt = tau * minrt`.
The distribution is parameterized by the shape (`mu` in brms, alpha in
standard Gamma) and scale (`sigma` in brms, beta in standard Gamma), as
well as `tau`, and `minrt` (which is typically fixed).

## Usage

``` r
rt_gamma(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

rt_gamma_lpdf_expose()

rt_gamma_stanvars()

log_lik_rt_gamma(i, prep)

posterior_predict_rt_gamma(i, prep, ...)

posterior_epred_rt_gamma(prep)
```

## Arguments

- link_mu, link_sigma, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Value

A
[`brms::custom_family`](https://paulbuerkner.com/brms/reference/custom_family.html)
object.
