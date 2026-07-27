# Ex-Gaussian Model (Classical Parameterization)

Density, random generation, and `brms` custom family for the Ex-Gaussian
distribution, using the "classical" parameterization familiar to
experimental psychologists, in which `mu` and `sigma` are the mean and
SD of the Gaussian component alone, and `tau` is the mean of the
exponential component (the tail). This is unlike `brms`'s built-in
[`exgaussian()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family, in which `mu` indexes the mean of the *entire* distribution
(Gaussian + exponential components combined).

Functions:

- `rrt_exgaussian()`: Simulates random draws from the Ex-Gaussian
  distribution.

- `drt_exgaussian()`: Computes the density (likelihood) of the
  Ex-Gaussian distribution.

- `rt_exgaussian()`: Creates a
  [`brms::custom_family()`](https://paulbuerkner.com/brms/reference/custom_family.html)
  for use in `brms` models.

## Usage

``` r
rrt_exgaussian(n, mu = 0.5, sigma = 0.1, tau = 0.2)

drt_exgaussian(x, mu = 0.5, sigma = 0.1, tau = 0.2, log = FALSE)

rt_exgaussian(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_tau = "softplus"
)

rt_exgaussian_lpdf_expose()

rt_exgaussian_stanvars()

log_lik_rt_exgaussian(i, prep)

posterior_predict_rt_exgaussian(i, prep, ...)

posterior_epred_rt_exgaussian(prep)
```

## Arguments

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- mu:

  Mean of the Gaussian component. Must be positive. Range: (0, Inf).

- sigma:

  SD of the Gaussian component. Must be positive. Range: (0, Inf).

- tau:

  Mean of the exponential component (the tail). Must be positive. Range:
  (0, Inf).

- x:

  Vector of quantiles (observed reaction times).

- log:

  Logical; if TRUE, probabilities p are given as log(p).

- link_mu, link_sigma, link_tau:

  Character of the type of link used to model the ex-Gaussian
  parameters. Defaults to `"softplus"` for all three (see Details).

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Value

A
[`brms::custom_family`](https://paulbuerkner.com/brms/reference/custom_family.html)
object.

## Details

The Ex-Gaussian distribution is the sum of an independent Normal
(Gaussian) random variable with mean `mu` and SD `sigma`, and an
Exponential random variable with mean `tau` (rate `1 / tau`). Unlike
`brms`'s built-in
[`exgaussian()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family - in which `mu` indexes the mean of the *entire* distribution
(Gaussian + exponential components combined) - here `mu` is the mean of
the Gaussian component alone, so that the mean of the full distribution
is `mu + tau`.

This distinction matters because changes in the Gaussian location (`mu`)
and changes in the exponential tail (`tau`) can offset one another at
the level of the overall mean, so effects estimated on `brms`'s default
`mu` can lead to different (and potentially incorrect) inferences than
effects estimated on this classical `mu`.

In the `brms` custom family (`rt_exgaussian()`), all three parameters
use a `"softplus"` link (`log(1 + exp(x))`) by default, rather than
`"identity"` or `"log"`. `mu` must remain strictly positive since it
represents the (unobserved) center of the Gaussian component. An
`"identity"` link would allow the linear predictor to cross zero or go
negative, which is invalid for all three parameters. A `"log"` link
would enforce positivity too, but its curvature explodes as the linear
predictor departs from zero, producing extreme gradients and making
priors/sampling harder to calibrate, especially for `mu` and `tau`,
which are already both on the RT scale (seconds) and can take
comparatively large values. `"softplus"` is positive-constrained like
`"log"` but behaves almost linearly (`softplus(x) ~ x`) away from zero,
making it easier to specify weakly-informative priors directly on the RT
scale while still guaranteeing valid, strictly positive parameters.

## References

- Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation
  of the ex-Gaussian and shifted Wald parameters: A diffusion model
  analysis. *Psychonomic Bulletin & Review*, *16*(5), 798–817.
  [doi:10.3758/PBR.16.5.798](https://doi.org/10.3758/PBR.16.5.798)

## Examples

``` r
# Simulate 1000 RTs
rts <- rrt_exgaussian(1000, mu = 0.5, sigma = 0.1, tau = 0.2)
# hist(rts, breaks = 50, main = "Simulated Ex-Gaussian RTs", xlab = "Reaction Time")
```
