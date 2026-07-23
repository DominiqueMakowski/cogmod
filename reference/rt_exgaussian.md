# Ex-Gaussian Model (Classical Parameterization) for brms

Provides the necessary functions to use the "classical" parameterization
of the Ex-Gaussian distribution as a custom family in `brms`. Unlike
`brms`'s built-in
[`exgaussian()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family - in which `mu` indexes the mean of the *entire* distribution
(Gaussian + exponential components combined) - this parameterization
follows the convention familiar to experimental psychologists, in which
`mu` and `sigma` are the mean and SD of the Gaussian component alone,
and `tau` is the mean of the exponential component (the tail). The mean
of the full distribution is `mu + tau`.

## Usage

``` r
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

This distinction matters because changes in the Gaussian location (`mu`)
and changes in the exponential tail (`tau`) can offset one another at
the level of the overall mean, so effects estimated on `brms`'s default
`mu` can lead to different (and potentially incorrect) inferences than
effects estimated on this classical `mu`.

By default, all three parameters use a `"softplus"` link
(`log(1 + exp(x))`), rather than `"identity"` or `"log"`. `mu` must
remain strictly positive since it represents the (unobserved) center of
the Gaussian component. An `"identity"` link would allow the linear
predictor to cross zero or go negative, which is invalid for all three
parameters. A `"log"` link would enforce positivity too, but its
curvature explodes as the linear predictor departs from zero, producing
extreme gradients and making priors/sampling harder to calibrate,
especially for `mu` and `tau`, which are already both on the RT scale
(seconds) and can take comparatively large values. `"softplus"` is
positive-constrained like `"log"` but behaves almost linearly
(`softplus(x) ~ x`) away from zero, making it easier to specify
weakly-informative priors directly on the RT scale while still
guaranteeing valid, strictly positive parameters.
