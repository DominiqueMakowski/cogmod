# Starting values that keep the sampler out of the flat regions

Builds an `init` argument for
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), for
the custom families whose default starting point is a bad one.

**This is not a tuning knob to reach for when sampling looks bad.** For
[`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
and
[`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
the usual default is actively harmful, and the failure is silent: the
chain does not error, it simply never moves.

## Usage

``` r
cogmod_inits(formula, data, jitter = 0.25, ...)
```

## Arguments

- formula:

  The model formula, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Must carry the family, i.e. be built with
  `brms::bf(..., family = cogmod_gamma())`.

- data:

  The data, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- jitter:

  SD of the noise added on the unconstrained scale, so that chains start
  at different points. Set to `0` for identical starts.

- ...:

  Passed to
  [`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
  and
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html),
  for arguments such as `data2`.

## Value

A function of one argument, suitable for `brms::brm(init = )`. Each call
returns a named list of starting values, one for every parameter the
generated Stan program declares.

## Details

Two regions of the parameter space have no gradient to escape on, and
both are easy to start inside.

`ndt` **too large.** `brms` initialises on the *unconstrained* scale, so
`init = 0` puts `ndt` at `exp(0) = 1` second. For sub-second reaction
times that is above nearly every observation, so every response is
attributed to the outlier component, the decision parameters drop out of
the density entirely, and their gradient is exactly zero. The chain is
stuck where it started.

**Shape below 1.** For
[`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
and
[`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
the density is unbounded at `ndt` whenever the shape is below 1, and
`init = 0` on a `softplus` link starts the shape at `log(2) = 0.69` -
inside that region.

A single scalar `init` cannot avoid both, because they pull in opposite
directions: `ndt = exp(c)` wants `c` around `-1.6`, while
`shape = softplus(c)` wants `c` above `1.9`. That is why `init = 0`
fails and why no other single number fixes it - the parameters have to
be set separately, which is what this function does.

## How it works

The Stan parameter declarations are read off
[`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
for the model you are actually fitting, and their dimensions off
[`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html),
rather than reconstructed from the formula. That is what makes it robust
to `0 + Intercept`, interactions, group-level terms and smooths:
whatever `brms` decided to call things, and however large it decided to
make them, that is what is matched.

**Every** declared parameter is given a value, not only the ones this
function has an opinion about. That is deliberate: CmdStan prints
`Init values were only set for a subset of parameters` and lists the
rest whenever the list is partial, which is noise for a list that is
partial on purpose. The parameters with no family-specific target -
regression slopes, standardized group-level effects, group-level SDs,
spline coefficients - get generic values that are at least as good as
Stan's own `U(-2, 2)`: slopes and standardized effects start at zero,
positive parameters just above their lower bound, bounded ones at the
midpoint, and Cholesky factors at the identity.

The family-specific values go to the intercept of each distributional
parameter, and to any dpar left out of the formula (which `brms`
declares as a plain auxiliary parameter). Values are set on whichever
scale the parameter lives on: the link scale for an intercept - using
the links on the family, so `cogmod_gamma(link_mu = "log")` is
honoured - and the natural scale for an auxiliary parameter. Under
`0 + Intercept` there is no separate intercept parameter, so the
coefficient named `Intercept` inside the `b` vector gets it instead.

Each chain gets a different draw. The noise is added on the
*unconstrained* scale - additive for a free parameter, multiplicative
for a positive one, on the logit scale for a doubly bounded one - so a
jittered value can never land outside its own bounds, and the chains
still start dispersed enough for `Rhat` to mean something.

`ndt` starts deliberately **small** (0.1 s, a third of its prior
median). The two errors are not symmetric: too small merely means the
shift has to grow, which the gradient will do, whereas too large removes
the gradient altogether.

## Supported families

The family is read off `formula`, so build it with
`brms::bf(..., family = cogmod_gamma())`. Every family built on the
direct `ndt` + `poutlier` parameterization is covered -
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
[`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md),
[`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md),
[`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
[`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md),
[`cogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md),
[`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
[`cogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
[`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
[`cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
[`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md),
[`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
and, for the choice-and-RT models,
[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
and
[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) -
plus
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md),
whose three parameters are all on the RT scale behind a `softplus` link
and so are equally badly served by starting at `log(2)`.

## See also

[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md),
[`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)

## Examples

``` r
d <- data.frame(RT = rcogmod_gamma(50, ndt = 0.3, poutlier = 0.02))
f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_gamma())
inits <- cogmod_inits(f, d)
inits(1)
#> $Intercept
#> [1] 1.477312
#> 
#> $Intercept_sigma
#> [1] -1.124586
#> 
#> $Intercept_ndt
#> [1] -2.195298
#> 
#> $Intercept_poutlier
#> [1] -3.861294
#> 

# \donttest{
# Fitting needs cmdstanr, which lives outside CRAN - see the package website.
if (requireNamespace("cmdstanr", quietly = TRUE) &&
    !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
  m <- brms::brm(f,
    data = d, prior = cogmod_priors(f, d),
    stanvars = cogmod_stanvars(f), init = cogmod_inits(f, d),
    backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
  )
}
# }
```
