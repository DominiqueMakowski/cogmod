# Priors that make a cogmod posterior proper

Fills in weakly informative priors for the parameters `brms` would
otherwise leave flat, for the model you are actually fitting.

**For
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
this is not a convenience.** `brms` assigns a flat, improper prior to
the intercept of any custom-family parameter it does not recognise,
which there means both `ndt` and `poutlier`. The likelihood has two
directions in which it is exactly flat: `poutlier` toward 1, where every
response is attributed to the outlier component and `mu`, `sigma` and
`ndt` drop out of the density altogether; and `ndt` toward 0, where the
model reduces to an unshifted LogNormal and the gradient with respect to
`log(ndt)` vanishes. Flat prior plus infinite flat region is an improper
posterior. The fit does not fail loudly - it returns intercepts around
`1e14` with `Rhat` near 2 and an effective sample size of about 5.

The second direction has nothing to do with the mixture; it is inherent
to putting a positive shift on a log link, which is why a prior on
`poutlier` alone is not enough.

## Usage

``` r
cogmod_priors(formula, data, ...)
```

## Arguments

- formula:

  The model formula, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Must carry the family, i.e. be built with
  `brms::bf(..., family = cogmod_lognormal())`.

- data:

  The data, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).

- ...:

  Passed to
  [`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
  and
  [`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html),
  for arguments such as `data2` or `knots`.

## Value

A `brmsprior` object, to pass to `brms::brm(prior = )`.

## Details

The function starts from
[`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
for the model in hand, edits the rows it knows how to improve, and
returns the result of
[`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html).
Three things follow.

Every row comes from the model `brms` is going to build, so a prior
matching no parameter is impossible by construction: `0 + Intercept`
formulas, interactions, group-level terms and smooths are all handled
because none of them are guessed at.

The return value is a `brmsprior` object: the same class returned by
[`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
and
[`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html).
It prints as a table, including the defaults that `brms` will use, but
can be combined with [`c()`](https://rdrr.io/r/base/c.html) like any
other `brms` prior specification.

Passing it through
[`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html)
means a malformed specification errors here, with the offending row in
view, rather than deep inside
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Setting your own priors

Combine it with
[`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
entries. Set `replace = TRUE` to replace one of these defaults, or omit
it to add a prior for a different slot:

    priors <- c(
      cogmod_priors(f, df),
      brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
      replace = TRUE
    )

## Supported families

The family is read off `formula`, so build it with
`brms::bf(..., family = cogmod_lognormal())`. Every family built on the
direct `ndt` + `poutlier` parameterization is edited:
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
[`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md)
and, for the choice-and-RT models,
[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
and
[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md).
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
is edited too, although it is not built on that parameterization - see
its own section below. Any other family, or a formula carrying none, is
passed through: you get a message and `brms`'s own defaults, unchanged,
so the call is always safe to leave in a script.

What gets set, on the link scale (`log` for `ndt`, `logit` for
`poutlier`, `identity` for `shape`):

|  |  |  |  |
|----|----|----|----|
| class | `ndt` | `poutlier` | `shape` |
| `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-1.2, 0.2)` | `normal(-5, 1)` | `normal(0, 0.5)` |
| `b` (slopes) | `normal(0, 0.2)` | `normal(0, 0.2)` | `normal(0, 0.2)` |
| `sd`, `sds` | `exponential(1)` | `exponential(1)` | `exponential(1)` |

`normal(-1.2, 0.2)` puts `ndt` at roughly 170 to 300 ms. Like everything
else in these families it is stated in **seconds**, which is the unit
the package expects throughout. `poutlier` is a proportion and does not
move: `normal(-5, 1)` is centred at about 0.7% and puts roughly 95% of
its mass between 0.1% and 5%. That is where the empirical estimates
sit - pooling across four lexical-decision megastudies, Miller (2024)
puts the outlier proportion below 0.5% and argues that the 5-10% assumed
in most simulation work is unrealistically large.

`shape` exists only for
[`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md).
`normal(0, 0.5)` is centred on the LogNormal shape and keeps the sampler
clear of `sigma * shape >= 1`, where the decision density becomes
unbounded at the shift and the likelihood with it.

A family may add rows of its own where its likelihood has a flat
direction that `brms` would leave improper. Seven do:

- [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md):
  `dof`. The LogNormal is only reached as `dof -> Inf`, and it is
  approached slowly - at `dof = 100` the density still differs from the
  LogNormal by 66% in the tail - so above about `dof = 30` the
  likelihood has almost stopped moving, and a `log` link puts that
  region at plus infinity. The prior fences the other end too: a
  Student-t is symmetric on the log scale, so a small `dof` piles mass
  just above `ndt` as well as in the slow tail, which is what `poutlier`
  is for. `normal(1.8, 0.7)` centres `dof` at 6 with 95% between 1.5 and
  24.

- [`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md):
  `tau`, the mean of the exponential residual stage. It is a length of
  time in seconds behind a `softplus` link, which `brms` has no way to
  know, and it shares a ridge with `ndt` - both delay the response, and
  only the shape of the leading edge tells them apart.
  `normal(-1.5, 0.7)` is the same statement
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  gets for the same quantity: roughly 55 to 630 ms.

- [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md):
  `sigmabias` and `boundary`. As the start-point range approaches zero
  the LBA converges to the recinormal and the likelihood stops depending
  on it, and a `softplus` link reaches zero only at minus infinity.
  Without those rows `sigmabias` runs off - `softplus(-10.4)`, `Rhat`
  1.69.

- [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md):
  `sigmabias` and `boundary`, for exactly the reason
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  needs them - the two enter the threshold only through the sum
  `b = boundary + sigmabias` and trade off along a ridge worth a handful
  of log units, and the `softplus` link reaches zero only at minus
  infinity. The drift rates are left alone: a zero-drift accumulator
  still finishes, and still wins sometimes, so there is no plateau at
  the bottom of that direction the way there is for
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)'s
  `nuone`.

- [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md):
  `sigmazero`, `sigmaone`, `sigmabias` and `boundary`. The last two
  share the threshold ridge above; the two drift SDs are worse. The
  evidence scale of an LBA is arbitrary - multiply the drifts, their
  SDs, the start-point range and the threshold by any `c > 0` and every
  finishing time is unchanged - so the likelihood is *exactly* constant
  along that ray. These priors make the posterior proper; only fixing
  one SD in the formula (`sigmazero = 1`) identifies the scale. See
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md).

- [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md):
  `sigmadrift`, `sigmabias` and `sigmandt`, the three between-trial
  variability parameters. Each has a floor at zero that its link reaches
  only at minus infinity, and the likelihood stops changing well before
  then. They are also the parameters a DDM is least able to recover, so
  these priors are deliberately tighter than the rest. Fixing the ones a
  design cannot identify (`sigmadrift = 0` in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)) is
  usually better than estimating them behind a prior.

- [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md):
  `sigmadrift` and `sigmandt`, for the same reasons as the DDM's, with
  `sigmandt` taking the DDM's prior for the same quantity verbatim. Both
  should be fixed at zero unless the design can identify them;
  `sigmandt` in particular needs a lot of data, a strong prior, or both.

- [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md):
  `nuone`, `sigmazero` and `sigmaone`. Push an accumulator's rate far
  enough down and it never finishes first, so the density depends on it
  only through a survival term that has already saturated at 1; past
  about `nuone = -6` the log-likelihood is exactly constant, and that
  accumulator's `sigma` is unidentified along with it. `mu` - which is
  `nuzero` - has the mirror-image plateau, but it is the response's own
  intercept and `brms` already gives it a proper `student_t` default, so
  it is left alone. Model a rarely-chosen option and it is worth
  mirroring the `nuone` prior onto it.

All of them are the same failure as `ndt` and `poutlier`: an infinite
flat region under a flat prior. See
[`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md),
[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md),
[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
and
[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md).

Note that no prior is set on the shape of
[`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
or
[`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
although a shape below 2 makes their `ndt` gradient unbounded. That
region is reached because the *likelihood* prefers it, by around 100 log
units on the data in `vignette("rt_models")`, so a prior weak enough to
be a sensible default cannot move the posterior out of it - only bias
it.
[`?rcogmod_weibull`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
sets out what to do instead.

## The ex-Gaussian

[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
has neither `ndt` nor `poutlier`, but all three of its parameters are
lengths of time in **seconds**, and `brms` has no way to know that.
`sigma` and `tau` sit behind a `softplus` link; `mu` is on `identity`.
All three intercepts are set.

|  |  |  |
|----|----|----|
| class | `sigma` | `tau` |
| `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-2.3, 0.7)` | `normal(-1.5, 0.7)` |
| `b` (slopes) | `normal(0, 0.5)` | `normal(0, 0.5)` |
| `sd`, `sds` | `exponential(1)` | `exponential(1)` |

On the `softplus` scale `normal(-2.3, 0.7)` puts `sigma` - the SD of the
Gaussian component - between roughly 25 and 330 ms with a median of 96
ms, and `normal(-1.5, 0.7)` puts `tau` - the mean of the exponential
tail - between roughly 55 and 630 ms with a median of 201 ms. Both cover
the range these parameters occupy across the usual simple- and choice-RT
tasks and are wide enough not to fight data that disagrees.

`tau` arrives from `brms` flat, so it is filled like any other
unrecognised custom dpar. `sigma` arrives with a **non-empty** default,
`student_t(3, 0, 2.5)`, because `brms` recognises the name from its own
families - centred on `softplus(0) = 0.69 s`, a Gaussian SD wider than
most whole RT distributions. That one is overridden rather than filled,
the same treatment `shape` and an omitted `ndt` get above.

`mu` gets `normal(0.4, 0.25)` on its own intercept - 95% of the mass
between -0.09 and 0.89 s. It is not left to `brms`, whose
`student_t(3, 0, 2.5)` is a fair statement about a location on a
`softplus` link (median 0.69 s) but not on the `identity` link `mu`
uses, where it is centred on zero seconds and puts a Gaussian centre of
-2 s on a par with one of +2 s. The prior does **not** exclude negative
values: `mu` is a location, and for fast heavily-tailed data the
Gaussian component genuinely belongs near or below zero with `tau`
carrying the mass. Only the intercept is set - the response's slopes are
the effects being estimated, and are left to `brms`. Note that `mu` is
the centre of the Gaussian component **alone**, so the mean of the
distribution it implies is `mu + tau`; see
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md).

## Parameters left out of the formula

Writing `ndt ~ 1` and omitting `ndt` entirely are not the same thing to
`brms`, and the difference matters here. A dpar that appears in
[`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) -
even as `~ 1` - gets a linear predictor, so it is estimated on the
**link** scale and reported under *Regression Coefficients* as
`ndt_Intercept`. A dpar left out is declared as a plain auxiliary
parameter instead, the same mechanism as `sigma` for
[`gaussian()`](https://rdrr.io/r/stats/family.html), estimated on the
**natural** scale with no link and reported under *Further
Distributional Parameters*.

Both forms are filled, with priors on whichever scale the parameter
actually lives on:

|  |  |  |
|----|----|----|
| dpar | in [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) (link scale) | omitted (natural scale) |
| `ndt` | `normal(-1.2, 0.2)` | `lognormal(-1.2, 0.2)` |
| `poutlier` | `normal(-5, 1)` | `exponential(100)` |
| `shape` | `normal(0, 0.5)` | `normal(0, 0.5)` |
| `sigmabias`, `boundary` ([`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md), [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md), [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
| `sigmazero`, `sigmaone` ([`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
| `sigmadrift` ([`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)) | `normal(0, 1)` | `lognormal(-1, 0.75)` |
| `sigmabias` ([`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)) | `normal(-2, 1)` | `beta(1, 5)` |
| `sigmandt` ([`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md), [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)) | `normal(-3, 1)` | `lognormal(-3, 1)` |
| `sigmadrift` ([`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
| `nuone` ([`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)) | `normal(0.7, 1.5)` | `normal(0.7, 1.5)` |
| `sigmazero`, `sigmaone` ([`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
| `dof` ([`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md)) | `normal(1.8, 0.7)` | `lognormal(1.8, 0.7)` |
| `tau` ([`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md)) | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` |
| `mu` ([`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)) | `normal(0.4, 0.25)` | \- (always modelled) |
| `sigma` ([`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)) | `normal(-2.3, 0.7)` | `lognormal(-2.3, 0.7)` |
| `tau` ([`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)) | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` |

The `ndt` pair describes the same belief twice: `lognormal` is just
`normal` on the log scale, written for the untransformed parameter. The
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
pairs do the same to within a rounding error, because `softplus(x)` and
`exp(x)` agree to three figures for the `x` below `-2` that both
parameters live at: `lognormal(-2.3, 0.7)` has median 0.100 against
`softplus(-2.3) = 0.096`.

If the data were trimmed before fitting, tighten this rather than
removing it: `normal(-7, 0.5)` asserts essentially no contamination
while keeping the density positive below `ndt`. Fixing `poutlier = 0`
outright reinstates the hard min-RT boundary that the component exists
to remove. See the *Trimmed data* section of
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md).

`poutlier` is deliberately *not* the same belief twice. Leaving it out
of the formula is itself information - you either trimmed the data
already or do not expect outliers - so the omitted form puts its **mode
at zero**, which a logit-scale prior cannot do at any location. The
centre is unchanged: `exponential(100)` has median 0.0069 against
`plogis(-5) = 0.0067`. It is still a prior rather than a constraint, so
a genuine spike of fast responses will still pull the rate up; to switch
the parameter off entirely, trim the data and it will simply sit near
zero.

Two rows override a **non-empty** `brms` default rather than filling an
empty one. `brms` recognises the name `ndt` from its own shifted
families and supplies `uniform(0, min_Y)` - precisely the min-RT bound
that
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)'s
parameterization exists to remove, reimposed silently and with a warning
about an upper bound on an unbounded parameter. It recognises `sigma`
too, and gives
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)'s
a half `student_t(3, 0, 2.5)`, whose median of 1.9 s is a Gaussian SD
wider than most whole RT distributions. An omitted `poutlier` is left
flat over `[0, 1]` by `brms`, which is proper but puts half its mass
above 0.5, and an omitted `shape` or `tau` is flat over the whole real
line, which is not proper at all.

Slope and group-level priors are deliberately narrow. On a log or a
logit link a flat slope prior is not as harmless as it looks, and a
group-level SD with no prior can wander far enough for individual groups
to reach the flat regions above even when the population intercept is
well behaved.

## References

- Miller, J. (2024). Estimating the proportions and latencies of
  reaction time outliers: A pooling method and case study of lexical
  decision tasks. *Behavior Research Methods*, *56*(7), 7280-7306.
  [doi:10.3758/s13428-024-02419-y](https://doi.org/10.3758/s13428-024-02419-y)

## See also

[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md),
[`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)

## Examples

``` r
d <- data.frame(RT = rcogmod_lognormal(50, ndt = 0.3, poutlier = 0.02))
f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
  family = cogmod_lognormal()
)
cogmod_priors(f, d)
#>                   prior     class coef group resp     dpar nlpar lb ub tag
#>  student_t(3, 0.8, 2.5) Intercept                                         
#>       normal(-1.2, 0.2) Intercept                      ndt                
#>           normal(-5, 1) Intercept                 poutlier                
#>    student_t(3, 0, 2.5) Intercept                    sigma                
#>   source
#>  default
#>  default
#>  default
#>  default

# Replace a default, or append a prior for another parameter.
priors <- c(
  cogmod_priors(f, d),
  brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
  replace = TRUE
)
```
