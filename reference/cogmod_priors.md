# Priors that make a cogmod posterior proper

Fills in weakly informative priors for the parameters `brms` would
otherwise leave flat, or give a default meant for another
parameterization, for the model you are actually fitting.

**For the families with a non-decision time this is not a convenience.**
`brms` gives the intercept of a custom-family parameter it does not
recognise a flat, improper prior, and on the `ndt` + `poutlier` mixture
the likelihood is exactly flat in two directions: `poutlier` toward 1,
where every response belongs to the outlier component and the decision
parameters drop out of the density, and `ndt` toward 0, where the
gradient with respect to `log(ndt)` vanishes. Flat prior plus infinite
flat region is an improper posterior, and the fit does not fail loudly:
it returns intercepts around `1e14` with `Rhat` near 2 and an effective
sample size of about 5. Several families have flat directions of their
own (see *Why these parameters*).

## Usage

``` r
cogmod_priors(formula, data, ..., warmstart = NULL, prior_scale = 3)
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

- warmstart:

  Optional. A previous fit to centre the priors on: a `brmsfit`, a
  `cogmod_warmstart` object, its data frame, or the path to a CSV of it,
  as
  [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
  takes. `NULL` (the default) leaves the priors where this function put
  them. **Changes the posterior** - see *Centring the priors on a
  previous fit*.

- prior_scale:

  The prior SD, as a multiple of the source's posterior SD. Only used
  with `warmstart`. Larger is weaker; 1 would use the source's posterior
  as the prior.

## Value

A `brmsprior` object, to pass to `brms::brm(prior = )`.

## Details

The function starts from
[`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
for the model in hand, edits the rows it has a better prior for, and
returns the result of
[`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html):
a `brmsprior` object, printed as a table that includes the defaults
`brms` will use. Because every row comes from the model `brms` will
build, no prior can miss its parameter - `0 + Intercept`, interactions,
group-level terms and smooths need no special handling - and a malformed
specification errors here rather than inside
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html).

## Setting your own priors

Combine the result with
[`brms::prior()`](https://paulbuerkner.com/brms/reference/set_prior.html)
entries. `replace = TRUE` replaces one of these defaults; without it the
entry is added for another slot:

    priors <- c(
      cogmod_priors(f, df),
      brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
      replace = TRUE
    )

## Supported families

The family is read off `formula`, so build it with
`brms::bf(..., family = cogmod_lognormal())`. Covered: every family with
an `ndt` + `poutlier` component -
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
and the choice-and-RT families
[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
and
[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) -
plus
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md),
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md)
and the rating-scale families
[`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
[`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md)
and
[`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md).
Any other family, or a formula carrying none, gets a message and
`brms`'s defaults unchanged, so the call is safe to leave in a script.

## What gets set

For each parameter in the table under *The priors*:

- its **intercept** gets the prior in the table;

- its **slopes** get `normal(0, 0.2)` on the link scale, or the wider
  prior in the table. On a log or logit link a flat slope prior is not
  harmless: it can carry a condition into the flat regions while the
  intercept stays well behaved;

- its **group-level SDs** get `exponential(1)`, for the same reason, and
  the correlations keep `brms`'s LKJ;

- its **smooth terms** follow *Smooth terms*, as every smooth does.

Parameters not in the table, such as
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)'s
`sigma`, keep `brms`'s defaults, and so does the response's own `mu`:
its slopes are the effects being estimated, and `brms`'s `student_t()`
default is a fair intercept prior on most links. Two exceptions.
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
and
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md)
give `mu` an intercept prior, because there it is a time in seconds on
an `identity` link. And in the race families,
[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
and
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md),
`mu` is the rate of one accumulator, so its group-level SD gets the
`exponential(1)` the other accumulator's gets and swapping the `dec()`
coding does not change what the priors say.

A row `brms` left flat is filled. A row `brms` already filled is
replaced only where its default was written for a parameter of the same
name in one of its own families, and is wrong here: `uniform(0, min_Y)`
on an `ndt` left out of
[`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html), the
min-RT bound the outlier component exists to remove; `gamma(0.01, 0.01)`
on
[`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)'s
`shape`, which must be free to go negative; and `student_t(3, 0, 2.5)`
and the like on the `sigma` of
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
and
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md),
the `shape` of
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md)
and the `phi` of
[`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md)
and
[`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md) -
a Gaussian SD wider than most whole RT distributions, or on
[`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md)'s
`log` link a precision reaching 2853.

## ndt and poutlier

`ndt` is in **seconds**, like everything else in these families.
`normal(-1.2, 0.5)` on its `log` link centres it on 0.30 s with 95% of
the mass between 0.11 and 0.80 s: wide enough for slow populations and
demanding responses, and still a fence against `ndt -> 0`.

`poutlier` is a proportion. `normal(-6, 1)` on its `logit` link centres
it on 0.25%, with 95% between 0.035% and 1.7%, which is where the
empirical estimates sit: Miller (2024) puts the outlier proportion below
0.5% across four lexical-decision megastudies, and the untrimmed data of
the [outliers
article](https://dominiquemakowski.github.io/cogmod/articles/outliers.html)
agree. If the data were trimmed before fitting, tighten the prior
(`normal(-7, 0.5)`) rather than fixing `poutlier = 0`, which brings the
hard min-RT boundary back; see the *Trimmed data* section of
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md).

## The priors

On the link scale for a parameter in
[`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) - for
its intercept, or under `0 + Intercept` its coefficient named
`Intercept` - and on the natural scale for one left out (see *Parameters
left out of the formula*).

|  |  |  |  |  |
|----|----|----|----|----|
| Parameter | Families | In [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) | Left out | Slopes |
| `ndt` | all with `ndt` | `normal(-1.2, 0.5)` | `lognormal(-1.2, 0.5)` | `normal(0, 0.2)` |
| `poutlier` | all with `ndt` | `normal(-6, 1)` | `exponential(280)` | `normal(0, 0.2)` |
| `sigmabias` | [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md), [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md) | `normal(0, 1)` | `lognormal(-0.35, 0.75)` | `normal(0, 0.5)` |
| `dof` | [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md) | `normal(1.8, 0.7)` | `lognormal(1.8, 0.7)` | `normal(0, 0.3)` |
| `shape` | [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md) | `normal(0, 0.5)` | `normal(0, 0.5)` | `normal(0, 0.2)` |
| `boundary` | [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md), [`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md), [`cogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md) | `normal(0.5, 1)` | `lognormal(0, 0.6)` | `normal(0, 0.5)` |
| `sigmadrift` | [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
| `sigmandt` | [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md), [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) | `normal(-3, 1)` | `lognormal(-3, 1)` | `normal(0, 0.5)` |
| `tau` | [`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md), [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md), [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md) | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` | `normal(0, 0.5)` |
| `sigmabias`, `boundary` | [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md), [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md), [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
| `sigmazero`, `sigmaone` | [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md), [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
| `nuone` | [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md) | `normal(0.7, 1.5)` | `normal(0.7, 1.5)` | `normal(0, 0.5)` |
| `driftone` | [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md) | `normal(3, 2)` | `lognormal(1, 0.75)` | `normal(0, 2)` |
| `driftone` | [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md) | `normal(1, 2)` | `normal(1, 2)` | `normal(0, 1.5)` |
| `boundary` | [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) | `normal(1, 1)` | `lognormal(0.27, 0.5)` | `normal(0, 0.5)` |
| `bias` | [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) | `logistic(0, 1)` | `brms`'s `beta(1, 1)` | `normal(0, 0.5)` |
| `sigmadrift` | [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) | `normal(0, 1)` | `lognormal(-1, 0.75)` | `normal(0, 0.5)` |
| `sigmabias` | [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) | `normal(-2, 1)` | `beta(1, 5)` | `normal(0, 0.5)` |
| `mu` | [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md), [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md) | `normal(0.4, 0.25)` | \- | `brms`'s |
| `sigma` | [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md), [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md) | `normal(-2.3, 0.7)` | `lognormal(-2.3, 0.7)` | `normal(0, 0.5)` |
| `shape` | [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md) | `normal(0, 0.5)` | `lognormal(0, 0.5)` | `normal(0, 0.5)` |
| `confright`, `confleft` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md) | `normal(0, 1)` | `beta(2, 2)` | `normal(0, 0.5)` |
| `precright`, `precleft` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md) | `normal(2, 1.5)` | `lognormal(0.7, 0.7)` | `normal(0, 0.5)` |
| `phi` | [`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md) | `normal(2, 1.5)` | `lognormal(0.7, 0.7)` | `normal(0, 0.5)` |
| `phi` | [`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md) | `normal(0.7, 0.8)` | `lognormal(0.7, 0.8)` | `normal(0, 0.5)` |
| `pex` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md), [`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md) | `normal(-2, 1)` | `beta(2, 12)` | `normal(0, 0.5)` |
| `bex` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md), [`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md) | `normal(0, 1)` | `beta(2, 2)` | `normal(0, 0.5)` |
| `pmid` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md) | `normal(-2.5, 1)` | `exponential(9)` | `normal(0, 0.5)` |
| `pzero` | [`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md) | `normal(-2.5, 1)` | `exponential(9)` | `normal(0, 0.5)` |

## Parameters left out of the formula

Writing `ndt ~ 1` and leaving `ndt` out are not the same model to
`brms`. A parameter in
[`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html), even
as `~ 1`, gets a linear predictor and lives on its **link** scale,
reported as `ndt_Intercept` under *Regression Coefficients*. One left
out is declared as a plain auxiliary parameter, like `sigma` in
[`gaussian()`](https://rdrr.io/r/stats/family.html), and lives on its
**natural** scale, under *Further Distributional Parameters*, so it
needs a prior on that scale. The two columns of the table describe the
same belief wherever they can: exactly on a `log` or `identity` link,
where `lognormal` is `normal` on the log scale, and matched on the
median elsewhere.

The exception is a probability whose omission is itself a statement.
Leaving `poutlier` out says the data were trimmed or are not expected to
contain outliers, and leaving `pmid` or `pzero` out says the event is
not expected. Those get an exponential, with its **mode at zero**, which
no logit-scale prior can have, and the same median: `exponential(280)`
against `plogis(-6) = 0.0025`, and `exponential(9)` against
`plogis(-2.5) = 0.076`. It is still a prior, so a genuine spike of fast
responses still pulls the rate up.

## Why these parameters

Each family's own page has the full account. In short, a parameter is in
the table for one of five reasons:

- **Its link puts a limit at infinity, and the likelihood is flat long
  before it** - the `ndt -> 0` failure again. That is the start-point
  ranges (`sigmabias`) and the LBA and RDM thresholds that share a ridge
  with them (`boundary`), the between-trial variabilities (`sigmadrift`,
  `sigmandt`), and
  [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md)'s
  `dof`, which reaches the LogNormal only as `dof -> Inf`. Where the
  design cannot identify a range or a variability, fixing it in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  (`sigmabias = 0`) is better than estimating it behind a prior.

- **An accumulator that stops winning.** Push a race accumulator's rate
  down far enough and the likelihood stops depending on it, or on its
  SD: `nuone` and the SDs in
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  `driftone` in
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md).
  Under `dec(Error)` that is the error accumulator, and wherever errors
  are rare its posterior sits in the lower tail of the prior: the fence
  working, not a prior at odds with the data. `mu` has the mirror
  plateau; if the rarely chosen option is on `mu`, copy the prior onto
  it.

- **An arbitrary scale.**
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)'s
  likelihood is exactly constant along a ray that rescales every
  evidence parameter. The priors on `sigmazero` and `sigmaone` make the
  posterior proper; only fixing one of them in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  (`sigmazero = 1`) identifies the scale.

- **A ridge with another parameter**: `tau` with `ndt` in
  [`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md),
  `shape` with `mu` in
  [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md),
  `dof` with `poutlier` at the fast end of
  [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md);
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)'s
  `shape`, whose prior also keeps the sampler clear of
  `sigma * shape >= 1`, where the density is unbounded at the shift.

- **Identified, but `brms` leaves it flat or gives it a default for
  another parameterization**: the Wald and DDM thresholds, the DDM's
  `bias`, and the ex-Gaussian and rating-scale parameters (below).

No prior is set on the shape of
[`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
or
[`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
although a shape below 2 makes their `ndt` gradient unbounded: the
likelihood itself prefers that region, by around 100 log units on the
data in the [RT models
article](https://dominiquemakowski.github.io/cogmod/articles/rt_models.html),
so a prior weak enough to be a default could only bias it.
[`?rcogmod_weibull`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
sets out what to do instead.

## The ex-Gaussian and the GEG

[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
has no `ndt`, but `mu`, `sigma` and `tau` are all lengths of time in
seconds, which `brms` cannot know. On the `softplus` link, `sigma` (the
Gaussian SD) gets 25 to 330 ms with a median of 96 ms, and `tau` (the
mean of the exponential tail) 55 to 630 ms with a median of 201 ms: the
range they occupy across the usual simple- and choice-RT tasks. `mu` is
on an `identity` link, where `brms`'s default, a `student_t()` with a
scale of 2.5 s, spreads the Gaussian's centre over several seconds
either way, and under `0 + Intercept` it gets no prior at all.
`normal(0.4, 0.25)` puts 95% of it between -0.09 and 0.89 s and
deliberately allows negative values: for fast, heavily tailed data the
Gaussian component belongs near or below zero, with `tau` carrying the
mass. It is the centre of the Gaussian alone, so the mean of the
distribution is `mu + tau`.
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md)
takes the same rows, and a `shape` centred on the ex-Gaussian it nests
(`shape = 1`, with 95% between 0.37 and 2.7): `shape` and `mu` sit on a
ridge, and the prior holds it.

## The rating-scale families

[`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
[`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md)
and
[`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md)
have no infinite flat region of the RT kind, but every parameter arrives
either flat or with a `brms` default meant for another parameterization.

- **Point masses.** `pmid` is the probability of exactly the midpoint,
  and `pzero` of an extra category outside the scale. When the event is
  absent from the data the likelihood rises all the way to zero, and
  under a flat prior the logit runs off to minus infinity: `poutlier`'s
  failure, and `poutlier`'s treatment. `normal(-2.5, 1)` centres it on
  8%, with 95% below 37%. To switch one off, fix it in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  (`pmid = 0`); a response exactly on the midpoint is then impossible,
  and the fit will not start.

- **Extremes.** `pex`, the total probability of a 0 or a 1, is *not*
  given a mode at zero: extreme responding is a response style, not a
  contaminant. `normal(-2, 1)` centres it on 12%, with 95% between 2%
  and 49%. `bex`, which end the extremes favour, is centred on
  symmetric, which also keeps it away from 0 and 1, where one gate
  closes and every response at that end has zero probability.

- **Beta shapes.** `confright` and `confleft`, the Beta means of each
  side, are centred on the middle of their half. The precisions are the
  awkward ones: below about 1 the Beta is U-shaped and unbounded at both
  ends, and a large one concentrates it sharply. `normal(2, 1.5)` on the
  `softplus` link spans 0.33 to 4.95, leaving 17% of its mass below 1,
  since a polarising item really can produce a U.
  [`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md)'s
  `phi` is on a `log` link: `normal(0.7, 0.8)`, 0.42 to 9.7 with a
  median of 2.

`mu` is left to `brms`, whose `student_t(3, 0, 2.5)` is the standard
choice on a `logit` link.

## Smooth terms

A smooth's SD, `sds`, is not in link units. `brms` writes the penalised
part of an [`s()`](https://paulbuerkner.com/brms/reference/s.html) or
[`t2()`](https://paulbuerkner.com/brms/reference/s.html) term as
`sds * Zs %*% zs` with `zs ~ N(0, 1)`, and how far one unit of `sds`
moves the linear predictor depends on the basis `mgcv` builds: about
0.45 link units for `s(x)`, 6 for `s(x, bs = "cr")` and 0.1 for
`t2(x, z, k = c(5, 5), bs = c("cr", "cr"))`, almost regardless of the
data. So each smooth term gets `exponential(rate)`, with the rate set to
that term's own scale (read off
[`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)):
`exponential(1)` on the wiggle measured in link units, median 0.69 and
95% below 3, whatever the basis. A term whose scale cannot be computed
gets `exponential(1)`. This applies to every smooth, the response's
included, and it matters however many rows there are: a smooth's SD is
informed by the handful of coefficients behind each penalty, as a
group-level SD is by the number of groups.

The unpenalised part of a smooth - the linear trend of an `s(x)`, the
marginal trends of a
[`t2()`](https://paulbuerkner.com/brms/reference/s.html) - is reported
under class `b`, where it would take its parameter's slope prior; but
`mgcv` sets the scale of those columns, not the user. On the parameters
in the table, each such column gets `normal(0, s / rms)` instead, which
is the slope prior `normal(0, s)` stated on the trend at a typical row,
`rms` being the column's root mean square.

## Centring the priors on a previous fit

`warmstart` takes what
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
takes - a `brmsfit`, the table
[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) makes of
one, or the path to a CSV of it - and replaces each prior with
`normal(median, prior_scale * sd)`, from that fit's posterior for the
same parameter. It applies to the population-level intercepts and
coefficients, the group-level SDs and the parameters left out of the
formula; correlations keep their LKJ. A parameter the source never saw
keeps the prior it would get anyway, and
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)'s
[`print()`](https://rdrr.io/r/base/print.html) says how many there are.

**This changes the posterior.** Unlike
[`cogmod_inv_metric()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
and
[`cogmod_step_size()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md),
which only change how the sampler moves, a prior is part of the model.
**If the source was fitted to data the new model also contains, it
counts those data twice**: a pilot on 5 of 10 participants used to
centre the priors for all 10 uses those 5 twice, and the intervals come
out too narrow by an amount nothing in the output shows. Use it when the
source is an independent data set - last year's sample, another lab's,
another session - or for a deliberately sequential analysis.

`prior_scale = 1` makes the source's posterior the prior. The default of
3 has nine times its variance: enough to keep the sampler in the right
region and out of the flat directions, weak enough for data that
disagree to win.

    # last year's sample as the prior for this year's
    priors <- cogmod_priors(f, df, warmstart = fit_2025)
    priors <- cogmod_priors(f, df, warmstart = fit_2025, prior_scale = 6)  # weaker

## References

- Miller, J. (2024). Estimating the proportions and latencies of
  reaction time outliers: A pooling method and case study of lexical
  decision tasks. *Behavior Research Methods*, *56*(7), 7280-7306.
  [doi:10.3758/s13428-024-02419-y](https://doi.org/10.3758/s13428-024-02419-y)

## See also

[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md),
[`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md),
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)

## Examples

``` r
d <- data.frame(RT = rcogmod_lognormal(50, ndt = 0.3, poutlier = 0.02))
f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
  family = cogmod_lognormal()
)
cogmod_priors(f, d)
#>                   prior     class coef group resp     dpar nlpar lb ub tag
#>  student_t(3, 0.8, 2.5) Intercept                                         
#>       normal(-1.2, 0.5) Intercept                      ndt                
#>           normal(-6, 1) Intercept                 poutlier                
#>    student_t(3, 0, 2.5) Intercept                    sigma                
#>  lognormal(-0.35, 0.75) sigmabias                                 0       
#>   source
#>  default
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
