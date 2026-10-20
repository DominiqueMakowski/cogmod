# Starting values that keep the sampler out of the flat regions

Builds an `init` argument for
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), for
the custom families whose default starting point is a bad one.

**This is not a tuning knob to reach for when sampling looks bad.** For
the families with a non-decision time the default start is actively
harmful, and the failure is silent: the chain does not error, it simply
never moves.

## Usage

``` r
cogmod_inits(formula = NULL, data = NULL, jitter = NULL, warmstart = NULL, ...)
```

## Arguments

- formula:

  The model formula, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Must carry the family, i.e. be built with
  `brms::bf(..., family = cogmod_gamma())`. May be left `NULL` only when
  `warmstart` is a `brmsfit`, whose formula is then used.

- data:

  The data, as passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html). May
  be left `NULL` only when `warmstart` is a `brmsfit`, whose data are
  then used.

- jitter:

  SD of the noise added on the unconstrained scale, so that chains start
  at different points. One number sets the population-level blocks
  (intercepts and slopes; for a slope, divided by how far its design
  columns reach) and a fifth of it the group-level and smooth blocks
  (`z_*`, `zs_*`, `sd_*`, `sds_*`), where a unit of noise is multiplied
  through a scale and a design column, once per group or basis function.
  Two numbers set the two tiers, population first. `0` gives identical
  starts. `NULL` (the default) means 0.25, or 0.05 with a `warmstart`,
  whose values come from a converged posterior.

- warmstart:

  A previous fit to start from instead of the family's generic values: a
  `brmsfit`, a
  [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
  object, the data frame
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) makes
  of one, or the path to a CSV file of it. The starting values are its
  posterior means, mapped onto this model by parameter name (a pilot on
  fewer participants included; see
  [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md));
  whatever it has no value for keeps the value this function would give
  it anyway.

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

- **`ndt` too large.** `brms` initialises on the *unconstrained* scale,
  so `init = 0` puts `ndt` at `exp(0) = 1` second, above nearly every
  sub-second response. Every response is then attributed to the outlier
  component, the decision parameters drop out of the density, and their
  gradient is exactly zero.

- **Shape below 1.** For
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
  and
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
  the density is unbounded at `ndt` whenever the shape is below 1, and
  `init = 0` starts a `softplus` shape at `log(2) = 0.69`, inside that
  region.

No single scalar `init` avoids both - `ndt = exp(c)` wants `c` around
`-1.6`, `shape = softplus(c)` wants it above `1.9` - so the parameters
have to be set separately, which is what this function does.

## How it works

The parameter declarations are read off
[`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
for the model you are fitting, and their dimensions off
[`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html),
so `0 + Intercept`, interactions, group-level terms and smooths are
matched however `brms` names and sizes them. **Every** declared
parameter gets a value, which also spares you CmdStan's message about a
partial list:

- **The family's values** go to the intercept of each distributional
  parameter, on its link scale - through the links on the family, so
  `cogmod_gamma(link_mu = "log")` is honoured - and, on the natural
  scale, to any parameter left out of the formula, which `brms` declares
  as a plain auxiliary parameter. Under `0 + Intercept` the value goes
  to the coefficient named `Intercept`.

- **`ndt`** starts at half the first percentile of the observed response
  times: below essentially every response, but on the data's own scale.
  The two errors are not symmetric: a start too small only means the
  shift has to grow, which the gradient does, while one too large
  removes the gradient.

- **The error accumulator** of
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  and
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  starts slower than the correct one (`driftone = 1` against `mu = 3`,
  `nuone = 0` against `mu = 0.7`): few data sets are a 50/50 race, and
  for the Wald and the LBA a start too fast costs far more than one too
  slow. Each family's page has the details.

- **Everything else** gets a generic value at least as good as Stan's
  own `U(-2, 2)`: zero for slopes and standardized effects, just above
  the bound for a positive parameter, the midpoint for a doubly bounded
  one, the identity for a Cholesky factor. A smooth's SD (`sds`) starts
  where its smooth is nearly flat, at 0.1 link units of wiggle, a
  different value on every basis; see *Smooth terms* in
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md).

Each chain gets a different draw. The noise is added on the
unconstrained scale, so a value can never leave its bounds. For a
regression coefficient it is divided by the largest row norm of the
coefficient block's design columns, wherever that exceeds 1, so that at
no row does the SD of what reaches the linear predictor exceed `jitter`:
a slope on age in years would otherwise start with an SD of six link
units at the oldest participant.

## Supported families

The family is read off `formula`, so build it with
`brms::bf(..., family = cogmod_gamma())`. The families are those of
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md),
listed in its *Supported families* section: every family with an `ndt` +
`poutlier` component;
[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
and
[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md),
whose parameters are lengths of time that the default start puts at 0 s
or 0.69 s; and the rating-scale families
[`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
[`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md)
and
[`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md),
where it describes a scale nobody uses - half the responses on the
midpoint or at an endpoint, or a U-shaped Beta. Any other family is an
error.

## See also

[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md),
[`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md),
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)

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
#> [1] -2.349416
#> 
#> $Intercept_poutlier
#> [1] -3.861294
#> 

# The rating-scale families are covered too. `pmid` starts at 0.05 rather
# than the logit origin's 0.5, which would put half of every response
# exactly on the midpoint of the scale.
r <- data.frame(y = rcogmod_choco(50, pmid = 0.05))
g <- brms::bf(y ~ 1, pmid ~ 1, family = cogmod_choco())
cogmod_inits(g, r, jitter = 0)(1)
#> $Intercept
#> [1] 0
#> 
#> $confright
#> [1] 0.5
#> 
#> $precright
#> [1] 2
#> 
#> $confleft
#> [1] 0.5
#> 
#> $precleft
#> [1] 2
#> 
#> $pex
#> [1] 0.1
#> 
#> $bex
#> [1] 0.5
#> 
#> $Intercept_pmid
#> [1] -2.944439
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
