# cogmod 0.2.1

## Breaking changes

* The **outlier component is now a half Normal with a fixed scale of 0.2 s**,
  where it was a half Student-t with 3 degrees of freedom and a user-supplied
  scale. Two things changed, for one reason.

  The Student-t's tail was heavier than every decision density in the package,
  so far-out slow responses were eventually explained better by the outlier
  component than by the model: against a shifted LogNormal at `poutlier =
  0.02`, a 5 s response was attributed to it with probability 0.86 and the
  crossover sat at 3.86 s, with `ndt` pulled up behind it. `vignette("outliers")`
  already flagged this as a defect. A Gaussian never gets there - the same
  responsibility is 0.000 out to 30 s - and it costs nothing where the
  component is actually needed, because `exp(-x^2 / 2s^2)` kills the far tail
  at any scale: at 0.2 s it holds 76% of its peak density at 0.15 s and 46% at
  0.25 s, against 85% and 66% for the half-t. The slow tail now belongs to the
  decision family, which is what `cogmod_loggamma()`'s `shape` and
  `cogmod_invgaussian()`'s `sigmadrift` are for.

  A welcome side effect: the `poutlier -> 1` degenerate mode is now thousands
  of log-likelihood units below the sensible one rather than hundreds, and
  `ndt` and the decision parameters no longer drop out of the density there,
  because a half Normal cannot explain a slow response at all. The mode still
  has infinite volume in `poutlier` itself, so `cogmod_priors()` is still not
  optional.

* **`minrt` is removed** from every family, density, RNG, `*_stanvars()` and
  `*_lpdf_expose()`. The package works in **seconds**, full stop. The
  equivariance `minrt` bought in the likelihood was already fictional end to
  end: `cogmod_priors()` shifted only the `ndt` prior with it, while the
  `sigmandt` prior of `cogmod_ddm()`, the `sigmadrift` prior of
  `cogmod_invgaussian()` and the `mu` priors are stated in seconds outright -
  and `cogmod_priors()` is not optional. Making the assumption explicit costs
  one argument from about twenty signatures and removes a whole class of
  misconfiguration.

  Calls passing `minrt` now fail with R's usual `unused argument` error. Data
  in another unit fails **silently**, as it always did when `minrt` was left at
  its default: the outlier component's log-density at `RT = 400` is about
  `-2e6`, so it contributes nothing, the mixture collapses to the unmixed
  shifted family, `poutlier` goes to zero and `ndt` is pinned by the fastest
  observed response. Nothing errors and the chains still initialise. Divide by
  1000 before fitting.

  `cogmod_priors()` accordingly gives `ndt` a fixed `normal(-1.2, 0.2)`, and
  `cogmod_inits()` starts it at 0.1 s.

* `cogmod_invgaussian()` gains **`sigmadrift`**, the between-trial SD of the
  drift rate, so the Wald can produce the long right tails empirical RT
  distributions have. Described in full under 0.2.0 below.

## Performance

* **`cogmod_rdm()` samples about 1.5x faster.** It was the most expensive
  likelihood in the package - roughly 5.9 us per observation per gradient,
  against 2.6 for `cogmod_lba2()` and 0.8 for `cogmod_lnr()` - and two thirds of
  that was avoidable. The Stan survival evaluated the normal CDF at `alpha` and
  at `beta` twice each, once inside `log_g()` and once again for the terms that
  follow it, and it assembled its six signed terms through a helper returning a
  `vector[2]`, which allocates on the autodiff stack once per term per
  observation per leapfrog step.

  `log_g()` now takes `log Phi(u)` from its caller, which brings the survival
  from six normal-CDF evaluations to four, and the six terms are grouped into
  the two differences that are individually monotone in the threshold, so each
  is one `log_diff_exp()` of known sign and nothing on the hot path returns a
  vector. The maths is unchanged and the grouped form cancels *less*: against
  the R implementation it agrees to 6e-11 where the term-by-term form reached
  5e-9. Gradients are unchanged to the precision finite differences can resolve.

* **`posterior_predict()` is about 2x faster for the choice+RT families.** brms
  calls it once per observation, so anything done per call is paid thousands of
  times, and for every family except `cogmod_ddm()` the sampling itself was the
  small part: two `data.frame()` constructions and an `as.matrix()` came to
  roughly 63% of the call, against 13% for the actual draws. The registry's
  `rng` entries now return a bare `list(rt, response)` and `.rchoice()` can
  return a matrix directly, so the prediction path builds no data frame at all.
  A `cbind()` costs about a seventeenth of the `data.frame()` it replaces.

  `rcogmod_lnr()`, `rcogmod_rdm()`, `rcogmod_lba2()` and `rcogmod_ddm()` still
  return a data frame, with the same draws for the same seed - only the internal
  path changed. Measured on the `vignette("decision_making")` models at 20
  draws: LNR 41.5 -> 20.7 us per observation-draw, LBA 32.9.

  `cogmod_ddm()` gains least (560 -> 474 us) because its cost is elsewhere:
  `brms::rwiener()` samples one draw per call, and roughly 85% of that is
  per-call setup which cannot amortise when every posterior draw has its own
  parameters. Predicting from a DDM remains an order of magnitude dearer than
  from the races.

# cogmod 0.2.1

## New features

* New family **`cogmod_exwald()`** ([Schwarz, 2001](https://doi.org/10.3758/bf03195403)):
  the decision time is a Wald convolved with an exponential residual stage of
  mean `tau` - the mechanistic counterpart of `cogmod_exgaussian()`, whose first
  stage is a descriptive Gaussian instead, with `tau` meaning the same thing in
  both. The mean exists and is `ndt + boundary / mu + tau`, so
  `posterior_epred()` returns a number.

  The density has two branches, **both exact**. Where `mu^2 > 2 / tau` the
  convolution collapses to a closed form in the Wald CDF; below that - which is
  the common regime, since at a drift of 3 and a threshold of 0.5 the closed form
  needs `tau > 0.22 s` - the same expression continues analytically through the
  Faddeeva function, giving
  `g * exp(-(boundary - mu * t)^2 / (2 * t)) * Re[w(z)]`. The exponent is the
  Wald's own, so nothing overflows, and the branches meet exactly at
  `mu^2 = 2 / tau`. Across a grid spanning the usual RT region the density
  integrates to 1 to within 8e-12, the mean is right to 3e-9, and the relative
  step across the branch seam is 5e-8.

  Note there is deliberately no `sigmadrift`: it and `tau` both fatten the right
  tail and are very hard to tell apart, and `cogmod_invgaussian()` is where the
  drift-variability route lives. `ndt` and `tau` also share a ridge - both delay
  the response, and only the shape of the leading edge separates them - so
  `cogmod_priors()` gives `tau` the same `normal(-1.5, 0.7)` it gives
  `cogmod_exgaussian()`. Fixing `ndt = 0` in `bf()` recovers Schwarz's own model.

* New family **`cogmod_bisa()`**, the Birnbaum-Saunders or fatigue-life
  distribution ([Birnbaum & Saunders,
  1969](https://doi.org/10.2307/3212003)): a first-passage-time model in which
  evidence arrives in **discrete cycles and only ever towards the boundary** -
  what is random is the size of each increment, never its sign. Summing those
  increments and applying the CLT, then treating the cycle count as continuous,
  gives the first-crossing time.

  It is parameterized mechanistically, as `mu` (drift) and `boundary`
  (threshold), so it sits directly alongside `cogmod_invgaussian()` with the
  parameters meaning the same thing and only the mode of accumulation differing.
  Fixing the per-cycle SD at 1 is the same convention that fixes the Wald's
  diffusion coefficient, and it makes `(mu * t - boundary) / sqrt(t)` **exactly**
  standard normal - the usual `(1 / a) * (sqrt(t / b) - sqrt(b / t))` written in
  these parameters, with `b = boundary / mu` and `a = 1 / sqrt(mu * boundary)`.
  The map between the two is a bijection, so nothing is given up.

  Everything is then closed form, and the density is the Wald's own tilted by
  `(mu * t + boundary) / (2 * boundary)` - one sign apart from it. That tilt
  makes the family an **equal mixture of an inverse Gaussian and its
  length-biased twin**, so at the same `(mu, boundary)` it is slower and more
  dispersed than the Wald (mean 0.222 s against 0.167, SD 0.184 against 0.136 at
  `mu = 3, boundary = 0.5`), while keeping the same exponential-order right
  tail. `E[T] = ndt + boundary / mu + 1 / (2 * mu^2)` is always finite, so
  `posterior_epred()` returns a number, and the median is exactly
  `ndt + boundary / mu`. There is no `sigmadrift`: the extra dispersion comes
  from the mixture structure at no cost in parameters, and drift variability is
  what `cogmod_invgaussian()` is for.

  It is also the cheapest first-passage density in the package - one log and one
  square, no branch and no special function - and `rcogmod_bisa()` is one normal
  draw per observation, exact, with no rejection step.

* New family **`cogmod_logstudent()`**: `log(RT - ndt)` follows a Student-t, a
  robust LogNormal that varies kurtosis where `cogmod_loggamma()` varies skew.
  The heavy tail absorbs slow contaminants into the likelihood rather than into
  a mixture component, which matters because the `poutlier` component is a half
  Normal and by construction cannot explain a slow response. At `dof = 5` the
  probability of a decision time beyond 5 s is about five orders of magnitude
  larger than the matching LogNormal's.

  The degrees of freedom are called **`dof`**, not `nu`: `cogmod_lnr()` already
  spends `nuzero`/`nuone` on drift rates, and `brms` recognises the name `nu`
  and supplies defaults of its own for it.

  Two things to know. **The mean does not exist** for any finite `dof`, so
  `posterior_epred()` errors rather than returning a number; the median is exact
  at `ndt + exp(mu)`. And **the density is unbounded at `ndt`** - integrable, so
  the posterior stays proper, but the likelihood has no maximum, which is one
  more reason `cogmod_priors()` is not optional. A Student-t is also symmetric
  on the log scale, so a small `dof` fattens the fast tail as well as the slow
  one and competes with `poutlier`; `cogmod_priors()` centres `dof` at 6 with
  95% of its mass between 1.5 and 24 to keep that in check.

* `cogmod_priors()` now supports **`cogmod_exgaussian()`**, where before it
  returned the `brms` defaults with a message. `sigma` and `tau` are both
  lengths of time in seconds behind a `softplus` link, which `brms` has no way
  to know: `tau` arrives flat, and `sigma` arrives with the
  `student_t(3, 0, 2.5)` that `brms` supplies because it recognises the *name* -
  a Gaussian SD centred on 0.69 s modelled, 1.9 s omitted, wider than most whole
  RT distributions. They now get `normal(-2.3, 0.7)` and `normal(-1.5, 0.7)` on
  the link scale (roughly 25-330 ms for `sigma`, 55-630 ms for `tau`), and the
  matching `lognormal` when the dpar is left out of `bf()` altogether. `mu` is
  deliberately untouched: it is the response's own intercept and the `brms`
  default is already proper and reasonable there.

* `cogmod_invgaussian()` gains **`sigmadrift`**, the between-trial SD of the
  drift rate. Above zero, each trial draws its own drift from a
  `Normal(mu, sigmadrift)` truncated at zero, which is what lets the Wald
  reach the long right tails empirical RT distributions have. Marginalising
  over that draw is a Gaussian integral, so the density stays closed form and
  costs two normal CDFs; `sigmadrift = 0` gives back the previous density
  exactly, not approximately.

  The truncation is what keeps the density proper: a single-boundary
  accumulator handed a negative drift never terminates, so an untruncated
  Normal would leave up to a third of the mass unaccounted for.
  `cogmod_ddm()`'s `sigmadrift` needs no such truncation, a diffusion between
  two boundaries always absorbing at one of them.

  It is fixed the same way as the `cogmod_ddm()` variability parameters -
  writing `sigmadrift = 0` in `bf()` pins it and recovers the classic Wald,
  while leaving it out of `bf()` *estimates* it. Fixing it is the better
  default: `sigmadrift` and `poutlier` both fatten the right tail and are only
  weakly distinguishable, and `cogmod_priors()` gives `sigmadrift` a
  deliberately informative prior where it is estimated. Note that with
  `sigmadrift > 0` the density decays as `t^-2` and **the mean does not
  exist**, so `posterior_epred()` returns `Inf`; summarise
  `posterior_predict()` draws instead.

  Two consequences for existing code. The `drift`/`boundary`/`ndt`/`poutlier`
  functions gained an argument, so `sigmadrift = 0` now sits between `ndt` and
  `poutlier` in the signatures of `rcogmod_invgaussian()`,
  `dcogmod_invgaussian()` and `pcogmod_invgaussian()` (positional calls that
  passed `poutlier` fourth need updating; named calls are unaffected). And a
  formula that does not mention `sigmadrift` at all now estimates it rather
  than fitting the fixed-drift Wald.

* New `cogmod_stanvars()`: the third of the three generics that take the model
  rather than the family, alongside `cogmod_priors()` and `cogmod_inits()`. It
  reads the family off the formula and returns that family's Stan code, so the
  family is named once - in `bf()` - instead of three times:

  ```r
  f <- bf(RT ~ Condition, ndt ~ Condition, family = cogmod_lognormal(minrt = 0.25))
  brm(f, data = df,
      prior    = cogmod_priors(f, df),
      init     = cogmod_inits(f, df),
      stanvars = cogmod_stanvars(f))
  ```

  This is not only tidier. `minrt` is baked into the generated Stan code as a
  literal, because a Stan function cannot see the data block, so
  `cogmod_lognormal(minrt = 0.25)` fitted with `cogmod_lognormal_stanvars()` runs
  happily against an outlier component the family does not describe.
  `cogmod_stanvars()` takes `minrt` off the family, and the two cannot
  disagree. The per-family `<family>_stanvars()` functions are unchanged.

* `cogmod_inits()` now supports `cogmod_exgaussian()`, whose three parameters are
  all on the RT scale behind a `softplus` link and so are equally badly served
  by the default start at `log(2) = 0.69` s - which makes the Gaussian SD alone
  wider than most whole RT distributions.

* `cogmod_inits()` now returns a value for **every** parameter the Stan program
  declares, not only the ones it has an opinion about, so CmdStan no longer
  prints `Init values were only set for a subset of parameters` and lists the
  rest. Regression slopes and standardized group-level effects start at zero,
  group-level and spline SDs just above zero, Cholesky factors at the identity;
  all are at least as good a starting point as Stan's own `U(-2, 2)`.

  Two related fixes come with it. Links are now read off the family rather than
  the registry, so `cogmod_gamma(link_mu = "log")` is honoured. And the jitter is
  applied on the unconstrained scale - additive when free, multiplicative for a
  positive parameter, on the logit scale for a bounded one - so a jittered
  start can no longer land outside its own bounds, which it could previously
  for a dpar left out of the formula and estimated on the natural scale.

* New `cogmod_inits()`: starting values for the families that estimate `ndt`
  directly. `brms` initialises on the unconstrained scale, so `init = 0` puts
  `ndt` at `exp(0) = 1` second - above most sub-second RTs, which leaves every
  response attributed to the outlier component and the decision parameters with
  no gradient at all. For `cogmod_gamma()` and `cogmod_weibull()` it *also* puts the
  shape at `softplus(0) = 0.69`, below the 1 at which the density becomes
  unbounded at the shift. No single scalar avoids both, since the two pull in
  opposite directions.

  On 1500 simulated Gamma trials (true shape 3, true `ndt` 0.25), `init = 0`
  left the shape stuck at its starting value with `Rhat` 2.3 and an ESS of 3
  after 306 s; `cogmod_inits()` recovered shape 3.23 and `ndt` 0.227 with
  `Rhat` 1.01 in 28 s. An informative prior on the shape did not rescue
  `init = 0` - a prior cannot move a chain whose gradient is zero.

  Parameter names are read off `brms::make_stancode()` for the model actually
  being fitted, so `0 + Intercept`, interactions, group-level terms and smooths
  are all handled; anything it does not recognise is left to Stan.

* New `cogmod_loggamma()` family: a shifted Log-Gamma model for reaction times,
  equivalently a shifted generalized gamma. `log(RT - ndt)` follows a
  location-scale log-gamma with location `mu`, scale `sigma` and shape `shape`,
  and `ndt` / `poutlier` / `minrt` work exactly as in `cogmod_lognormal()`.

  `shape` is unconstrained, with `shape = 0` in the interior: it recovers
  `cogmod_lognormal()` exactly, `shape = sigma` the shifted Gamma, `shape = 1` the shifted
  Weibull and `shape = -1` the shifted inverse Weibull. Fitting it is therefore a
  way of testing whether the LogNormal shape is adequate - an interval for `shape`
  covering 0 says it is. Negative `shape` gives a power-law right tail.

  Note the boundary at `sigma * shape >= 1`, where the decision density becomes
  unbounded at `ndt` and the likelihood with it; `cogmod_priors()` sets
  `normal(0, 0.5)` on the `shape` intercept to keep well clear of it.

  **Fit this family with `init = 0`.** The prior keeps the posterior clear of
  that boundary but not the starting point: `brms` initialises from `U(-2, 2)`
  on the unconstrained scale, so about 15% of chains start with
  `sigma * shape >= 1`, fall into the spike at `ndt` and never finish. `init = 0`
  starts every chain at `shape = 0`, the LogNormal, and removes the problem.

* `with_outliers()`, `without_outliers()`, `p_outlier()` and `cogmod_priors()`
  now work on `cogmod_loggamma()` as well as `cogmod_lognormal()`.

* `cogmod_stanvars()` now **warns when the evidence scale is left free** for
  `cogmod_lba1()` and `cogmod_lba2()`. Both have a likelihood that is *exactly*
  constant along the ray that multiplies the drift rates, their SDs, the
  start-point range and the threshold offset by a common factor - verified to
  machine precision, not merely near-flat - so nothing in the data can pick a
  point on it. The failure is quiet rather than loud: the fit converges,
  `pp_check()` looks right, and only the individual parameter estimates are
  meaningless, being whatever the priors happen to say about that direction.

  Fixing any one member of the ray in `bf()` pins it and silences the warning -
  `sigmazero = 1` conventionally, as `rtdists` and `EMC2` both do, but
  `boundary = 1` or `sigmabias = 0.5` work as well. Note that leaving a
  parameter *out* of `bf()` does not fix it: `brms` declares it as a free
  auxiliary parameter and the ray stays exactly as free, which is the case the
  warning mostly exists to catch. `cogmod_rdm()` and `cogmod_ddm()` are quiet
  by construction, their unit diffusion coefficient having pinned the scale
  already.

## Bug fixes

* `cogmod_ddm()` no longer reports `Non-finite gradient` during warmup or a
  Pathfinder search, and no longer collects the divergent transitions that come
  with it. Stan's classic 4-parameter `wiener_lpdf()` - much the fastest of the
  three Wiener densities Stan offers, and the one this family used whenever the
  three between-trial variability parameters were zero - returns `-inf` in two
  regions, and hands back **NaN** partial derivatives when it does. A NaN
  partial is not made harmless by the mixture weight on it being zero:
  reverse-mode multiplies the (zero) adjoint into the stored partial, and
  `0 * NaN` is `NaN`, so a single trial in one of those regions turns the
  gradient of the whole model to NaN, and Stan rejects the proposal.

  The two regions are the alternating small-time series losing its sum to
  cancellation - which depends only on the rescaled decision time
  `tau = t / boundary^2`, not on the drift or the scale separately, and sets in
  below `tau = 6.6e-4` - and the density underflowing to zero, which a cheap
  leading-term estimate detects. Those calls now go to Stan's `sv`-capable
  density instead, which works in log space throughout and stays finite, with
  finite gradients, down to log-densities of `-1e7`. The two agree to `1e-13`
  where the paths meet, so there is no step in the likelihood, and the fast path
  still handles the overwhelming majority of evaluations.

  On the 2000-trial fit in `?cogmod_ddm`, over three seeds: 16 `Non-finite
  gradient` reports per Pathfinder run became 0, and 47-420 divergent
  transitions per NUTS run became 0. Sampling is about 35% slower per iteration
  and roughly two to a hundred times *better* per effective sample.

  The same guard is applied to the general 7-parameter form, used when
  `sigmabias` or `sigmandt` is nonzero, which fails the same way. There it
  returns the `-inf` that form would have returned anyway, but as a constant,
  which carries no partial derivatives.

* The test suite runs in a third of the time (2472 s to 968 s on Windows).
  Every family's Stan `lpdf` now goes into **one** model, compiled once per
  session, instead of nine separate `*_lpdf_expose()` compilations; the
  factorial parameter sweeps are thinned to subsets that still cover every level
  of every factor; and the six `brms::brm()` model fits are behind
  `COGMOD_TEST_SLOW`, which the CI workflow sets. Setting it locally runs them:

  ```r
  Sys.setenv(COGMOD_TEST_SLOW = "true"); devtools::test()
  ```

  The shared model also removes a trap the RDM tests had worked around with a
  cache of their own: `expose_functions()` fails on a model `cmdstan_model()`
  returns pre-compiled, so a second call in the same session errors rather than
  reusing the first.

* `cogmod_priors()` now covers `cogmod_lnr()`'s `nuone`, `sigmazero` and
  `sigmaone`, which it previously left flat. Push an accumulator's rate down far
  enough and it stops finishing first ever; the density then depends on it only
  through the loser's survival term, which has already saturated at 1. Past
  about `nuone = -6` the log-likelihood is *exactly* constant, and that
  accumulator's `sigma` is unidentified along with it. The outlier component
  makes this reachable rather than hypothetical: it floors the trials the
  retreating accumulator can no longer explain, so the plateau is there even
  when both responses are observed.

  `nuone` gets `normal(0.7, 1.5)` on its identity link, the two sigmas
  `normal(0, 1)` on softplus (`lognormal(-0.7, 0.75)` when omitted from `bf()`
  and so on the natural scale), and all three `normal(0, 0.5)` on slopes. `mu` -
  which is `nuzero` - has the mirror-image plateau, but it is the response's own
  intercept and `brms` already gives it a proper `student_t` default, so it is
  left alone; if you model a rarely-chosen option it is worth mirroring the
  `nuone` prior onto it by hand.

  `cogmod_inits()` already covered this family and is unchanged.

* `cogmod_priors()` now covers `cogmod_lba1()`'s `sigmabias` and `boundary`,
  which it previously left on `brms`'s flat default. As the start-point range
  approaches zero the LBA converges smoothly to the recinormal, so the
  likelihood stops depending on `sigmabias` altogether - and a `softplus` link
  reaches zero only at minus infinity. That is a flat prior over an infinite
  flat region, the same improper posterior the function already exists to
  prevent for `ndt` and `poutlier`, and it failed just as quietly: on the
  4285-trial fit in `vignette("rt_models")`, `sigmabias` for one condition ran
  off to `softplus(-10.4) = 3e-05` with `Rhat` 1.69 and an effective sample size
  of 6, while every other parameter looked healthy. With the priors in place the
  same fit gives `Rhat` 1.02, an effective sample size of 387, no
  maximum-treedepth hits, and a finite estimate. `boundary` is covered too,
  since `b = sigmabias + boundary` puts the two on the same ridge.

  Both get `normal(0, 1)` on the link scale, `lognormal(-0.7, 0.75)` when the
  dpar is omitted from `bf()` and so lives on the natural scale, and
  `normal(0, 0.5)` on slopes - wider than the blanket `normal(0, 0.2)` the other
  dpars take, because the point is to fence off zero rather than to shrink
  effects. Families can now declare such rows in the registry, so the next one
  that needs them does not need a special case.

* The `?rcogmod_weibull` and `?rcogmod_gamma` notes about the shape were
  understated. They said the density is unbounded at `ndt` for a shape below 1,
  which is true, but the threshold that matters in practice is **2**: below it
  the derivative of the log-likelihood with respect to `ndt` is unbounded at
  every observation, so the posterior stays proper while the sampler grinds. On
  the data in `vignette("rt_models")` the Weibull shape comes out at 1.4, `ndt`
  lands inside the dense left edge of the data, the step size collapses to 0.005
  against 0.19 for `cogmod_lognormal()`, and mean treedepth goes from 3.9 to
  8.1 - 19x the gradient evaluations and 19x the wall time, with `Rhat` 1.18 on
  `ndt`. The density is cheap; all of the cost is geometry. Both help pages now
  set out the three regimes and say to prefer `cogmod_loggamma()`, which nests
  the Weibull at `shape = 1`, when the shape comes out below 2.

  No prior is set on those shapes, deliberately, and none is set on `ndt`
  either. Every obvious remedy was tried on that fit and measured:

  - `normal(2.4, 0.4)` on the shape (softplus scale, 95% of its mass above a
    shape of 1.9) moved the posterior shape by 0.01. The likelihood prefers the
    low-shape corner by around 100 log units; the prior contributes 5.
  - `normal(-1.25, 0.05)` on `ndt`, centred below the fastest bulk response,
    left the posterior 6.5 prior SDs away at essentially its unconstrained
    value. The `ndt` likelihood has a posterior SD of 0.003, some fifteen times
    sharper than that prior. The attempt cost 4% divergent transitions against
    0.5%, 16% of iterations at maximum treedepth against 7%, `Rhat` 1.43 against
    1.18, and a slightly worse `loo`.
  - Fixing `ndt` at the fastest observed response does remove the problem, by
    removing the parameter - and reinstates the min-RT bound this
    parameterization exists to remove. On these data the fastest response is
    71 ms, which is not a decision.

  Note what is *not* wrong: `ndt` and the shape are jointly identified, sharply
  (posterior SD on `ndt` of 3 ms), so this is not two parameters trading off
  with nothing to separate them and pinning one is not the missing ingredient.

  The slow sampling and a poor fit turn out to be the same fact. Across the ten
  families fitted in `vignette("rt_models")` the Weibull comes **last** by
  `loo` - 196 elpd (SE 21) behind `cogmod_loggamma()`, and 95 behind the next
  worst. What the sampler struggles with is the model contorting itself to
  represent a left edge it cannot otherwise reach. Where the shape comes out
  above 2 the family is fine, as `cogmod_gamma()` is on these same data at 2.2;
  a shape below 2 is best read as the model asking for a different one.

  Under the older `ndt = tau * min(RT)` parameterization the same singularity
  was damped by the logit Jacobian vanishing as `tau` approached 1, which is why
  it only became visible once `ndt` was estimated directly.

* `cogmod_priors()` no longer leaves `brms` warning that a global `b` prior
  "will not be used in the model as all related coefficients have individual
  priors already". It set both the blanket row and the row for the coefficient
  named `Intercept`, which is right when there are slopes beside the intercept
  but leaves the blanket row covering nothing under `ndt ~ 0 + Intercept`. The
  pair is now checked in both directions.

* `cogmod_lognormal()` and `cogmod_loggamma()` also carried hand-written R
  copies of the shared mixture - their own `r*()`, `d*()`, `.prepare_*()`,
  `log_lik_*()`, `posterior_predict_*()` and `posterior_epred_*()`, about 420
  lines reproducing what the other seven families delegate. All of it now
  delegates too, verified bit-identical to the code it replaced before the
  change was made.

* `cogmod_lognormal()` and `cogmod_loggamma()` generated their Stan code from
  hand-written copies of the shared mixture rather than from the registry every
  other shifted family uses. The copies had drifted: the LogNormal one described
  its outlier scale as `0.8 * minrt` when it is `minrt`, and neither carried the
  folded normalising constant that makes the outlier term ~1.4x cheaper per
  gradient evaluation. Both now come from the one generator, verified against
  the R density to machine precision.

* `cogmod_priors()` now also covers dpars left out of `bf()` entirely. `brms`
  declares those as plain auxiliary parameters - class `"<name>"` with an empty
  `dpar`, on the natural scale with no link - so matching on `dpar` alone missed
  them and they kept `brms`'s own defaults. Those defaults are actively wrong
  here: `uniform(0, min_Y)` on `ndt` reimposes the very min-RT bound the
  parameterization exists to remove, `gamma(0.01, 0.01)` on `shape` has support
  only on the positives and would silently truncate away the inverse-Weibull
  half of the family, and `poutlier` was left flat over `[0, 1]` with half its
  mass above 0.5. All three are now replaced, with priors on the natural scale:
  `lognormal(-1.2, 0.2)`, `normal(0, 0.5)` and `exponential(100)`.

  The `poutlier` prior for the omitted case has its mode at **zero**, unlike its
  modelled counterpart: leaving it out of the formula is taken to mean the data
  were trimmed or no outliers are expected. Its median is unchanged.

## Breaking changes

* **Every family is renamed to a single `cogmod_*` scheme.** `rt_lognormal()` is
  `cogmod_lognormal()`, `lnr()` is `cogmod_lnr()`, `ddm()` is `cogmod_ddm()`, and
  so on; the two LBAs are told apart by their accumulator count, so `rt_lba()`
  (no choice) is `cogmod_lba1()` and `lba()` (two choices) is `cogmod_lba2()`.
  Every derived function follows its family: `rrt_lognormal()` is
  `rcogmod_lognormal()`, `rlnr()` is `rcogmod_lnr()`,
  `rt_lognormal_stanvars()` is `cogmod_lognormal_stanvars()`, and likewise for
  the densities, the `*_lpdf_expose()` and the `brms` post-processing hooks.

  Two things were wrong with the old names. `lba()`, `lnr()`, `ddm()`, `rdm()`
  and `choco()` are generic enough to collide with anything else attached, and
  the `rt_` prefix said "reaction time" on families that are not all RT-only.
  One prefix fixes both, and `cogmod_` tab-completes the whole package.

  **The old names all still work**, as exact synonyms rather than wrappers - the
  `brms` hooks included, so a model fitted before the rename can still be
  summarised. They are undocumented and will be removed in a future release; see
  `?"cogmod-deprecated"` for the full table.

* **The `bs` dpar is renamed `boundary`** in `cogmod_invgaussian()`,
  `cogmod_lba1()`, `cogmod_lba2()`, `cogmod_rdm()` and `cogmod_ddm()`. This is not
  cosmetic: `brms` names the unpenalized spline coefficients of a model `bs`, so
  any of those families combined with a smooth failed to compile with
  `Identifier "bs" is already in use`. There is no alias for a dpar name - update
  `bs ~ ...` to `boundary ~ ...` in formulas, `bs = ` to `boundary = ` in the
  `r*()`/`d*()` calls, and expect `boundary_Intercept` where you had
  `bs_Intercept` in the output.

  Fits made before this change cannot be post-processed, because their draws
  carry the old dpar name; refit them.

* **`cogmod_lnr(link_nuzero = )` is now `link_mu`.** `brms` requires the first
  distributional parameter of a custom family to be called `mu`, so `mu` is what
  the formula uses and what the argument now matches; `nuzero` remains the name
  of the quantity, in the prose and in `rcogmod_lnr()`/`dcogmod_lnr()`.

* **The `r*()` and `d*()` functions of every shifted family now validate their
  decision parameters**, against the same bounds the generated Stan code checks.
  Only `cogmod_lognormal()` and `cogmod_loggamma()` did this before; the other
  seven silently returned zero density or `NaN` draws for, say, a negative
  `sigma`. `r*()` errors, `d*()` warns and returns zero - so a single bad
  posterior draw still cannot abort a whole call.

* **`cogmod_invgaussian()`, `cogmod_gamma()`, `cogmod_invgamma()`, `cogmod_weibull()`,
  `cogmod_invweibull()` and `cogmod_logweibull()` now use the same parameterization as
  `cogmod_lognormal()`**: `tau` and the `minrt` *dpar* are gone, replaced by `ndt`
  estimated directly (log link) plus `poutlier`, with `minrt` carried on the
  family as a constant. Every one of them gains an outlier component, and
  `with_outliers()`, `without_outliers()`, `p_outlier()` and `cogmod_priors()`
  now work on all eight families.

  Update models as for `cogmod_lognormal()`: replace `tau ~ ...` with `ndt ~ ...`,
  drop `minrt = min(df$RT)`, and add `poutlier ~ 1`. `ndt` coefficients are on
  the log scale.

  The `d*()`/`r*()` functions gain `poutlier` and `minrt` arguments, and
  `rcogmod_gamma()`, `dcogmod_gamma()`, `rcogmod_invgamma()`, `dcogmod_invgamma()`,
  `rcogmod_weibull()`, `dcogmod_weibull()`, `rcogmod_invweibull()`, `dcogmod_invweibull()`,
  `rcogmod_logweibull()` and `dcogmod_logweibull()` are new - those families previously
  had no R-level density or RNG at all. `pcogmod_invgaussian()` now returns the
  mixture CDF.

  Note that `cogmod_gamma()` and `cogmod_weibull()` have an unbounded density at `ndt`
  whenever their shape falls below 1, which the outlier component cannot repair;
  fit them with `init = 0`. `cogmod_loggamma()` nests both and lets the data choose
  the shape instead.

* **`cogmod_lba1()`** moves to the same parameterization: `tau` and the `minrt` dpar
  are replaced by `ndt` (log link) plus `poutlier`, with `minrt` a family
  constant. `rcogmod_lba1()` and `dcogmod_lba1()` gain `poutlier` and `minrt`.
  `posterior_epred_cogmod_lba1()` now explains *why* there is no expectation rather
  than calling it prohibitive: the decision time is `(b - U(0, A)) / drift` with
  a drift truncated at zero, whose density is positive at 0, so `E[1 / drift]`
  diverges and the mean does not exist.

* All nine shifted families are now generated from a single internal registry,
  so the Stan code, the R density, the RNG, the likelihood and the predictions
  cannot drift apart. Verified family by family: Stan agrees with R to machine
  precision, each density integrates to 1, and each RNG reproduces its own
  density.

* **`cogmod_lnr()` moves to the same `ndt` + `poutlier` parameterization as the
  RT-only families.** `tau` and the `minrt` *dpar* are gone; `ndt` is estimated
  directly, on the log link, and `minrt` is a constant carried on the family
  object rather than something `tau` is scaled by. `with_outliers()`,
  `without_outliers()`, `p_outlier()`, `cogmod_priors()` and `cogmod_inits()`
  all work on it now.

  Because the LNR produces a **choice as well as a time**, its outlier
  component is not quite the RT-only one: the contaminant guesses uniformly
  over the two response options in addition to drawing an RT from the same
  half Student-t, so the mixture is `poutlier * (1/2) * g(t) + (1 - poutlier)
  * f_k(t - ndt)`. The `1/2` is what keeps the joint density summing to one
  over both responses - without it the total comes to `1 + poutlier`.

  Update models as for `cogmod_lognormal()`: replace `tau ~ ...` with
  `ndt ~ ...`, drop `minrt = min(df$RT)` from `bf()`, and add `poutlier ~ 1`.
  Fit with `init = cogmod_inits(f, df)` rather than `init = 0` - on the log
  link, `init = 0` starts `ndt` at `exp(0) = 1` second, above nearly every
  sub-second RT, which leaves every response attributed to the outlier
  component and the race parameters with no gradient at all.

  `rcogmod_lnr()` and `dcogmod_lnr()` gain `poutlier` and `minrt` arguments.
  Verified against the RT-only families' checklist: the joint density sums to
  one over both responses and integrates to one over time (with and without
  the `1/K` term, to confirm it is load-bearing), the Stan `cogmod_lnr_lpdf`
  agrees with the R density to machine precision, and a simulated fit recovers
  `ndt` well above the fastest observed response.

* **`cogmod_rdm()` moves to the same `ndt` + `poutlier` parameterization**, on
  the same shared machinery as `cogmod_lnr()`. `tau` and the `minrt` *dpar* are
  gone: the family's dpars are now `mu`, `driftone`, `sigmabias`, `boundary`,
  `ndt`, `poutlier`, `ndt` is estimated directly on the log link, and `minrt` is
  a constant carried on the family object. `with_outliers()`,
  `without_outliers()`, `p_outlier()`, `cogmod_priors()`, `cogmod_inits()` and
  `cogmod_stanvars()` all work on it now.

  Update models as for `cogmod_lnr()`: replace `tau ~ ...` with `ndt ~ ...`,
  drop `minrt = min(df$RT)` from `bf()`, and add `poutlier ~ 1`. Fit with
  `init = cogmod_inits(f, df)` rather than `init = 0` or `init = 0.5`.

  `cogmod_priors()` also supplies the `sigmabias` / `boundary` priors that
  `?cogmod_rdm` previously told you to write by hand: the two enter the model
  only through the sum `b = boundary + sigmabias` and trade off along a ridge
  worth a handful of log units, which under a flat prior and a `softplus` link
  is an improper posterior. The drift rates are left flat on purpose - unlike
  `cogmod_lnr()`'s `nuone`, a drift pushed to zero does not produce a plateau,
  because a driftless accumulator still finishes and still wins sometimes.

  `rcogmod_rdm()` and `dcogmod_rdm()` gain `poutlier` and `minrt` arguments,
  and `pcogmod_rdm()` gains `poutlier` and `minrt` too - its CDF is now the
  mixture's, and still keeps the far-tail survival in log space. `dcogmod_rdm()`
  keeps its `response = NULL` marginal, which is now the sum of the two
  defective mixture densities. Invalid parameters now warn and return a zero
  density rather than erroring, matching the other mixture families. **A drift
  of exactly zero is still accepted** - it is the one closed lower bound in
  either registry, because driftless Brownian motion still reaches the
  threshold with probability one.

  Verified against the checklist: the joint density sums to one over both
  responses and integrates to one over time for several parameter sets and
  `poutlier` values (including as the start-point range shrinks to `1e-6`), the
  Stan `cogmod_rdm_lpdf` agrees with the R density across the parameter grid,
  and a simulated fit recovers `ndt = 0.256` against a true `0.25` - some 200
  times the fastest observed response, which the old `tau * minrt` bound could
  not have expressed.

* **`cogmod_lba2()` moves to the same `ndt` + `poutlier` parameterization.**
  Its dpars are now `mu`, `driftone`, `sigmazero`, `sigmaone`, `sigmabias`,
  `boundary`, `ndt`, `poutlier`. Update models as for `cogmod_lnr()`: replace
  `tau ~ ...` with `ndt ~ ...`, drop `minrt = min(df$RT)` from `bf()`, and add
  `poutlier ~ 1`.

  Two long-standing bugs in the density came out with it.

  **The density was not normalised.** A normal drift rate can come out negative,
  and such an accumulator never reaches the threshold, so a trial on which
  *both* drifts are negative produces no response at all. `rcogmod_lba2()` has
  always resampled until at least one is positive - but `dcogmod_lba2()` never
  divided by the probability of that event, so the density integrated to the
  probability rather than to one. At drift rates of `0.5` and `0.2` with SDs of
  `1.5` it came to `0.83`, and the simulated choice proportion was `0.562`
  against an integral of `0.468`. Because the shortfall depends on the
  parameters, it biased estimates rather than merely offsetting the likelihood.
  Both the R and the Stan densities now condition on the event the process is
  conditioned on.

  **The `(1 / A)` cancellation `cogmod_lba1()` was fixed for was still here.**
  The defective density divides `drift * (Phi(z2) - Phi(z1)) + sigma * (phi(z1)
  - phi(z2))` by the start-point range, and both differences vanish linearly in
  it; the loser's survival was computed as `1 - CDF`, which cancels the same way.
  Relative error reached 2.7% at `sigmabias = 1e-5` and 320% at `1e-7`. Both now
  go through the kernels `cogmod_lba1()` already uses (`.lba_dens_over_A()`, and
  a new `.lba_surv_raw()` that takes the survival directly rather than as
  `1 - CDF`), which the two families now share in R and in Stan.

  The `.Machine$double.eps` floor is gone too. It turned every RT below the
  point where the density becomes representable into a log-density of exactly
  `-36.04` - a constant the model never produced, with a gradient of zero. Where
  the density really has underflowed the log-density is now `-Inf`, and the
  outlier component is what keeps the *mixture* finite there.

  `rcogmod_lba2()` now imposes the positive-drift condition **exactly**, by
  sampling which accumulator is positive and then the truncated normals, instead
  of a rejection loop. Its `max_iter` argument is therefore gone, along with the
  fallback that forced a drift positive with `abs()` - and so drew from the
  wrong distribution - whenever the loop ran out.

  Note that the evidence scale of an LBA is **arbitrary**: multiply the drifts,
  their SDs, the start-point range and the threshold by any `c > 0` and every
  finishing time is unchanged, so the likelihood is exactly constant along that
  ray. `cogmod_priors()` now fences all four positive parameters, which makes
  the posterior proper, but only fixing one SD in the formula
  (`sigmazero = 1` in `bf()`) identifies the scale. This was true before and is
  now documented.

* **`cogmod_ddm()` moves to the same `ndt` + `poutlier` parameterization**, and
  **`sigmatau` is renamed `sigmandt`**. Its dpars are now `mu`, `boundary`,
  `bias`, `sigmadrift`, `sigmabias`, `sigmandt`, `ndt`, `poutlier`.

  `sigmatau` was the between-trial range of the non-decision time expressed as a
  fraction of `minrt` (`st0 = sigmatau * minrt`). With `tau` and the `minrt`
  dpar both gone it was named after a parameter that no longer exists and scaled
  by a constant the user no longer sees, so it is now **`sigmandt`**, which is
  `st0` itself, in the same unit as the data, on a log link - `ndt` remains the
  lower bound of the resulting Uniform. Update models by replacing
  `tau ~ ...` with `ndt ~ ...`, dropping `minrt = min(df$RT)`, adding
  `poutlier ~ 1`, and rewriting any `sigmatau` term as `sigmandt` in seconds.

  All three between-trial variability parameters remain legitimately **zero**,
  and fixing them in the formula (`sigmadrift = 0`) still recovers the classic
  4-parameter DDM. `cogmod_priors()` now supplies priors for all three:
  each has a floor at zero that its link reaches only at minus infinity, and the
  likelihood stops changing well before then, which under `brms`'s flat default
  is an improper posterior.

  `posterior_epred_cogmod_ddm()` keeps its closed form - the DDM is the one
  choice family here with a usable one - now using `ndt` directly and blending
  in the outlier component when `predict_outliers` is set, like the RT-only
  families. `rcogmod_ddm()` and `dcogmod_ddm()` no longer take `...` for
  `brms::rwiener()`/`brms::dwiener()`, which the shared mixture machinery cannot
  forward; set `options(wiener_backend = )` instead.

  Both R and Stan evaluate the decision component at a non-decision time of
  zero, which `dwiener()`, `rwiener()` and `wiener_lpdf()` all refuse. Since the
  Wiener density depends on the time and the non-decision time only through
  their difference, both offset the pair by the same `1e-10`, and the Stan
  literal is generated from the R constant so the two cannot drift apart.

  All four choice+RT families are now on the direct `ndt` + `poutlier`
  parameterization, and `tau` + `minrt` is gone from the package.

## Bug fixes

* **`dcogmod_lba1()` was wrong for small start-point ranges.** Its density is built
  from `drift * (Phi(z2) - Phi(z1)) + sigma * (phi(z1) - phi(z2))` divided by the
  start-point range `A`, and both differences vanish linearly in `A` - so
  evaluating them directly and then dividing lost every significant digit once
  `A` was small. The old code also floored the bracket at `1e-10`, which turned
  that underflow into a spurious density floor spread over the whole line. The
  result: the density stopped integrating to one below about `A = 0.1` and was
  outright divergent below `A = 0.01`.

  Both differences are now computed stably - a Taylor expansion in
  `delta = A / (sigma * t)` below `1e-4`, and tail-aware differencing above it -
  and the floor is gone. The density integrates to 1 from `A = 2` down to
  `A = 1e-8`, and converges to the recinormal (LATER) limit at the expected
  first-order rate. The same fix is in the Stan code, which matches the R
  density across 960 parameter combinations.

## Bug fixes (parameterization)

* `posterior_epred_cogmod_logweibull()` returned `exp(mu + sigma * 0.5772)`, which is
  the *geometric* mean of the decision time - the exponential of `E[log(RT)]` -
  rather than its mean. It now returns `exp(mu) * gamma(1 - sigma)`, and `Inf`
  where `sigma >= 1` and no mean exists.

* `posterior_epred_cogmod_invweibull()` now returns `Inf` where the Frechet shape is
  `<= 1` and the mean does not exist, instead of a finite but meaningless value.

* `cogmod_lognormal()` no longer uses `tau` and `minrt`. Non-decision time is now
  estimated directly as `ndt`, in seconds, through a log link, and the family
  gains `poutlier`, the proportion of trials generated by an outlier process
  rather than by the decision process.

  The old parameterization set `ndt = tau * minrt` with `tau` in `(0, 1)` and
  `minrt` injected as a constant, which capped `ndt` at an order statistic of
  the sample. Because that cap was shared across the whole dataset, a
  non-decision time larger than the fastest observed response was
  inexpressible - so any condition or participant whose true `ndt` exceeded the
  global minimum RT could not be recovered, and the misfit surfaced instead as
  spurious effects on the other parameters.

  What makes the direct parameterization tractable is the outlier component: a
  fixed half Student-t (scale `0.4`, `3` df) mixed in with weight `poutlier`
  keeps the density positive below `ndt`, turning the hard min-RT boundary into a
  finite cost and leaving the log-density smooth. No bound is taken from the
  data. The half-t is flat at the origin, so the fastest responses are not
  starved of density; its tails are heavy enough that supplying RT in
  milliseconds degrades rather than underflowing to zero; and its mean is finite,
  which `posterior_epred()` requires.

  Models must be updated: replace `tau ~ ...` with `ndt ~ ...`, drop
  `minrt = min(df$RT)`, and add `poutlier ~ 1`. Note that `ndt` coefficients
  are on the log scale, so `exp()` them for seconds.

* The outlier component is scaled by `minrt`, the fastest reaction time that
  could plausibly be a real decision, used directly as the half-t scale with no
  conversion factor. It defaults to `0.3` seconds, where the conditional
  accuracy functions in the *Outliers* article show responses sitting at chance
  across three paradigms. 61% of the component falls below `minrt` whatever
  value is chosen.

  It is a judgement about the task rather than a statistic: nothing is read off
  the sample, and `ndt` is not bounded by it.

  This matters because the shifted LogNormal is scale-equivariant on its own -
  multiply every response by 1000 and `ndt` comes back multiplied by 1000 - and
  a component pinned to the second would be the one thing breaking that. With
  millisecond data and a fixed scale, the outlier component sits many orders of
  magnitude below the decision density everywhere in the data, `poutlier`
  collapses toward zero and `ndt` reverts to being pinned by the fastest
  observed response, silently and without a warning. Setting `minrt` in the unit
  of the data makes the likelihood exactly equivariant instead, so millisecond
  data need `minrt = 300`.

  `minrt` is a constant, never estimated, and is deliberately **not** a dpar:
  `brms` has no notion of a default for one, so a dpar omitted from the formula
  is estimated rather than defaulted. It is carried on the family, and the
  matching Stan constant comes from `cogmod_lognormal_stanvars()`, which also
  accepts the family itself so the two cannot drift apart:

  ```r
  fam <- cogmod_lognormal(minrt = 0.3)
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  brms::brm(f, data = df, prior = cogmod_priors(f, df),
            stanvars = cogmod_lognormal_stanvars(fam))
  ```

  `rcogmod_lognormal()` and `dcogmod_lognormal()` gain a matching `minrt` argument.
  Models fitted before this change fall back to the default, so their
  predictions are unaffected.

* The `wagenmakers2008` dataset has been removed. Those data were supplied by
  the original authors for distribution in the `rtdists` package specifically,
  and no open licence covers them, so redistributing them here was not
  appropriate. They remain available as `rtdists::speed_acc`; `rtdists` is in
  `Suggests`, and the vignettes and paper now reconstruct the same subset from
  it (`!censor & response != "error" & rt <= 2`), so previously reported
  results are unchanged.

* The experimental confidence signal detection model (`rconf_sdt()`,
  `dconf_sdt()`, `conf_sdt_stanvars()`, `conf_sdt_custom_family()` and its
  `brms` methods) is no longer exported. It was not ready, and the code is
  commented out in `R/conf_sdt.R` and `R/conf_sdt_brms.R` pending a rework.

* **Priors are now required.** `brms` assigns a flat, improper prior to the
  intercept of any custom-family parameter it does not recognise, which here
  means both `ndt` and `poutlier`, and the likelihood has two flat directions
  that a flat prior turns into an improper posterior: `poutlier` toward 1, where
  every response is attributed to the outlier component and `mu`, `sigma` and
  `ndt` drop out of the density altogether; and `ndt` toward 0, where the model
  reduces to an unshifted LogNormal and the gradient with respect to `log(ndt)`
  vanishes. The second is inherent to putting a positive shift on a log link and
  has nothing to do with the mixture, which is why a prior on `poutlier` alone
  is not enough. Symptom: intercepts around `1e14`, `Rhat` near 2 and an
  effective sample size of about 5, with no error raised. `cogmod_priors()`
  fills the gap.

## New features

* `cogmod_priors(formula, data)` fills in every prior `brms` would otherwise
  leave flat - for `cogmod_lognormal()`, the `ndt` and `poutlier` rows. It starts
  from `brms::get_prior()` for the model in hand rather than guessing, so
  `0 + Intercept` formulas, interactions, group-level terms and smooths are all
  handled and a prior matching no parameter is impossible by construction. The
  result is passed through `brms::validate_prior()`, so a malformed
  specification errors there with the offending row in view, and the return
  value is the complete prior table with a `source` column marking each row as
  `user` or `default` - print it to see exactly what the model will be fitted
  with. To change one, edit the row; `c()` will not work for a slot the table
  already covers, because `brms` rejects two priors for the same slot.

  The family is read off the formula, so build it with
  `brms::bf(..., family = cogmod_lognormal())`. Any other family, or a formula
  carrying none, gets a message and the `brms` defaults unchanged, so the call
  is always safe to leave in a script.

* `p_outlier()` returns the posterior probability that each trial came from the
  outlier component rather than the decision process - the mixture
  responsibility, averaged over draws. Responses below `ndt` come out at 1 and
  those in the bulk near 0, but the probability rises again in the far slow tail,
  where the half-t has heavier tails than the LogNormal; that is the mechanism
  behind the advice to filter implausibly slow responses before fitting. The
  responsibility is computed on the log scale, so it stays finite in tails where
  both components underflow to zero. It returns `rt` and `p_outlier` only; the
  `fast` column was a marginal median split that ignored any grouping in the
  model and is gone.

* `posterior_predict()` and `posterior_epred()` for `cogmod_lognormal()` describe
  the **decision process alone** by default, as if `poutlier` were zero. For
  visualising effects the outlier component is a nuisance that pulls expected
  values toward its own mean (0.441 s) and adds a spike of implausibly fast
  draws; it is also a fixed regularizer rather than a claim about how guesses
  are distributed, so simulating from it means simulating from something the
  model does not assert. The likelihood is unaffected and is always the full
  mixture, so `posterior_predict()` and `log_lik()` no longer describe the same
  distribution - a hand-rolled LOO-PIT check should use `with_outliers()`.

* `with_outliers()` restores the fitted mixture for prediction, and
  `without_outliers()` returns to the default. The main use for the former is
  `pp_check()`: on untrimmed data the decision-only predictive has no fast spike
  to match the one in the data, which reads as misfit. The same flag can be set
  up front with `cogmod_lognormal(predict_outliers = TRUE)`.

  It is carried on the family rather than passed as an argument because `brms`
  does not forward extra arguments to a custom family's prediction methods -
  `posterior_epred` reaches the method with `prep` and nothing else - and
  `insight`, `modelbased` and `marginaleffects` inherit that. Carrying it on the
  object is what makes it work through all of them.

* New *Outliers* article validating the outlier-mixture specification against
  the Illusion Game dataset, cross-checked against the lexical decision data of
  Wagenmakers et al. (2008) and the brightness discrimination of Ratcliff and
  Rouder (1998), both from `rtdists`. It also derives recommended priors for all
  four parameters from the rates it measures, and explains why they matter more
  under this parameterization than the last: `ndt` is no longer bounded above by
  construction, and `poutlier` has a degenerate region near 1 where `ndt` becomes
  unidentified.

# cogmod 0.1.0

## Models for subjective scales

* Beta-Gate (`cogmod_betagate()`, `rcogmod_betagate()`, `dcogmod_betagate()`), a reparametrised
  ordered beta model.
* Discrete Beta (`cogmod_betadiscrete()` and friends) for Likert-type responses.
* Choice-Confidence, CHOCO (`cogmod_choco()`, `rcogmod_choco()`, `dcogmod_choco()`) for bipolar
  scales.
* Signal detection with confidence ratings (`conf_sdt()`).

## Models for decision making

* Lognormal race (`cogmod_lnr()`), linear ballistic accumulator (`cogmod_lba2()`), drift
  diffusion with optional across-trial variability (`cogmod_ddm()`), and the racing
  diffusion model (`cogmod_rdm()`).

## Models for reaction times alone

* A consistently parametrised set of shifted, right-skewed response
  distributions: `cogmod_lognormal()`, `cogmod_invgaussian()`, `cogmod_gamma()`,
  `cogmod_invgamma()`, `cogmod_weibull()`, `cogmod_logweibull()`, `cogmod_invweibull()`,
  `cogmod_exgaussian()` and `cogmod_lba1()`.

## Data

* `wagenmakers2008`, lexical decision data from Wagenmakers et al. (2008).
* `badlm`, a simulated dataset in which two conditions share a mean reaction
  time while differing in shift, spread and tail weight.
