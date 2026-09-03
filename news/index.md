# Changelog

## cogmod 0.3.1

### New features

- **[`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  gains `sigmandt`**, the between-trial range of the non-decision time
  (`st0`): each trial’s non-decision time is drawn from
  `Uniform(ndt, ndt + sigmandt)`, so `ndt` becomes its lower bound,
  exactly as in
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md).
  Spreading the shift turns the Wald density into a difference of two
  CDFs and the CDF into a difference of two integrated CDFs, both closed
  form at a fixed drift, so the parameter costs a few normal CDFs per
  observation, works with `cens()` unchanged, and rides the existing
  drift quadrature when `sigmadrift > 0` too. It is on a `log` link with
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)’s
  prior for the same quantity. **It is hard to estimate and should be
  fixed at zero for most applications** (`sigmandt = 0` in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html));
  it shares the leading edge of the distribution with `ndt` and
  `poutlier`, and should only be freed with a lot of data, a strong
  prior, or both. As with `sigmadrift`, leaving it out of
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  *estimates* it, so existing Wald formulas that do not mention it now
  fit one more parameter unless they add `sigmandt = 0`, and fits made
  before this version cannot be post-processed with it; the vignette
  models were refit.
  [`rcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  [`dcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  and
  [`pcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  take `sigmandt` right after `sigmadrift`, so a `poutlier` passed by
  position moves along one.

- **Censored reaction times: `brms`’s `cens()` works on the RT-only
  families.** `bf(rt | cens(error) ~ ...)` scores an error trial as a
  *right-censored correct response*: its RT is a lower bound on when the
  correct process would have finished, so it contributes that process’s
  survival rather than its density. On
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  this is the censored shifted Wald of [Miller et
  al. (2018)](https://doi.org/10.1177/0146621617710465), the `cswald`
  model of [`bmm`](https://github.com/popov-lab/bmm). Here it is not a
  family but a construction, so the same formula works on every RT-only
  family with a closed-form CDF:
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md),
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`cogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
  [`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md),
  [`cogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md),
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  and
  [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md).
  Left- and interval-censoring come with it, and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html) -
  hence [`loo()`](https://mc-stan.org/loo/reference/loo.html) - honours
  all three, which `brms` leaves to a custom family’s own method.
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  predicts the latent, uncensored RT, as `brms` does for its own
  families.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  and
  [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)
  refuse `cens()` on the families that cannot take it, and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  warns above 20% censored trials. What the model is for, what it
  assumes and the one check to run before using it are in
  [`?rcogmod_invgaussian`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  and the *Censored Shifted Wald* section of `vignette("rt_models")`.

  Under the hood every censorable family gets a `<family>_lcdf` and a
  `<family>_lccdf` beside its `_lpdf`, generated from two new registry
  slots so a family cannot drift out of step with itself. The survivals
  are written as survivals - never as `log(1 - exp(lcdf))` - and the
  half Normal outlier’s through `std_normal_lcdf(-z)` rather than
  `std_normal_lccdf(z)`, which is `-inf` from 1.66 s on: the two places
  `bmm`’s implementation broke. With `sigmadrift > 0` the Wald CDF has
  no closed form and is taken by 64-point Gauss-Legendre quadrature over
  the drift, in R and Stan alike off one node table.

- **`pcogmod_*()` for every censorable family.**
  [`pcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  [`pcogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md),
  [`pcogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`pcogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`pcogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`pcogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
  [`pcogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md),
  [`pcogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md)
  and
  [`pcogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  join
  [`pcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  with `lower.tail` and `log.p`. The upper tail is computed *as* the
  upper tail rather than as `1 - CDF`; these are the R side of the Stan
  `_lcdf`/`_lccdf` pair, and the tests hold the two to each other.

- **[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now checks the response before returning.** Everything in this package
  is stated in seconds and none of it is unit-equivariant - the `ndt`
  prior means 170-300 ms, `.POUTLIER_SCALE` is 0.2 s - but `brms` fills
  its own defaults from the data, so a column of milliseconds produces a
  model whose two halves silently describe different quantities. It
  compiles, it samples, it converges, and the estimates are meaningless.
  The check catches that and the handful of other mistakes with the same
  character.

  It **stops** where the offending rows would make the fit impossible or
  wrong in a way `Stan` cannot report: a non-positive reaction time
  under a family that places no density below `ndt`; a response outside
  `[0, 1]` for
  [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md)
  or
  [`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md);
  a non-integer rating for
  [`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md);
  a non-numeric response; and a third level in `dec()`, which the choice
  families would otherwise fold silently into option 1, since their Stan
  code tests `dec == 0` and takes the else branch for everything else.

  It **warns** about the rest: a median implying milliseconds, `NA`s,
  and either tail running past what `poutlier` can absorb. The tails are
  judged as proportions rather than counts, because the outlier
  component is *supposed* to produce the occasional fast response -
  `rcogmod_lognormal(200, ndt = 0.2, poutlier = 0.02)` puts one at 81
  ms - and a count-based test fires on the package’s own generator. Over
  20000 draws the component sends 0.8% of responses below 0.1 s at
  `poutlier = 0.02` and 1.9% at 0.05, the top of the default prior, so
  the warning sits at 5%.

  Families with neither `ndt` nor `poutlier` -
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  and
  [`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md) -
  are exempted from the tail checks, and a non-positive response is a
  warning rather than an error there, their support being the whole real
  line. A formula or family the check cannot read is passed through
  untouched, so `brms`’s own error is what the user sees.

### Performance

- **The R-side DDM density no longer goes through
  [`brms::dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html).**
  [`dcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md),
  and with it
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md)
  and every other post-processing method that evaluates the likelihood
  in R, now use a vectorised Navarro and Fuss (2009) series written in
  log space. The 4-parameter density is about eight times cheaper per
  element and agrees with
  [`brms::dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
  to `1e-12` on the log scale. The 7-parameter density, which evaluates
  that series 625 times per observation under Gauss-Legendre quadrature,
  goes from about 5 ms to about 0.6 ms per draw-observation - a LOO over
  4000 draws of 500 trials drops from close to three hours to about
  twenty minutes. Both now return a finite log-density in the far tails
  where
  [`dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
  returns `log(0)`. The Stan likelihood is unchanged. `RWiener` is still
  needed by the test suite, which uses
  [`dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html) as
  the reference.

### Breaking changes

- **[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  now truncates each drift rate at zero**, the convention of `rtdists`
  (`posdrift = TRUE`), `DMC`, `EMC2` and `ggdmc`. Previously the pair of
  drifts was conditioned on at least one being positive and a losing
  accumulator was allowed a negative rate, which it kept forever. The
  two are different models of the same race wherever a drift is small
  relative to its SD: densities up to about 40-50% apart in the tails,
  choice probabilities a few percentage points apart. The change makes
  `cogmod` LBA estimates directly comparable with those packages and
  with the literature built on them, and
  [`dcogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  now reproduces
  [`rtdists::dLBA()`](https://rdrr.io/pkg/rtdists/man/LBA.html) at the
  same parameter values. It also simplifies the sampler, which draws
  each drift from its truncated Normal rather than splitting the
  conditional law into cases. **Fits made with earlier versions cannot
  be post-processed with this one**, and their estimates are not
  comparable with new ones at low drift rates; the vignette model was
  refit. The loser’s survival is computed as
  `P(v > 0, unfinished) / P(v > 0)` from whichever tail keeps its
  digits, so the density stays accurate for a loser with a strongly
  negative drift, where both quantities are tiny, and in the far tail,
  where `1 - CDF` would cancel.
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  is unaffected: with one accumulator the two conventions coincide.

  The truncation has one cost, and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now covers it: once an accumulator rarely wins, its `drift` and
  `sigma` are identified only through `|drift| / sigma^2` (the truncated
  Normal converges to an Exponential along that ray), so a flat prior
  lets the drift run off - the vignette’s error accumulator sat at `-12`
  with an interval of `-23` to `-6.5`.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  therefore puts `normal(1, 2)` on `driftone` and `normal(0, 1.5)` on
  its slopes, the treatment
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)’s
  `nuone` already had. Existing formulas that left `driftone` to `brms`
  get this prior on their next
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  call.

## cogmod 0.3.0

### New features

- New family
  **[`cogmod_geg()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_geg.md)**,
  the Generalised Ex-Gaussian of [Marmolejo-Ramos et
  al. (2023)](https://doi.org/10.1007/s11571-022-09813-2): the
  ex-Gaussian with its CDF raised to a power,
  `F_GEG(x) = F_EG(x)^shape`. The construction is Durrans’ alpha-power
  family, so the density is `shape * F_EG(x)^(shape - 1) * f_EG(x)` and
  stays closed form - Stan ships both `exp_mod_normal_lpdf` and
  `exp_mod_normal_lcdf`, so the whole likelihood is three lines and
  costs one extra CDF.

  `shape = 1` is
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  **exactly**, not approximately, in R and in Stan alike, so
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  between the two is like-for-like.

  What it buys is shape. Sweeping `sigma` and `tau` across the values RT
  data occupy, the ex-Gaussian spans skewness 0 to 2 and excess kurtosis
  0 to 6; freeing `shape` widens that to roughly -0.4 to 4.8 and 0
  to 35. In particular the GEG can be **negatively skewed**, which the
  ex-Gaussian cannot be at any parameter value.

  What it costs is interpretability, and specifically the property the
  ex-Gaussian is normally reported for. **The mean is no longer
  `mu + tau`** - at `mu = 0.4`, `tau = 0.2` it runs 0.31 at
  `shape = 0.2` and 1.15 at `shape = 20` - and it has no closed form, so
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  integrates numerically and is the one generic materially slower here
  than for a closed-form family. On a full data set, summarise
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  draws instead.

  `shape` is also badly confounded with `mu`: fitted by maximum
  likelihood to the lexical-decision data used in the vignettes, the two
  correlate about -0.98 at the optimum, and the other estimates move
  with it (on one condition `mu` goes 0.429 to 0.508, `sigma` 0.051 to
  0.037, `tau` 0.119 to 0.162). `shape` re-slices the same bulk-and-tail
  split rather than adding an independent axis.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  therefore gives it a deliberately informative `normal(0, 0.5)` on the
  `log` link - centred on `shape = 1`, the ex-Gaussian - and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  starts it there.

  Use it when fit is the point. When the `mu`/`tau` decomposition is the
  point,
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  is the family to fit; when a better-fitting descriptive family is the
  point,
  [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md)
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  decouple skew from tail weight with parameters that stay
  interpretable.

- New vignette, `vignette("performance")`, covering how to speed up
  sampling for the choice+RT families: approximating a fit first with
  Pathfinder before committing to full MCMC, chain and within-chain
  (multithreading) parallelization, spreading short warmup runs across
  an HPC job array and recombining them with
  [`brms::combine_models()`](https://paulbuerkner.com/brms/reference/combine_models.html),
  and amortized inference (e.g. BayesFlow) as a longer-term direction.

- **[`pcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  gains a `response` argument**, giving the *defective* CDF
  `P(RT <= q, choice = response)` rather than only the RT distribution
  marginally over the choice. It does not reach one - its limit is the
  probability of that response, which `pcogmod_rdm(Inf, response = k)`
  gives - and the two of them sum back to the marginal CDF. This is what
  a defective-CDF or quantile-probability plot of a race model needs.

  `response` goes *after* `poutlier` in the signature, unlike in
  [`dcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`pcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md),
  which put it before. `poutlier` is a model parameter and `response` is
  not, so every parameter now comes first and the argument stays where a
  positional call already expects it - only a positional `lower.tail` or
  `log.p`, in eighth place or later, is affected.

  There is no closed form for it, unlike the marginal CDF, where the
  race survival factorises as `S0 * S1`. It is quadrature over the
  defective density, so it is accurate to about `1e-8` rather than to
  machine precision, and about ten times slower per element - 200 points
  take 0.6 s against 0.07 s for the marginal. `lower.tail = FALSE`
  integrates the upper side directly rather than subtracting, so the
  defective survival stays accurate into the tail.

- **New
  [`qcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)**,
  the quantile function, inverting
  [`pcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  by root-finding - marginally, or per response. `scale_p = TRUE` reads
  `p` as a fraction of the chosen response’s own probability, so
  `p = 0.5` is that response’s median; that is the form a
  quantile-probability plot wants.

- **[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  now accepts a start-point range of exactly zero.** `bias = 0`
  (`sigmabias` in the `brms` family) is a model, not a degenerate
  parameter: both accumulators start at `0` on every trial and the race
  is between two plain Walds - equation 2 of Tillman et al. (2020),
  which is the limit the density already took.
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  and
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  have always allowed it and the RDM excluding it was an inconsistency.
  Accepted in R and in Stan alike; a negative range is still rejected.

  This does not make the `sigmabias` direction any better identified.
  The `softplus` link still reaches zero only at minus infinity, so the
  flat `sigmabias -> 0` ridge is as long as it ever was and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  still fences it.

### Bug fixes

- **A missing, infinite or zero-length reaction time no longer aborts a
  density.** Of the sixteen mixture families, four threw
  `missing value where TRUE/FALSE needed` on `dcogmod_*(NA_real_)` -
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md),
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md),
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) -
  and two threw on `Inf`:
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  and
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md).
  The rest returned `0`. The split was not principled: those are the
  densities whose cores branch on a comparison rather than being
  arithmetic all the way down, and `any(NA)` is `NA` while `if (NA)` is
  an error. All sixteen now return `0`, which is what the shared
  machinery already wrote for those rows and could never reach. A
  missing **parameter**, and a missing **response**, are handled the
  same way and for the same reason.

  This matters beyond tidiness: one bad entry used to take the whole
  vector down with it, so a single `NA` in a column of reaction times
  aborted the call instead of costing that row.

  [`dcogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  still returns `NA` for `NA`. It is a plain density with no outlier
  component, so it follows
  [`dnorm()`](https://rdrr.io/r/stats/Normal.html) rather than the
  mixture convention.

- **Zero-length input gives a zero-length answer.**
  `dcogmod_*(numeric(0))` returned a length-1 value - or, in the same
  four families, threw - because the shared preparation recycled the
  empty vector up to the parameters and `rep_len(numeric(0), 1)` is
  `NA`. It now returns `numeric(0)`, as every `d`/`p`/`q` function in
  base R does. A zero-length *parameter* alongside a real quantile is
  now rejected rather than silently becoming a vector of `NA`s.

- **`pcogmod_ddm(q, response = k, lower.tail = FALSE)` was not a
  survival.** It returned `1 - P(RT <= q, choice = k)`, which is
  `P(RT > q OR choice != k)`; at `q = Inf` that gave the probability of
  the *other* response rather than zero. It is now the defective
  survival `P(RT > q, choice = k)`, so the two tails add to that
  response’s own probability rather than to one. The marginal
  (`response = NULL`) is unchanged, and so is the lower tail in both
  forms. Same convention as
  [`pcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md).

- `pcogmod_invgaussian(NA_real_)` threw rather than returning `NA`.

- `pcogmod_rdm(q, lower.tail = TRUE)` returned `NaN` instead of `0` at
  `q <= 0` when `poutlier > 0`. Both mixture components are exactly
  `log(1)` there, and the mixture of them landed 2e-17 *above* zero
  rather than on it, which `log(1 - exp(.))` cannot take. Reached in
  practice by
  [`qcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
  whose root search starts at zero.

- `pcogmod_rdm(NA)` returned `1` rather than `NA`, the mixture helper
  mapping a missing value onto a log-survival of `-Inf`.

- The direct lower-tail branch of `.pwald()` had no small-`bias` case;
  both of its branches divide by the start-point range. Unreachable
  before, since `bias = 0` was rejected.

### Breaking changes

- **The pre-rename names are gone.** `rt_lognormal()`, `lnr()`, `ddm()`,
  `rdm()`, `choco()`, `betagate()`, `betadiscrete()`, `lba()`,
  `rt_lba()` and every function derived from them - 140 exports in all -
  were kept as synonyms through 0.2.1 and are removed here. Use the
  `cogmod_*` name: `rt_lognormal()` is
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  `rrt_lognormal()` is
  [`rcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  `rt_lognormal_stanvars()` is
  [`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  and likewise for the densities, the `*_lpdf_expose()` and the `brms`
  post-processing hooks.

  The synonyms existed so that a model fitted before the rename could
  still be summarised, since `brms` looks up `log_lik_<family>()` by the
  name stored on the fit. That window closes with the first CRAN
  release: there is no released version to be compatible with, and a fit
  made under an old name can be brought forward by setting
  `fit$family$name` to the `cogmod_*` one. Refitting is the safer route,
  as several parameterizations changed in 0.2.1 as well.

- **[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)’s
  `mu` is now on an `identity` link** and is unbounded, where it was on
  `softplus` with a lower bound of zero. `sigma` and `tau` are
  unchanged.

  `mu` is the **location** of the Gaussian component, not a scale. The
  convolution is well defined for any real value - the density
  integrates to one at `mu = 0` and below - and the Stan `lpdf` has
  always agreed, checking only `sigma` and `tau`. The old bound lived in
  `.prepare_exgaussian()` and in the family declaration alone, so the R
  functions were refusing inputs the sampler would happily fit.

  Two things were wrong with constraining it. Interpretability, which is
  most of the point of the ex-Gaussian: behind `softplus` a coefficient
  is not in seconds, and the conversion factor moves with the
  intercept - the local slope is 0.33 at `mu = 0.4` s, 0.39 at 0.5 s and
  0.63 at 1 s, so the same effect reads as a different number depending
  on where the intercept sits. And fidelity: for fast, heavily-tailed
  data the Gaussian component genuinely belongs near or below zero with
  `tau` carrying the mass, and forcing `mu > 0` distorts the `mu`/`tau`
  split in exactly the cases where that decomposition is the quantity
  being estimated. `identity` also matches every other implementation -
  `brms`’s own
  [`exgaussian()`](https://paulbuerkner.com/brms/reference/brmsfamily.html),
  `retimes`, and the estimates in the literature - so fitted values are
  now directly comparable.

  **What to change.** Coefficients on `mu` are now in seconds and are
  not comparable to values from an earlier fit; refit rather than
  reinterpret. Pass `cogmod_exgaussian(link_mu = "softplus")` to keep
  the old behaviour.

- **[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now sets a prior on
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)’s
  `mu` intercept**, `normal(0.4, 0.25)`. It previously left `mu` to
  `brms`, whose `student_t(3, 0, 2.5)` was a fair statement on the
  softplus scale (median 0.69 s) but on `identity` is centred on zero
  seconds and rates a Gaussian centre of -2 s as plausible as one of
  +2 s. The prior deliberately does not exclude negative values. Only
  the intercept is set; the response’s slopes are the effects being
  estimated and are left alone.

- The **outlier component is now a half Normal with a fixed scale of 0.2
  s**, where it was a half Student-t with 3 degrees of freedom and a
  user-supplied scale. Two things changed, for one reason.

  The Student-t’s tail was heavier than every decision density in the
  package, so far-out slow responses were eventually explained better by
  the outlier component than by the model: against a shifted LogNormal
  at `poutlier = 0.02`, a 5 s response was attributed to it with
  probability 0.86 and the crossover sat at 3.86 s, with `ndt` pulled up
  behind it. `vignette("outliers")` already flagged this as a defect. A
  Gaussian never gets there - the same responsibility is 0.000 out to 30
  s - and it costs nothing where the component is actually needed,
  because `exp(-x^2 / 2s^2)` kills the far tail at any scale: at 0.2 s
  it holds 76% of its peak density at 0.15 s and 46% at 0.25 s, against
  85% and 66% for the half-t. The slow tail now belongs to the decision
  family, which is what
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)’s
  `shape` and
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)’s
  `sigmadrift` are for.

  A welcome side effect: the `poutlier -> 1` degenerate mode is now
  thousands of log-likelihood units below the sensible one rather than
  hundreds, and `ndt` and the decision parameters no longer drop out of
  the density there, because a half Normal cannot explain a slow
  response at all. The mode still has infinite volume in `poutlier`
  itself, so
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  is still not optional.

- **`minrt` is removed** from every family, density, RNG, `*_stanvars()`
  and `*_lpdf_expose()`. The package works in **seconds**, full stop.
  The equivariance `minrt` bought in the likelihood was already
  fictional end to end:
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  shifted only the `ndt` prior with it, while the `sigmandt` prior of
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md),
  the `sigmadrift` prior of
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  and the `mu` priors are stated in seconds outright - and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  is not optional. Making the assumption explicit costs one argument
  from about twenty signatures and removes a whole class of
  misconfiguration.

  Calls passing `minrt` now fail with R’s usual `unused argument` error.
  Data in another unit fails **silently**, as it always did when `minrt`
  was left at its default: the outlier component’s log-density at
  `RT = 400` is about `-2e6`, so it contributes nothing, the mixture
  collapses to the unmixed shifted family, `poutlier` goes to zero and
  `ndt` is pinned by the fastest observed response. Nothing errors and
  the chains still initialise. Divide by 1000 before fitting.

  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  accordingly gives `ndt` a fixed `normal(-1.2, 0.2)`, and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  starts it at 0.1 s.

- [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  gains **`sigmadrift`**, the between-trial SD of the drift rate, so the
  Wald can produce the long right tails empirical RT distributions have.
  Described in full under 0.2.0 below.

- New function
  **[`pcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)**,
  the diffusion’s cumulative distribution function - the package had
  [`pcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`pcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  but no DDM counterpart. `response = NULL` gives the RT distribution
  marginally over the choice, and `response = 0`/`1` the defective CDF
  that boundary carries.

  It is the series that
  [`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  now inverts to draw from, so it was written and validated anyway;
  exposing it costs nothing and it is the more accurate of the two
  available implementations. Against numerical integration of
  [`dcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  over a grid of 360 cells the worst deviation is 1e-15, where
  [`rtdists::pdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
  is out by up to 5e-4 (and `pdiffusion(Inf, ...)` by 1.4e-4). It covers
  the classic 4-parameter DDM plus the outlier component; the
  between-trial variability parameters would each need their own
  quadrature layer, so they are not arguments rather than being silently
  ignored.

### Performance

- **[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  samples about 1.5x faster.** It was the most expensive likelihood in
  the package - roughly 5.9 us per observation per gradient, against 2.6
  for
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  and 0.8 for
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md) -
  and two thirds of that was avoidable. The Stan survival evaluated the
  normal CDF at `alpha` and at `beta` twice each, once inside `log_g()`
  and once again for the terms that follow it, and it assembled its six
  signed terms through a helper returning a `vector[2]`, which allocates
  on the autodiff stack once per term per observation per leapfrog step.

  `log_g()` now takes `log Phi(u)` from its caller, which brings the
  survival from six normal-CDF evaluations to four, and the six terms
  are grouped into the two differences that are individually monotone in
  the threshold, so each is one `log_diff_exp()` of known sign and
  nothing on the hot path returns a vector. The maths is unchanged and
  the grouped form cancels *less*: against the R implementation it
  agrees to 6e-11 where the term-by-term form reached 5e-9. Gradients
  are unchanged to the precision finite differences can resolve.

- **[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  is about 2x faster for the choice+RT families.** brms calls it once
  per observation, so anything done per call is paid thousands of times,
  and for every family except
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  the sampling itself was the small part: two
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html) constructions
  and an [`as.matrix()`](https://rdrr.io/r/base/matrix.html) came to
  roughly 63% of the call, against 13% for the actual draws. The
  registry’s `rng` entries now return a bare `list(rt, response)` and
  `.rchoice()` can return a matrix directly, so the prediction path
  builds no data frame at all. A
  [`cbind()`](https://rdrr.io/r/base/cbind.html) costs about a
  seventeenth of the
  [`data.frame()`](https://rdrr.io/r/base/data.frame.html) it replaces.

  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  [`rcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md),
  [`rcogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  and
  [`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  still return a data frame, with the same draws for the same seed -
  only the internal path changed. Measured on the
  `vignette("decision_making")` models at 20 draws: LNR 41.5 -\> 20.7 us
  per observation-draw, LBA 32.9.

- **[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  samples its predictions 4-6x faster**, on top of the above. It was the
  slowest family to predict from by an order of magnitude, because
  [`brms::rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
  takes one parameter set per call and roughly 85% of that call is fixed
  setup - which cannot amortise when every posterior draw carries its
  own parameters. (rtdists’ marginal cost is 1.4 us per draw against a
  772 us fixed cost; RWiener’s sampler does not amortise at all, staying
  at 55-90 us per draw for any n.)

  Draws are now taken by inverting the CDF, which vectorises across
  parameter sets because every step acts on the whole vector at once.
  The large-time series splits into a part that depends on the time and
  a part that does not, so the latter is built once and each evaluation
  is a single [`exp()`](https://rdrr.io/r/base/Log.html) and a column
  sum; and the density falls out of the same
  [`exp()`](https://rdrr.io/r/base/Log.html), which makes a Newton step
  cost exactly what a bisection step costs. Eight Newton passes leave
  about 0.6% of draws for a bisection cleanup that cannot fail, since
  the bracket is valid by construction.

  Accurate to 1e-13 in log RT against a 60-step bisection, with no
  Kolmogorov-Smirnov failure against either
  [`brms::rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
  or
  [`rtdists::rdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
  over a wide parameter grid. The CDF underneath it agrees with
  numerical integration to 1e-9, where
  [`rtdists::pdiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
  is out by up to 5e-4. On the `vignette("decision_making")` models: DDM
  560 -\> 122 us per observation-draw, DDM-5 483 -\> 75.

### Documentation

- Every help page now documents what each of its functions returns,
  rather than only the random-generation function it is named after -
  the `brms` family object, the `stanvars`, and the shape of the
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  output, including the families whose
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  errors because the decision time has no finite mean.

- Examples that were commented out now run. The plots are live, the
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  formulas are built, and the `*_lpdf_expose()` and model-fitting
  snippets are in `\donttest{}` behind a check for `cmdstanr` and a
  CmdStan installation instead of `\dontrun{}`, so they execute wherever
  the toolchain is present.
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md),
  [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md)
  and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  gained real examples.

- `DESCRIPTION` cites the papers behind the models.

## cogmod 0.2.1

### New features

- **`sigmabias = 0` is now allowed in
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  and
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)**,
  and in the single-accumulator case it is the **recinormal**, or LATER,
  model of [Carpenter & Williams
  (1995)](https://doi.org/10.1038/377059a0): the accumulator starts at
  zero on every trial, so the decision time is `boundary / drift` and
  `1 / (RT - ndt)` is normally distributed, with `mu` and `sigma` the
  mean and SD of *promptness*. Zero was previously rejected as an
  invalid parameter; it is a nested model, and the bound is now closed
  for the same reason
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)’s
  `sigmadrift` and
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)’s
  three between-trial variabilities are.

  ``` r

  bf(rt ~ 1, sigmabias = 0, boundary = 1)  # free: mu, sigma
  ```

  Nothing about the density had to change for this to be exact - at
  `sigmabias = 0` the existing Taylor branch evaluates to
  `dnorm(b / t, mu, sigma) * b / t^2 / pnorm(mu / sigma)` to machine
  precision, in R and in Stan alike - so the two models share one
  likelihood and
  [`loo_compare()`](https://mc-stan.org/loo/reference/loo_compare.html)
  between them is like-for-like.

  This matters beyond nesting. Estimating `sigmabias` freely is
  treacherous precisely *because* the recinormal limit is smooth: the
  likelihood goes flat as the start-point range shrinks, and `softplus`
  reaches zero only at minus infinity, so
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  has to fence the direction off. Pinning it at zero removes the
  parameter instead, which is the honest option when the design cannot
  identify a start-point range.

  Note that **two** pins are now needed for the evidence scale, not one:
  scaling multiplies every member of the scale ray by a common constant
  and leaves zero at zero, so `sigmabias = 0` drops off the ray rather
  than pinning it.
  [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)
  says so explicitly when that is the only fix present.

- New family
  **[`cogmod_exwald()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exwald.md)**
  ([Schwarz, 2001](https://doi.org/10.3758/bf03195403)): the decision
  time is a Wald convolved with an exponential residual stage of mean
  `tau` - the mechanistic counterpart of
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md),
  whose first stage is a descriptive Gaussian instead, with `tau`
  meaning the same thing in both. The mean exists and is
  `ndt + boundary / mu + tau`, so
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  returns a number.

  The density has two branches, **both exact**. Where `mu^2 > 2 / tau`
  the convolution collapses to a closed form in the Wald CDF; below
  that - which is the common regime, since at a drift of 3 and a
  threshold of 0.5 the closed form needs `tau > 0.22 s` - the same
  expression continues analytically through the Faddeeva function,
  giving `g * exp(-(boundary - mu * t)^2 / (2 * t)) * Re[w(z)]`. The
  exponent is the Wald’s own, so nothing overflows, and the branches
  meet exactly at `mu^2 = 2 / tau`. Across a grid spanning the usual RT
  region the density integrates to 1 to within 8e-12, the mean is right
  to 3e-9, and the relative step across the branch seam is 5e-8.

  Note there is deliberately no `sigmadrift`: it and `tau` both fatten
  the right tail and are very hard to tell apart, and
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  is where the drift-variability route lives. `ndt` and `tau` also share
  a ridge - both delay the response, and only the shape of the leading
  edge separates them - so
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  gives `tau` the same `normal(-1.5, 0.7)` it gives
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md).
  Fixing `ndt = 0` in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  recovers Schwarz’s own model.

- New family
  **[`cogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md)**,
  the Birnbaum-Saunders or fatigue-life distribution ([Birnbaum &
  Saunders, 1969](https://doi.org/10.2307/3212003)): a
  first-passage-time model in which evidence arrives in **discrete
  cycles and only ever towards the boundary** - what is random is the
  size of each increment, never its sign. Summing those increments and
  applying the CLT, then treating the cycle count as continuous, gives
  the first-crossing time.

  It is parameterized mechanistically, as `mu` (drift) and `boundary`
  (threshold), so it sits directly alongside
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  with the parameters meaning the same thing and only the mode of
  accumulation differing. Fixing the per-cycle SD at 1 is the same
  convention that fixes the Wald’s diffusion coefficient, and it makes
  `(mu * t - boundary) / sqrt(t)` **exactly** standard normal - the
  usual `(1 / a) * (sqrt(t / b) - sqrt(b / t))` written in these
  parameters, with `b = boundary / mu` and
  `a = 1 / sqrt(mu * boundary)`. The map between the two is a bijection,
  so nothing is given up.

  Everything is then closed form, and the density is the Wald’s own
  tilted by `(mu * t + boundary) / (2 * boundary)` - one sign apart from
  it. That tilt makes the family an **equal mixture of an inverse
  Gaussian and its length-biased twin**, so at the same `(mu, boundary)`
  it is slower and more dispersed than the Wald (mean 0.222 s against
  0.167, SD 0.184 against 0.136 at `mu = 3, boundary = 0.5`), while
  keeping the same exponential-order right tail.
  `E[T] = ndt + boundary / mu + 1 / (2 * mu^2)` is always finite, so
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  returns a number, and the median is exactly `ndt + boundary / mu`.
  There is no `sigmadrift`: the extra dispersion comes from the mixture
  structure at no cost in parameters, and drift variability is what
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  is for.

  It is also the cheapest first-passage density in the package - one log
  and one square, no branch and no special function - and
  [`rcogmod_bisa()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_bisa.md)
  is one normal draw per observation, exact, with no rejection step.

- New family
  **[`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md)**:
  `log(RT - ndt)` follows a Student-t, a robust LogNormal that varies
  kurtosis where
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  varies skew. The heavy tail absorbs slow contaminants into the
  likelihood rather than into a mixture component, which matters because
  the `poutlier` component is a half Normal and by construction cannot
  explain a slow response. At `dof = 5` the probability of a decision
  time beyond 5 s is about five orders of magnitude larger than the
  matching LogNormal’s.

  The degrees of freedom are called **`dof`**, not `nu`:
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  already spends `nuzero`/`nuone` on drift rates, and `brms` recognises
  the name `nu` and supplies defaults of its own for it.

  Two things to know. **The mean does not exist** for any finite `dof`,
  so
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  errors rather than returning a number; the median is exact at
  `ndt + exp(mu)`. And **the density is unbounded at `ndt`** -
  integrable, so the posterior stays proper, but the likelihood has no
  maximum, which is one more reason
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  is not optional. A Student-t is also symmetric on the log scale, so a
  small `dof` fattens the fast tail as well as the slow one and competes
  with `poutlier`;
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  centres `dof` at 6 with 95% of its mass between 1.5 and 24 to keep
  that in check.

- [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now supports
  **[`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)**,
  where before it returned the `brms` defaults with a message. `sigma`
  and `tau` are both lengths of time in seconds behind a `softplus`
  link, which `brms` has no way to know: `tau` arrives flat, and `sigma`
  arrives with the `student_t(3, 0, 2.5)` that `brms` supplies because
  it recognises the *name* - a Gaussian SD centred on 0.69 s modelled,
  1.9 s omitted, wider than most whole RT distributions. They now get
  `normal(-2.3, 0.7)` and `normal(-1.5, 0.7)` on the link scale (roughly
  25-330 ms for `sigma`, 55-630 ms for `tau`), and the matching
  `lognormal` when the dpar is left out of
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  altogether. `mu` is deliberately untouched: it is the response’s own
  intercept and the `brms` default is already proper and reasonable
  there.

- [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  gains **`sigmadrift`**, the between-trial SD of the drift rate. Above
  zero, each trial draws its own drift from a `Normal(mu, sigmadrift)`
  truncated at zero, which is what lets the Wald reach the long right
  tails empirical RT distributions have. Marginalising over that draw is
  a Gaussian integral, so the density stays closed form and costs two
  normal CDFs; `sigmadrift = 0` gives back the previous density exactly,
  not approximately.

  The truncation is what keeps the density proper: a single-boundary
  accumulator handed a negative drift never terminates, so an
  untruncated Normal would leave up to a third of the mass unaccounted
  for.
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)’s
  `sigmadrift` needs no such truncation, a diffusion between two
  boundaries always absorbing at one of them.

  It is fixed the same way as the
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  variability parameters - writing `sigmadrift = 0` in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  pins it and recovers the classic Wald, while leaving it out of
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  *estimates* it. Fixing it is the better default: `sigmadrift` and
  `poutlier` both fatten the right tail and are only weakly
  distinguishable, and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  gives `sigmadrift` a deliberately informative prior where it is
  estimated. Note that with `sigmadrift > 0` the density decays as
  `t^-2` and **the mean does not exist**, so
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  returns `Inf`; summarise
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  draws instead.

  Two consequences for existing code. The
  `drift`/`boundary`/`ndt`/`poutlier` functions gained an argument, so
  `sigmadrift = 0` now sits between `ndt` and `poutlier` in the
  signatures of
  [`rcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  [`dcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  and
  [`pcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  (positional calls that passed `poutlier` fourth need updating; named
  calls are unaffected). And a formula that does not mention
  `sigmadrift` at all now estimates it rather than fitting the
  fixed-drift Wald.

- New
  [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md):
  the third of the three generics that take the model rather than the
  family, alongside
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md).
  It reads the family off the formula and returns that family’s Stan
  code, so the family is named once - in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) -
  instead of three times:

  ``` r

  f <- bf(RT ~ Condition, ndt ~ Condition, family = cogmod_lognormal(minrt = 0.25))
  brm(f, data = df,
      prior    = cogmod_priors(f, df),
      init     = cogmod_inits(f, df),
      stanvars = cogmod_stanvars(f))
  ```

  This is not only tidier. `minrt` is baked into the generated Stan code
  as a literal, because a Stan function cannot see the data block, so
  `cogmod_lognormal(minrt = 0.25)` fitted with
  [`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  runs happily against an outlier component the family does not
  describe.
  [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)
  takes `minrt` off the family, and the two cannot disagree. The
  per-family `<family>_stanvars()` functions are unchanged.

- [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  now supports
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md),
  whose three parameters are all on the RT scale behind a `softplus`
  link and so are equally badly served by the default start at
  `log(2) = 0.69` s - which makes the Gaussian SD alone wider than most
  whole RT distributions.

- [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  now returns a value for **every** parameter the Stan program declares,
  not only the ones it has an opinion about, so CmdStan no longer prints
  `Init values were only set for a subset of parameters` and lists the
  rest. Regression slopes and standardized group-level effects start at
  zero, group-level and spline SDs just above zero, Cholesky factors at
  the identity; all are at least as good a starting point as Stan’s own
  `U(-2, 2)`.

  Two related fixes come with it. Links are now read off the family
  rather than the registry, so `cogmod_gamma(link_mu = "log")` is
  honoured. And the jitter is applied on the unconstrained scale -
  additive when free, multiplicative for a positive parameter, on the
  logit scale for a bounded one - so a jittered start can no longer land
  outside its own bounds, which it could previously for a dpar left out
  of the formula and estimated on the natural scale.

- New
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md):
  starting values for the families that estimate `ndt` directly. `brms`
  initialises on the unconstrained scale, so `init = 0` puts `ndt` at
  `exp(0) = 1` second - above most sub-second RTs, which leaves every
  response attributed to the outlier component and the decision
  parameters with no gradient at all. For
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
  and
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
  it *also* puts the shape at `softplus(0) = 0.69`, below the 1 at which
  the density becomes unbounded at the shift. No single scalar avoids
  both, since the two pull in opposite directions.

  On 1500 simulated Gamma trials (true shape 3, true `ndt` 0.25),
  `init = 0` left the shape stuck at its starting value with `Rhat` 2.3
  and an ESS of 3 after 306 s;
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  recovered shape 3.23 and `ndt` 0.227 with `Rhat` 1.01 in 28 s. An
  informative prior on the shape did not rescue `init = 0` - a prior
  cannot move a chain whose gradient is zero.

  Parameter names are read off
  [`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html)
  for the model actually being fitted, so `0 + Intercept`, interactions,
  group-level terms and smooths are all handled; anything it does not
  recognise is left to Stan.

- New
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  family: a shifted Log-Gamma model for reaction times, equivalently a
  shifted generalized gamma. `log(RT - ndt)` follows a location-scale
  log-gamma with location `mu`, scale `sigma` and shape `shape`, and
  `ndt` / `poutlier` / `minrt` work exactly as in
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md).

  `shape` is unconstrained, with `shape = 0` in the interior: it
  recovers
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  exactly, `shape = sigma` the shifted Gamma, `shape = 1` the shifted
  Weibull and `shape = -1` the shifted inverse Weibull. Fitting it is
  therefore a way of testing whether the LogNormal shape is adequate -
  an interval for `shape` covering 0 says it is. Negative `shape` gives
  a power-law right tail.

  Note the boundary at `sigma * shape >= 1`, where the decision density
  becomes unbounded at `ndt` and the likelihood with it;
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  sets `normal(0, 0.5)` on the `shape` intercept to keep well clear of
  it.

  **Fit this family with `init = 0`.** The prior keeps the posterior
  clear of that boundary but not the starting point: `brms` initialises
  from `U(-2, 2)` on the unconstrained scale, so about 15% of chains
  start with `sigma * shape >= 1`, fall into the spike at `ndt` and
  never finish. `init = 0` starts every chain at `shape = 0`, the
  LogNormal, and removes the problem.

- [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`without_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md)
  and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now work on
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  as well as
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md).

- [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)
  now **warns when the evidence scale is left free** for
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  and
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md).
  Both have a likelihood that is *exactly* constant along the ray that
  multiplies the drift rates, their SDs, the start-point range and the
  threshold offset by a common factor - verified to machine precision,
  not merely near-flat - so nothing in the data can pick a point on it.
  The failure is quiet rather than loud: the fit converges,
  [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html)
  looks right, and only the individual parameter estimates are
  meaningless, being whatever the priors happen to say about that
  direction.

  Fixing any one member of the ray in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  pins it and silences the warning - `sigmazero = 1` conventionally, as
  `rtdists` and `EMC2` both do, but `boundary = 1` or `sigmabias = 0.5`
  work as well. Note that leaving a parameter *out* of
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  does not fix it: `brms` declares it as a free auxiliary parameter and
  the ray stays exactly as free, which is the case the warning mostly
  exists to catch.
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  are quiet by construction, their unit diffusion coefficient having
  pinned the scale already.

### Bug fixes

- [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  no longer reports `Non-finite gradient` during warmup or a Pathfinder
  search, and no longer collects the divergent transitions that come
  with it. Stan’s classic 4-parameter `wiener_lpdf()` - much the fastest
  of the three Wiener densities Stan offers, and the one this family
  used whenever the three between-trial variability parameters were
  zero - returns `-inf` in two regions, and hands back **NaN** partial
  derivatives when it does. A NaN partial is not made harmless by the
  mixture weight on it being zero: reverse-mode multiplies the (zero)
  adjoint into the stored partial, and `0 * NaN` is `NaN`, so a single
  trial in one of those regions turns the gradient of the whole model to
  NaN, and Stan rejects the proposal.

  The two regions are the alternating small-time series losing its sum
  to cancellation - which depends only on the rescaled decision time
  `tau = t / boundary^2`, not on the drift or the scale separately, and
  sets in below `tau = 6.6e-4` - and the density underflowing to zero,
  which a cheap leading-term estimate detects. Those calls now go to
  Stan’s `sv`-capable density instead, which works in log space
  throughout and stays finite, with finite gradients, down to
  log-densities of `-1e7`. The two agree to `1e-13` where the paths
  meet, so there is no step in the likelihood, and the fast path still
  handles the overwhelming majority of evaluations.

  On the 2000-trial fit in
  [`?cogmod_ddm`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md),
  over three seeds: 16 `Non-finite gradient` reports per Pathfinder run
  became 0, and 47-420 divergent transitions per NUTS run became 0.
  Sampling is about 35% slower per iteration and roughly two to a
  hundred times *better* per effective sample.

  The same guard is applied to the general 7-parameter form, used when
  `sigmabias` or `sigmandt` is nonzero, which fails the same way. There
  it returns the `-inf` that form would have returned anyway, but as a
  constant, which carries no partial derivatives.

- The test suite runs in a third of the time (2472 s to 968 s on
  Windows). Every family’s Stan `lpdf` now goes into **one** model,
  compiled once per session, instead of nine separate `*_lpdf_expose()`
  compilations; the factorial parameter sweeps are thinned to subsets
  that still cover every level of every factor; and the six
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  model fits are behind `COGMOD_TEST_SLOW`, which the CI workflow sets.
  Setting it locally runs them:

  ``` r

  Sys.setenv(COGMOD_TEST_SLOW = "true"); devtools::test()
  ```

  The shared model also removes a trap the RDM tests had worked around
  with a cache of their own:
  [`expose_functions()`](https://paulbuerkner.com/brms/reference/expose_functions.brmsfit.html)
  fails on a model
  [`cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
  returns pre-compiled, so a second call in the same session errors
  rather than reusing the first.

- [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now covers
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)’s
  `nuone`, `sigmazero` and `sigmaone`, which it previously left flat.
  Push an accumulator’s rate down far enough and it stops finishing
  first ever; the density then depends on it only through the loser’s
  survival term, which has already saturated at 1. Past about
  `nuone = -6` the log-likelihood is *exactly* constant, and that
  accumulator’s `sigma` is unidentified along with it. The outlier
  component makes this reachable rather than hypothetical: it floors the
  trials the retreating accumulator can no longer explain, so the
  plateau is there even when both responses are observed.

  `nuone` gets `normal(0.7, 1.5)` on its identity link, the two sigmas
  `normal(0, 1)` on softplus (`lognormal(-0.7, 0.75)` when omitted from
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) and
  so on the natural scale), and all three `normal(0, 0.5)` on slopes.
  `mu` - which is `nuzero` - has the mirror-image plateau, but it is the
  response’s own intercept and `brms` already gives it a proper
  `student_t` default, so it is left alone; if you model a rarely-chosen
  option it is worth mirroring the `nuone` prior onto it by hand.

  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  already covered this family and is unchanged.

- [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now covers
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)’s
  `sigmabias` and `boundary`, which it previously left on `brms`’s flat
  default. As the start-point range approaches zero the LBA converges
  smoothly to the recinormal, so the likelihood stops depending on
  `sigmabias` altogether - and a `softplus` link reaches zero only at
  minus infinity. That is a flat prior over an infinite flat region, the
  same improper posterior the function already exists to prevent for
  `ndt` and `poutlier`, and it failed just as quietly: on the 4285-trial
  fit in `vignette("rt_models")`, `sigmabias` for one condition ran off
  to `softplus(-10.4) = 3e-05` with `Rhat` 1.69 and an effective sample
  size of 6, while every other parameter looked healthy. With the priors
  in place the same fit gives `Rhat` 1.02, an effective sample size of
  387, no maximum-treedepth hits, and a finite estimate. `boundary` is
  covered too, since `b = sigmabias + boundary` puts the two on the same
  ridge.

  Both get `normal(0, 1)` on the link scale, `lognormal(-0.7, 0.75)`
  when the dpar is omitted from
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) and
  so lives on the natural scale, and `normal(0, 0.5)` on slopes - wider
  than the blanket `normal(0, 0.2)` the other dpars take, because the
  point is to fence off zero rather than to shrink effects. Families can
  now declare such rows in the registry, so the next one that needs them
  does not need a special case.

- The
  [`?rcogmod_weibull`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
  and
  [`?rcogmod_gamma`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
  notes about the shape were understated. They said the density is
  unbounded at `ndt` for a shape below 1, which is true, but the
  threshold that matters in practice is **2**: below it the derivative
  of the log-likelihood with respect to `ndt` is unbounded at every
  observation, so the posterior stays proper while the sampler grinds.
  On the data in `vignette("rt_models")` the Weibull shape comes out at
  1.4, `ndt` lands inside the dense left edge of the data, the step size
  collapses to 0.005 against 0.19 for
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  and mean treedepth goes from 3.9 to 8.1 - 19x the gradient evaluations
  and 19x the wall time, with `Rhat` 1.18 on `ndt`. The density is
  cheap; all of the cost is geometry. Both help pages now set out the
  three regimes and say to prefer
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md),
  which nests the Weibull at `shape = 1`, when the shape comes out below
  2.

  No prior is set on those shapes, deliberately, and none is set on
  `ndt` either. Every obvious remedy was tried on that fit and measured:

  - `normal(2.4, 0.4)` on the shape (softplus scale, 95% of its mass
    above a shape of 1.9) moved the posterior shape by 0.01. The
    likelihood prefers the low-shape corner by around 100 log units; the
    prior contributes 5.
  - `normal(-1.25, 0.05)` on `ndt`, centred below the fastest bulk
    response, left the posterior 6.5 prior SDs away at essentially its
    unconstrained value. The `ndt` likelihood has a posterior SD of
    0.003, some fifteen times sharper than that prior. The attempt cost
    4% divergent transitions against 0.5%, 16% of iterations at maximum
    treedepth against 7%, `Rhat` 1.43 against 1.18, and a slightly worse
    `loo`.
  - Fixing `ndt` at the fastest observed response does remove the
    problem, by removing the parameter - and reinstates the min-RT bound
    this parameterization exists to remove. On these data the fastest
    response is 71 ms, which is not a decision.

  Note what is *not* wrong: `ndt` and the shape are jointly identified,
  sharply (posterior SD on `ndt` of 3 ms), so this is not two parameters
  trading off with nothing to separate them and pinning one is not the
  missing ingredient.

  The slow sampling and a poor fit turn out to be the same fact. Across
  the ten families fitted in `vignette("rt_models")` the Weibull comes
  **last** by `loo` - 196 elpd (SE 21) behind
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md),
  and 95 behind the next worst. What the sampler struggles with is the
  model contorting itself to represent a left edge it cannot otherwise
  reach. Where the shape comes out above 2 the family is fine, as
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
  is on these same data at 2.2; a shape below 2 is best read as the
  model asking for a different one.

  Under the older `ndt = tau * min(RT)` parameterization the same
  singularity was damped by the logit Jacobian vanishing as `tau`
  approached 1, which is why it only became visible once `ndt` was
  estimated directly.

- [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  no longer leaves `brms` warning that a global `b` prior “will not be
  used in the model as all related coefficients have individual priors
  already”. It set both the blanket row and the row for the coefficient
  named `Intercept`, which is right when there are slopes beside the
  intercept but leaves the blanket row covering nothing under
  `ndt ~ 0 + Intercept`. The pair is now checked in both directions.

- [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  also carried hand-written R copies of the shared mixture - their own
  `r*()`, `d*()`, `.prepare_*()`, `log_lik_*()`, `posterior_predict_*()`
  and `posterior_epred_*()`, about 420 lines reproducing what the other
  seven families delegate. All of it now delegates too, verified
  bit-identical to the code it replaced before the change was made.

- [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  generated their Stan code from hand-written copies of the shared
  mixture rather than from the registry every other shifted family uses.
  The copies had drifted: the LogNormal one described its outlier scale
  as `0.8 * minrt` when it is `minrt`, and neither carried the folded
  normalising constant that makes the outlier term ~1.4x cheaper per
  gradient evaluation. Both now come from the one generator, verified
  against the R density to machine precision.

- [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now also covers dpars left out of
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  entirely. `brms` declares those as plain auxiliary parameters - class
  `"<name>"` with an empty `dpar`, on the natural scale with no link -
  so matching on `dpar` alone missed them and they kept `brms`’s own
  defaults. Those defaults are actively wrong here: `uniform(0, min_Y)`
  on `ndt` reimposes the very min-RT bound the parameterization exists
  to remove, `gamma(0.01, 0.01)` on `shape` has support only on the
  positives and would silently truncate away the inverse-Weibull half of
  the family, and `poutlier` was left flat over `[0, 1]` with half its
  mass above 0.5. All three are now replaced, with priors on the natural
  scale: `lognormal(-1.2, 0.2)`, `normal(0, 0.5)` and
  `exponential(100)`.

  The `poutlier` prior for the omitted case has its mode at **zero**,
  unlike its modelled counterpart: leaving it out of the formula is
  taken to mean the data were trimmed or no outliers are expected. Its
  median is unchanged.

### Breaking changes

- **Every family is renamed to a single `cogmod_*` scheme.**
  `rt_lognormal()` is
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  `lnr()` is
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  `ddm()` is
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md),
  and so on; the two LBAs are told apart by their accumulator count, so
  `rt_lba()` (no choice) is
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  and `lba()` (two choices) is
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md).
  Every derived function follows its family: `rrt_lognormal()` is
  [`rcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  `rlnr()` is
  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  `rt_lognormal_stanvars()` is
  [`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  and likewise for the densities, the `*_lpdf_expose()` and the `brms`
  post-processing hooks.

  Two things were wrong with the old names. `lba()`, `lnr()`, `ddm()`,
  `rdm()` and `choco()` are generic enough to collide with anything else
  attached, and the `rt_` prefix said “reaction time” on families that
  are not all RT-only. One prefix fixes both, and `cogmod_`
  tab-completes the whole package.

  **The old names all still work** in this version, as exact synonyms
  rather than wrappers - the `brms` hooks included, so a model fitted
  before the rename can still be summarised. (They were removed in
  0.3.0.)

- **The `bs` dpar is renamed `boundary`** in
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md),
  [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md),
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md).
  This is not cosmetic: `brms` names the unpenalized spline coefficients
  of a model `bs`, so any of those families combined with a smooth
  failed to compile with `Identifier "bs" is already in use`. There is
  no alias for a dpar name - update `bs ~ ...` to `boundary ~ ...` in
  formulas, `bs =` to `boundary =` in the `r*()`/`d*()` calls, and
  expect `boundary_Intercept` where you had `bs_Intercept` in the
  output.

  Fits made before this change cannot be post-processed, because their
  draws carry the old dpar name; refit them.

- **`cogmod_lnr(link_nuzero = )` is now `link_mu`.** `brms` requires the
  first distributional parameter of a custom family to be called `mu`,
  so `mu` is what the formula uses and what the argument now matches;
  `nuzero` remains the name of the quantity, in the prose and in
  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)/[`dcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md).

- **The `r*()` and `d*()` functions of every shifted family now validate
  their decision parameters**, against the same bounds the generated
  Stan code checks. Only
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  did this before; the other seven silently returned zero density or
  `NaN` draws for, say, a negative `sigma`. `r*()` errors, `d*()` warns
  and returns zero - so a single bad posterior draw still cannot abort a
  whole call.

- **[`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`cogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md)
  and
  [`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md)
  now use the same parameterization as
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)**:
  `tau` and the `minrt` *dpar* are gone, replaced by `ndt` estimated
  directly (log link) plus `poutlier`, with `minrt` carried on the
  family as a constant. Every one of them gains an outlier component,
  and
  [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`without_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md)
  and
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now work on all eight families.

  Update models as for
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md):
  replace `tau ~ ...` with `ndt ~ ...`, drop `minrt = min(df$RT)`, and
  add `poutlier ~ 1`. `ndt` coefficients are on the log scale.

  The `d*()`/`r*()` functions gain `poutlier` and `minrt` arguments, and
  [`rcogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`dcogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`rcogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`dcogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`rcogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`dcogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`rcogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
  [`dcogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
  [`rcogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md)
  and
  [`dcogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md)
  are new - those families previously had no R-level density or RNG at
  all.
  [`pcogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)
  now returns the mixture CDF.

  Note that
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md)
  and
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md)
  have an unbounded density at `ndt` whenever their shape falls below 1,
  which the outlier component cannot repair; fit them with `init = 0`.
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  nests both and lets the data choose the shape instead.

- **[`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)**
  moves to the same parameterization: `tau` and the `minrt` dpar are
  replaced by `ndt` (log link) plus `poutlier`, with `minrt` a family
  constant.
  [`rcogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  and
  [`dcogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  gain `poutlier` and `minrt`.
  [`posterior_epred_cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  now explains *why* there is no expectation rather than calling it
  prohibitive: the decision time is `(b - U(0, A)) / drift` with a drift
  truncated at zero, whose density is positive at 0, so `E[1 / drift]`
  diverges and the mean does not exist.

- All nine shifted families are now generated from a single internal
  registry, so the Stan code, the R density, the RNG, the likelihood and
  the predictions cannot drift apart. Verified family by family: Stan
  agrees with R to machine precision, each density integrates to 1, and
  each RNG reproduces its own density.

- **[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  moves to the same `ndt` + `poutlier` parameterization as the RT-only
  families.** `tau` and the `minrt` *dpar* are gone; `ndt` is estimated
  directly, on the log link, and `minrt` is a constant carried on the
  family object rather than something `tau` is scaled by.
  [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`without_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md),
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  all work on it now.

  Because the LNR produces a **choice as well as a time**, its outlier
  component is not quite the RT-only one: the contaminant guesses
  uniformly over the two response options in addition to drawing an RT
  from the same half Student-t, so the mixture is \`poutlier \* (1/2) \*
  g(t) + (1 - poutlier)

  - f_k(t -
    ndt)`. The`1/2`is what keeps the joint density summing to one over both responses - without it the total comes to`1 +
    poutlier\`.

  Update models as for
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md):
  replace `tau ~ ...` with `ndt ~ ...`, drop `minrt = min(df$RT)` from
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html),
  and add `poutlier ~ 1`. Fit with `init = cogmod_inits(f, df)` rather
  than `init = 0` - on the log link, `init = 0` starts `ndt` at
  `exp(0) = 1` second, above nearly every sub-second RT, which leaves
  every response attributed to the outlier component and the race
  parameters with no gradient at all.

  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  and
  [`dcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  gain `poutlier` and `minrt` arguments. Verified against the RT-only
  families’ checklist: the joint density sums to one over both responses
  and integrates to one over time (with and without the `1/K` term, to
  confirm it is load-bearing), the Stan `cogmod_lnr_lpdf` agrees with
  the R density to machine precision, and a simulated fit recovers `ndt`
  well above the fastest observed response.

- **[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  moves to the same `ndt` + `poutlier` parameterization**, on the same
  shared machinery as
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md).
  `tau` and the `minrt` *dpar* are gone: the family’s dpars are now
  `mu`, `driftone`, `sigmabias`, `boundary`, `ndt`, `poutlier`, `ndt` is
  estimated directly on the log link, and `minrt` is a constant carried
  on the family object.
  [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`without_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md),
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md),
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  and
  [`cogmod_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_stanvars.md)
  all work on it now.

  Update models as for
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md):
  replace `tau ~ ...` with `ndt ~ ...`, drop `minrt = min(df$RT)` from
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html),
  and add `poutlier ~ 1`. Fit with `init = cogmod_inits(f, df)` rather
  than `init = 0` or `init = 0.5`.

  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  also supplies the `sigmabias` / `boundary` priors that
  [`?cogmod_rdm`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  previously told you to write by hand: the two enter the model only
  through the sum `b = boundary + sigmabias` and trade off along a ridge
  worth a handful of log units, which under a flat prior and a
  `softplus` link is an improper posterior. The drift rates are left
  flat on purpose - unlike
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)’s
  `nuone`, a drift pushed to zero does not produce a plateau, because a
  driftless accumulator still finishes and still wins sometimes.

  [`rcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  and
  [`dcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  gain `poutlier` and `minrt` arguments, and
  [`pcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  gains `poutlier` and `minrt` too - its CDF is now the mixture’s, and
  still keeps the far-tail survival in log space.
  [`dcogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  keeps its `response = NULL` marginal, which is now the sum of the two
  defective mixture densities. Invalid parameters now warn and return a
  zero density rather than erroring, matching the other mixture
  families. **A drift of exactly zero is still accepted** - it is the
  one closed lower bound in either registry, because driftless Brownian
  motion still reaches the threshold with probability one.

  Verified against the checklist: the joint density sums to one over
  both responses and integrates to one over time for several parameter
  sets and `poutlier` values (including as the start-point range shrinks
  to `1e-6`), the Stan `cogmod_rdm_lpdf` agrees with the R density
  across the parameter grid, and a simulated fit recovers `ndt = 0.256`
  against a true `0.25` - some 200 times the fastest observed response,
  which the old `tau * minrt` bound could not have expressed.

- **[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  moves to the same `ndt` + `poutlier` parameterization.** Its dpars are
  now `mu`, `driftone`, `sigmazero`, `sigmaone`, `sigmabias`,
  `boundary`, `ndt`, `poutlier`. Update models as for
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md):
  replace `tau ~ ...` with `ndt ~ ...`, drop `minrt = min(df$RT)` from
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html),
  and add `poutlier ~ 1`.

  Two long-standing bugs in the density came out with it.

  **The density was not normalised.** A normal drift rate can come out
  negative, and such an accumulator never reaches the threshold, so a
  trial on which *both* drifts are negative produces no response at all.
  [`rcogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  has always resampled until at least one is positive - but
  [`dcogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  never divided by the probability of that event, so the density
  integrated to the probability rather than to one. At drift rates of
  `0.5` and `0.2` with SDs of `1.5` it came to `0.83`, and the simulated
  choice proportion was `0.562` against an integral of `0.468`. Because
  the shortfall depends on the parameters, it biased estimates rather
  than merely offsetting the likelihood. Both the R and the Stan
  densities now condition on the event the process is conditioned on.

  **The `(1 / A)` cancellation
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  was fixed for was still here.** The defective density divides \`drift
  \* (Phi(z2) - Phi(z1)) + sigma \* (phi(z1)

  - phi(z2))`by the start-point range, and both differences vanish linearly in it; the loser's survival was computed as`1 -
    CDF`, which cancels the same way. Relative error reached 2.7% at`sigmabias
    =
    1e-5`and 320% at`1e-7`. Both now go through the kernels`cogmod_lba1()`already uses (`.lba_dens_over_A()`, and a new`.lba_surv_raw()`that takes the survival directly rather than as`1 -
    CDF\`), which the two families now share in R and in Stan.

  The `.Machine$double.eps` floor is gone too. It turned every RT below
  the point where the density becomes representable into a log-density
  of exactly `-36.04` - a constant the model never produced, with a
  gradient of zero. Where the density really has underflowed the
  log-density is now `-Inf`, and the outlier component is what keeps the
  *mixture* finite there.

  [`rcogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
  now imposes the positive-drift condition **exactly**, by sampling
  which accumulator is positive and then the truncated normals, instead
  of a rejection loop. Its `max_iter` argument is therefore gone, along
  with the fallback that forced a drift positive with
  [`abs()`](https://rdrr.io/r/base/MathFun.html) - and so drew from the
  wrong distribution - whenever the loop ran out.

  Note that the evidence scale of an LBA is **arbitrary**: multiply the
  drifts, their SDs, the start-point range and the threshold by any
  `c > 0` and every finishing time is unchanged, so the likelihood is
  exactly constant along that ray.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now fences all four positive parameters, which makes the posterior
  proper, but only fixing one SD in the formula (`sigmazero = 1` in
  [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html))
  identifies the scale. This was true before and is now documented.

- **[`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  moves to the same `ndt` + `poutlier` parameterization**, and
  **`sigmatau` is renamed `sigmandt`**. Its dpars are now `mu`,
  `boundary`, `bias`, `sigmadrift`, `sigmabias`, `sigmandt`, `ndt`,
  `poutlier`.

  `sigmatau` was the between-trial range of the non-decision time
  expressed as a fraction of `minrt` (`st0 = sigmatau * minrt`). With
  `tau` and the `minrt` dpar both gone it was named after a parameter
  that no longer exists and scaled by a constant the user no longer
  sees, so it is now **`sigmandt`**, which is `st0` itself, in the same
  unit as the data, on a log link - `ndt` remains the lower bound of the
  resulting Uniform. Update models by replacing `tau ~ ...` with
  `ndt ~ ...`, dropping `minrt = min(df$RT)`, adding `poutlier ~ 1`, and
  rewriting any `sigmatau` term as `sigmandt` in seconds.

  All three between-trial variability parameters remain legitimately
  **zero**, and fixing them in the formula (`sigmadrift = 0`) still
  recovers the classic 4-parameter DDM.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  now supplies priors for all three: each has a floor at zero that its
  link reaches only at minus infinity, and the likelihood stops changing
  well before then, which under `brms`’s flat default is an improper
  posterior.

  [`posterior_epred_cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  keeps its closed form - the DDM is the one choice family here with a
  usable one - now using `ndt` directly and blending in the outlier
  component when `predict_outliers` is set, like the RT-only families.
  [`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  and
  [`dcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  no longer take `...` for
  [`brms::rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)/[`brms::dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html),
  which the shared mixture machinery cannot forward; set
  `options(wiener_backend = )` instead.

  Both R and Stan evaluate the decision component at a non-decision time
  of zero, which
  [`dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html),
  [`rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html) and
  `wiener_lpdf()` all refuse. Since the Wiener density depends on the
  time and the non-decision time only through their difference, both
  offset the pair by the same `1e-10`, and the Stan literal is generated
  from the R constant so the two cannot drift apart.

  All four choice+RT families are now on the direct `ndt` + `poutlier`
  parameterization, and `tau` + `minrt` is gone from the package.

### Bug fixes

- **[`dcogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
  was wrong for small start-point ranges.** Its density is built from
  `drift * (Phi(z2) - Phi(z1)) + sigma * (phi(z1) - phi(z2))` divided by
  the start-point range `A`, and both differences vanish linearly in
  `A` - so evaluating them directly and then dividing lost every
  significant digit once `A` was small. The old code also floored the
  bracket at `1e-10`, which turned that underflow into a spurious
  density floor spread over the whole line. The result: the density
  stopped integrating to one below about `A = 0.1` and was outright
  divergent below `A = 0.01`.

  Both differences are now computed stably - a Taylor expansion in
  `delta = A / (sigma * t)` below `1e-4`, and tail-aware differencing
  above it - and the floor is gone. The density integrates to 1 from
  `A = 2` down to `A = 1e-8`, and converges to the recinormal (LATER)
  limit at the expected first-order rate. The same fix is in the Stan
  code, which matches the R density across 960 parameter combinations.

### Bug fixes (parameterization)

- [`posterior_epred_cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md)
  returned `exp(mu + sigma * 0.5772)`, which is the *geometric* mean of
  the decision time - the exponential of `E[log(RT)]` - rather than its
  mean. It now returns `exp(mu) * gamma(1 - sigma)`, and `Inf` where
  `sigma >= 1` and no mean exists.

- [`posterior_epred_cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md)
  now returns `Inf` where the Frechet shape is `<= 1` and the mean does
  not exist, instead of a finite but meaningless value.

- [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  no longer uses `tau` and `minrt`. Non-decision time is now estimated
  directly as `ndt`, in seconds, through a log link, and the family
  gains `poutlier`, the proportion of trials generated by an outlier
  process rather than by the decision process.

  The old parameterization set `ndt = tau * minrt` with `tau` in
  `(0, 1)` and `minrt` injected as a constant, which capped `ndt` at an
  order statistic of the sample. Because that cap was shared across the
  whole dataset, a non-decision time larger than the fastest observed
  response was inexpressible - so any condition or participant whose
  true `ndt` exceeded the global minimum RT could not be recovered, and
  the misfit surfaced instead as spurious effects on the other
  parameters.

  What makes the direct parameterization tractable is the outlier
  component: a fixed half Student-t (scale `0.4`, `3` df) mixed in with
  weight `poutlier` keeps the density positive below `ndt`, turning the
  hard min-RT boundary into a finite cost and leaving the log-density
  smooth. No bound is taken from the data. The half-t is flat at the
  origin, so the fastest responses are not starved of density; its tails
  are heavy enough that supplying RT in milliseconds degrades rather
  than underflowing to zero; and its mean is finite, which
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  requires.

  Models must be updated: replace `tau ~ ...` with `ndt ~ ...`, drop
  `minrt = min(df$RT)`, and add `poutlier ~ 1`. Note that `ndt`
  coefficients are on the log scale, so
  [`exp()`](https://rdrr.io/r/base/Log.html) them for seconds.

- The outlier component is scaled by `minrt`, the fastest reaction time
  that could plausibly be a real decision, used directly as the half-t
  scale with no conversion factor. It defaults to `0.3` seconds, where
  the conditional accuracy functions in the *Outliers* article show
  responses sitting at chance across three paradigms. 61% of the
  component falls below `minrt` whatever value is chosen.

  It is a judgement about the task rather than a statistic: nothing is
  read off the sample, and `ndt` is not bounded by it.

  This matters because the shifted LogNormal is scale-equivariant on its
  own - multiply every response by 1000 and `ndt` comes back multiplied
  by 1000 - and a component pinned to the second would be the one thing
  breaking that. With millisecond data and a fixed scale, the outlier
  component sits many orders of magnitude below the decision density
  everywhere in the data, `poutlier` collapses toward zero and `ndt`
  reverts to being pinned by the fastest observed response, silently and
  without a warning. Setting `minrt` in the unit of the data makes the
  likelihood exactly equivariant instead, so millisecond data need
  `minrt = 300`.

  `minrt` is a constant, never estimated, and is deliberately **not** a
  dpar: `brms` has no notion of a default for one, so a dpar omitted
  from the formula is estimated rather than defaulted. It is carried on
  the family, and the matching Stan constant comes from
  [`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  which also accepts the family itself so the two cannot drift apart:

  ``` r

  fam <- cogmod_lognormal(minrt = 0.3)
  f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = fam)
  brms::brm(f, data = df, prior = cogmod_priors(f, df),
            stanvars = cogmod_lognormal_stanvars(fam))
  ```

  [`rcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  and
  [`dcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  gain a matching `minrt` argument. Models fitted before this change
  fall back to the default, so their predictions are unaffected.

- The `wagenmakers2008` dataset has been removed. Those data were
  supplied by the original authors for distribution in the `rtdists`
  package specifically, and no open licence covers them, so
  redistributing them here was not appropriate. They remain available as
  [`rtdists::speed_acc`](https://rdrr.io/pkg/rtdists/man/speed_acc.html);
  `rtdists` is in `Suggests`, and the vignettes and paper now
  reconstruct the same subset from it
  (`!censor & response != "error" & rt <= 2`), so previously reported
  results are unchanged.

- The experimental confidence signal detection model (`rconf_sdt()`,
  `dconf_sdt()`, `conf_sdt_stanvars()`, `conf_sdt_custom_family()` and
  its `brms` methods) is no longer exported. It was not ready, and the
  code is commented out in `R/conf_sdt.R` and `R/conf_sdt_brms.R`
  pending a rework.

- **Priors are now required.** `brms` assigns a flat, improper prior to
  the intercept of any custom-family parameter it does not recognise,
  which here means both `ndt` and `poutlier`, and the likelihood has two
  flat directions that a flat prior turns into an improper posterior:
  `poutlier` toward 1, where every response is attributed to the outlier
  component and `mu`, `sigma` and `ndt` drop out of the density
  altogether; and `ndt` toward 0, where the model reduces to an
  unshifted LogNormal and the gradient with respect to `log(ndt)`
  vanishes. The second is inherent to putting a positive shift on a log
  link and has nothing to do with the mixture, which is why a prior on
  `poutlier` alone is not enough. Symptom: intercepts around `1e14`,
  `Rhat` near 2 and an effective sample size of about 5, with no error
  raised.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  fills the gap.

### New features

- `cogmod_priors(formula, data)` fills in every prior `brms` would
  otherwise leave flat - for
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  the `ndt` and `poutlier` rows. It starts from
  [`brms::get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
  for the model in hand rather than guessing, so `0 + Intercept`
  formulas, interactions, group-level terms and smooths are all handled
  and a prior matching no parameter is impossible by construction. The
  result is passed through
  [`brms::validate_prior()`](https://paulbuerkner.com/brms/reference/validate_prior.html),
  so a malformed specification errors there with the offending row in
  view, and the return value is the complete prior table with a `source`
  column marking each row as `user` or `default` - print it to see
  exactly what the model will be fitted with. To change one, edit the
  row; [`c()`](https://rdrr.io/r/base/c.html) will not work for a slot
  the table already covers, because `brms` rejects two priors for the
  same slot.

  The family is read off the formula, so build it with
  `brms::bf(..., family = cogmod_lognormal())`. Any other family, or a
  formula carrying none, gets a message and the `brms` defaults
  unchanged, so the call is always safe to leave in a script.

- [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md)
  returns the posterior probability that each trial came from the
  outlier component rather than the decision process - the mixture
  responsibility, averaged over draws. Responses below `ndt` come out at
  1 and those in the bulk near 0, but the probability rises again in the
  far slow tail, where the half-t has heavier tails than the LogNormal;
  that is the mechanism behind the advice to filter implausibly slow
  responses before fitting. The responsibility is computed on the log
  scale, so it stays finite in tails where both components underflow to
  zero. It returns `rt` and `p_outlier` only; the `fast` column was a
  marginal median split that ignored any grouping in the model and is
  gone.

- [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  for
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  describe the **decision process alone** by default, as if `poutlier`
  were zero. For visualising effects the outlier component is a nuisance
  that pulls expected values toward its own mean (0.441 s) and adds a
  spike of implausibly fast draws; it is also a fixed regularizer rather
  than a claim about how guesses are distributed, so simulating from it
  means simulating from something the model does not assert. The
  likelihood is unaffected and is always the full mixture, so
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  and
  [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
  no longer describe the same distribution - a hand-rolled LOO-PIT check
  should use
  [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md).

- [`with_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md)
  restores the fitted mixture for prediction, and
  [`without_outliers()`](https://dominiquemakowski.github.io/cogmod/reference/with_outliers.md)
  returns to the default. The main use for the former is
  [`pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html):
  on untrimmed data the decision-only predictive has no fast spike to
  match the one in the data, which reads as misfit. The same flag can be
  set up front with `cogmod_lognormal(predict_outliers = TRUE)`.

  It is carried on the family rather than passed as an argument because
  `brms` does not forward extra arguments to a custom family’s
  prediction methods - `posterior_epred` reaches the method with `prep`
  and nothing else - and `insight`, `modelbased` and `marginaleffects`
  inherit that. Carrying it on the object is what makes it work through
  all of them.

- New *Outliers* article validating the outlier-mixture specification
  against the Illusion Game dataset, cross-checked against the lexical
  decision data of Wagenmakers et al. (2008) and the brightness
  discrimination of Ratcliff and Rouder (1998), both from `rtdists`. It
  also derives recommended priors for all four parameters from the rates
  it measures, and explains why they matter more under this
  parameterization than the last: `ndt` is no longer bounded above by
  construction, and `poutlier` has a degenerate region near 1 where
  `ndt` becomes unidentified.

## cogmod 0.1.0

### Models for subjective scales

- Beta-Gate
  ([`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md),
  [`rcogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md),
  [`dcogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md)),
  a reparametrised ordered beta model.
- Discrete Beta
  ([`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md)
  and friends) for Likert-type responses.
- Choice-Confidence, CHOCO
  ([`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
  [`rcogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
  [`dcogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md))
  for bipolar scales.
- Signal detection with confidence ratings (`conf_sdt()`).

### Models for decision making

- Lognormal race
  ([`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)),
  linear ballistic accumulator
  ([`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)),
  drift diffusion with optional across-trial variability
  ([`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)),
  and the racing diffusion model
  ([`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)).

### Models for reaction times alone

- A consistently parametrised set of shifted, right-skewed response
  distributions:
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md),
  [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md),
  [`cogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md),
  [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md),
  [`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md),
  [`cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md),
  [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md)
  and
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md).

### Data

- `wagenmakers2008`, lexical decision data from Wagenmakers et
  al. (2008).
- `badlm`, a simulated dataset in which two conditions share a mean
  reaction time while differing in shift, spread and tail weight.
