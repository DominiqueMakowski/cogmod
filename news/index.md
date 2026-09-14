# Changelog

## cogmod 0.3.2

### New features

- **[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  gains `sigmabias`**, the same between-trial start-point range as the
  LNR’s below: the decision time is the LogNormal multiplied by a
  `Uniform(1, 1 + sigmabias)` distance, which is the single-accumulator
  LBA with a LogNormal drift rate and its threshold offset pinned at 1.
  At `sigmabias = 0` - the default of
  [`rcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  [`dcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
  and
  [`pcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  placed after `ndt` and before `poutlier` as the Wald’s `sigmadrift`
  is, and the value to fix in the formula unless the design speaks to
  start-point variability - the family is the shifted LogNormal exactly
  as before, bit for bit and at the same cost. The CDF and survival
  `cens()` needs, the mean
  [`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
  reports (which gains a factor `1 + sigmabias / 2`) and the Stan
  functions all carry the range, and agree with quadrature over the
  start point to `1e-7` and with each other to `1e-10`. The kernels are
  shared with
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  which is now a race of two of these accumulators in code as well as in
  theory. As for the LNR, a formula that omits `sigmabias` estimates it,
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  fences the flat direction at zero, and fits made with earlier versions
  have to be refit; the two vignette LogNormal models were.
  [`cogmod_logstudent()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logstudent.md)
  and
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  do not get the parameter: the density needs a partial first moment of
  the rate distribution, which does not exist for a Student-t on the log
  scale and needs incomplete gamma functions for the log-Gamma.

- **[`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  gains `sigmabias`**, a between-trial start-point range: each
  accumulator now starts at `Uniform(0, sigmabias)` and runs to a
  threshold `1 + sigmabias` at its LogNormal rate, so its finishing time
  is the distance divided by the rate rather than the reciprocal of the
  rate alone. At `sigmabias = 0` - the default of
  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  and
  [`dcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md),
  and the value to fix in the formula unless the design can identify a
  start-point range - the family is the LNR exactly as before, bit for
  bit and at the same cost, since both the R and the Stan kernels take
  the plain lognormal branch there. Above zero it is the LBA with
  LogNormal drift rates ([Heathcote & Love,
  2012](https://doi.org/10.3389/fpsyg.2012.00292)), the model the LNR
  was introduced as a limit of: with a LogNormal *distance* as well as a
  LogNormal rate the two fold into one `sigma`, which is why the LNR
  never had a start-point parameter, whereas a Uniform distance leaves a
  shape the rate alone cannot produce. The threshold offset is pinned at
  1 rather than `sigma` at 1 because rescaling the evidence axis shifts
  `nu` and scales the range and the threshold but leaves a LogNormal
  rate’s `sigma` untouched, so `sigma` cannot pin the scale; `sigmabias`
  is therefore read in units of the threshold offset. The density is a
  difference of two normal CDFs and the survival one more, evaluated
  from a series below a start-point range of about `1e-4 * sigma` and
  from ratios of log-CDFs in the tails, and both agree with
  one-dimensional quadrature over the start point to `1e-7`.
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  fences the flat direction at zero the way it does for
  [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)’s
  `sigmabias`, and
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  starts it at 0.5. A formula that omits `sigmabias` now estimates it,
  as `brms` does with any dpar; fits made with earlier versions cannot
  be post-processed, because their family carries no `sigmabias`, and
  have to be refit. In
  [`rcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  and
  [`dcogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md)
  the argument sits between `ndt` (or `response`) and `poutlier`, where
  the DDM keeps its between-trial variabilities, so a call that passed
  `poutlier` by position needs it named. The vignette LNR model pins it
  at zero and was refit.

- **[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  gains `warmstart`**, which re-centres the priors on a previous fit:
  every population-level intercept and coefficient, every group-level
  SD, and every dpar left out of the formula gets
  `normal(median, prior_scale * sd)` from that fit’s posterior, matched
  parameter by parameter on the `class`, `dpar`, `coef` and `group` a
  [`get_prior()`](https://paulbuerkner.com/brms/reference/default_prior.html)
  row carries. No transformation is involved - the parameter a prior row
  is about is the parameter the source sampled, the `Intercept` prior
  being stated on the centred intercept in both - so the scales line up
  by construction.
  [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
  extracts the posterior median and SD alongside the means for this,
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) carries
  them in two new columns, and a table written before they existed still
  reads (the priors are then left alone, with a message). The
  correlations keep their LKJ and the standardized effects have no
  stated prior to change. **Unlike the rest of the warm start, this
  changes the posterior**, and if the source was fitted to data the new
  model also contains it double-counts it - a pilot on half the
  participants used to centre the priors for the fit on all of them uses
  that half twice. `prior_scale`, 3 by default, is what stands between a
  prior that only says roughly where the parameter lives and one that is
  the source’s posterior outright;
  [`?cogmod_priors`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  says when the argument is and is not legitimate.

### Bug fixes

- **[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  no longer freezes one chain in four on a cold start.**
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  started both drifts at 3; the error accumulator now starts at 1. The
  failure looked like a stuck chain - every transition at the maximum
  treedepth, step size a thousand times smaller than the other chains’,
  Rhat 1.5 to 2.9 - and it was traced to the very first warmup
  transition. A Wald density is thin on the fast side and flat on the
  slow side, so on data whose error drift is about 0.2 (the
  speed-accuracy data of `vignette("performance")`) the old start sat
  400 log-density units above the posterior. The first trajectory
  converted that into momentum along the flat `driftone` direction (the
  plateau
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  fences, where the likelihood no longer changes) and carried the chain
  from a link value of +3 to -20 in eight leapfrog steps; the step size
  then collapsed to 1e-5 within a dozen iterations and the metric
  windows that followed were estimated from a chain that no longer
  moved. Stan’s model methods found nothing numerical at the frozen
  position - the log-density and its gradient are smooth across the
  driftless branch of the survival function, and finite differences
  agree with autodiff - so the fix is where the chain starts, not what
  it computes. On an 800-trial mixed model with a 200-iteration warmup
  the old start froze a chain in 1 run in 4 with Stan’s default metric
  and in 5 of 6 when handed a metric adapted to the bulk (4 of 4 when
  the benchmark cell was rerun); the new start has done so in 0 of 6
  under the latter, the harder case, and each such fit ran in 3 minutes
  instead of 30. When the error accumulator really is as fast as the
  correct one the start is off by a factor of three on the cheap side,
  which costs a few dozen units and changes nothing.

- **[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  now starts `ndt` at half the first percentile of the observed response
  times instead of a fixed 0.1 s.** Same mechanism as the previous item,
  other cold start. The fixed value was a third of the prior median and
  safely below any ordinary data, which was the whole argument for it;
  on data whose non-decision time is 0.6 s it sat half a second low,
  every decision time looked far too long, a driftless race then fit
  better than a fast one, and the first trajectory of a cold
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  chain threw both drifts onto their flat regions - three chains in four
  on the benchmark’s `shifted` target when handed a metric adapted to
  the bulk. Half the first percentile is still below essentially every
  response, so the gradient the small start was protecting is intact,
  and it follows the scale of the data. Every `ndt` + `poutlier` family
  gets it. Starting values do not change a posterior, so no fitted model
  needs revisiting for this.

### Breaking changes

- **The default `ndt` prior is wider: `normal(-1.2, 0.5)` on the log
  scale, `lognormal(-1.2, 0.5)` for an omitted `ndt`, in place of the
  `0.2` SD.** The centre is unchanged at 0.30 s; 95% of the mass now
  sits between about 0.11 and 0.80 s instead of 0.20 to 0.44 s. The old
  SD put a 0.6 s non-decision time - not unusual for older participants
  or more demanding responses - 3.5 SDs from the centre, and on the
  warm-start benchmark’s `shifted` target (`vignette("performance")`),
  whose non-decision time is about 0.6 s,
  [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  produced divergent transitions in every run with that prior, the
  600-warmup reference included; the wider prior removed them (Rhat 1.06
  and a minimum ESS of 67 became 1.01 and 392 on the same seed). The
  prior’s job is to fence the `ndt -> 0` direction, where the likelihood
  goes flat and a flat prior would make the posterior improper; at 0.5
  it still does (0.01 s is 6.8 SDs out) without telling the data where
  in 0.1 to 0.8 s the non-decision time is. Fits with the default priors
  will move slightly, most where the data put `ndt` far from 0.30 s,
  which is where they should have been free to move. Every `ndt` +
  `poutlier` family is affected.

## cogmod 0.3.1

### New features

- **New
  [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)**
  turns what a previous fit’s warmup produced - the adapted inverse
  metric, the step size and the posterior means - into the `init`,
  `inv_metric` and `step_size` arguments of a new
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call, so
  that a refit, or the same model on more participants, can run a much
  shorter warmup. Stan adapts one variance per unconstrained parameter,
  and the function labels each with its Stan name (read off the
  generated program with the parser
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  already has) and joins the two models on those names: population-level
  entries carry over, a pilot participant’s standardized effects follow
  it by level name to its position in the bigger model, new participants
  take their effect’s average variance and start at zero, and anything
  without a counterpart gets Stan’s default variance and a generic
  start, with a count in [`print()`](https://rdrr.io/r/base/print.html).
  The standardized effects and Cholesky factors that `brms` drops from a
  saved fit are rebuilt from the `r_`, `sd_` and `cor_` it keeps.
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives a
  table of a few kilobytes that survives
  [`write.csv()`](https://rdrr.io/r/utils/write.table.html) and can be
  passed back as a file path, so a pilot fitted on a laptop can
  warm-start an array job on a cluster. On the other side,
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  gains a `warmstart` argument and the new
  [`cogmod_inv_metric()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
  and
  [`cogmod_step_size()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
  share its signature - the model’s formula and data, then the source,
  be it a fit, a table or a file - so that each argument of
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) has one
  helper and the table is mapped onto the model on the way. Whatever of
  `formula` and `data` is not given is taken from the source fit. Works
  for any `brms` model fitted with the `cmdstanr` backend and the
  diagonal metric. On a mixed LNR and a mixed DDM, a pilot on 4 of 8
  participants warm-started the full fit to about twice the effective
  draws per second of a cold start with the full warmup, and four to six
  times those of a cold start with the same short warmup; the starting
  values alone bought nothing, so the metric and step size are the
  product (`vignette("performance")`).

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
  this is the *simple* censored shifted Wald of [Miller et
  al. (2018)](https://doi.org/10.1177/0146621617710465), their Eq. 4,
  the `version = "simple"` of the `cswald` model in
  [`bmm`](https://github.com/popov-lab/bmm). It is not their
  competing-risks variant (Eq. 5, a race of two Wald accumulators with
  drifts `v` and `-v`, as implemented in `rtdists` and `bmm`’s
  `version = "crisk"`), which is a choice model rather than a censoring
  construction;
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  with `bias` fixed at 0.5 covers that ground. Here censoring is not a
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
  warns above 20% censored trials - a threshold that is exact for
  timeouts and omissions at any rate but lenient for commission errors,
  where the construction is biased well before it (see
  [`?rcogmod_invgaussian`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md)).
  What the model is for, what it assumes and the one check to run before
  using it are in
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

### Bug fixes

- **[`dcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  and everything built on it
  ([`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html),
  [`loo()`](https://mc-stan.org/loo/reference/loo.html),
  [`p_outlier()`](https://dominiquemakowski.github.io/cogmod/reference/p_outlier.md),
  the other R-side post-processing of
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  fits) are now accurate when the `sigmandt` range reaches down to fast
  decision times.** The R density integrates the non-decision time out
  with a fixed 25-node Gauss-Legendre rule, and when that range covers
  decision times from about zero up to the response - a fast response,
  or a wide `sigmandt` - the integrand holds the whole early peak of the
  first-passage density inside a sliver of it, which 25 nodes on the
  plain time scale cannot resolve. It is the defect reported against
  [`rtdists::ddiffusion()`](https://rdrr.io/pkg/rtdists/man/Diffusion.html)
  in [rtdists issue 28](https://github.com/rtdists/rtdists/issues/28),
  and it was here too: on the issue’s own example (`boundary = 0.5`,
  `drift = 0.5`, `bias = 0.3`, `sigmandt = 0.16`, decision time 0.16)
  the density was out by 5e-4, by 3% with `bias = 0.1` and
  `sigmandt = 0.2`, and by 50% with the start point almost on the
  responding boundary. The rule now runs over *log* decision time, from
  the point where the density is dead rather than from zero, so the peak
  is about one log unit wide wherever it sits and the same 25 nodes
  resolve it at any time scale, and the rule takes more nodes only when
  the log range is wide enough to need them: the worst error over the
  issue’s sweep is now below 1e-9, and below 4e-12 over a much broader
  grid, against a converged 1600-node rule. Densities with
  `sigmandt = 0` are unchanged to the last bit, and the common case -
  `sigmandt` well inside the response time, where the old rule was
  already accurate - gives the same values at the same cost. The Stan
  likelihood uses Stan’s own adaptive `wiener_lpdf()` and was never
  affected, so fitted models are unchanged; only their R-side
  post-processing moves, and only in that regime. The R-versus-Stan
  density test is tightened from a relative 1e-4 to 1e-5 accordingly,
  which is Stan’s own tolerance.

- **[`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md)
  no longer produces divergent transitions by the hundred on healthy
  posteriors.** The Stan log-survival of the losing accumulator formed
  its reflection term as `log_diff_exp(log R(b), log R(k))`, and for a
  response less than about half a millisecond above the non-decision
  time both normal CDFs in `R` round to exactly 1, so it evaluated
  `log_diff_exp(0, 0)`. The *value* is fine (`-Inf` for a term that
  really is negligible there, which is why the R-versus-Stan density
  tests never caught it), but its reverse-mode adjoint is `0 / 0`, and
  that `NaN` propagated into the gradient of every parameter. Stan
  reports a `NaN` gradient as a divergent transition, and because `ndt`
  is estimated a few milliseconds below the fastest responses, most
  trajectories crossed one of those windows: on the lexical decision
  data of the decision-making article, 900 trials from 6 participants,
  between a third and two thirds of the transitions were divergent -
  with population-level effects only or with participant intercepts,
  under `diag_e` or `dense_e` - while `Rhat` and the effective sample
  sizes said the posterior was fine, because it was. The difference is
  now assembled from three pieces that each stay away from the saturated
  end of the normal CDF, at the cost of two extra normal CDFs on early
  responses only. Values are unchanged to `1e-12` on the log scale; the
  same fits now run without a divergence (population-level) or with the
  handful the other families also show under `dense_e` with random
  effects. The other race families are unaffected. A gradient regression
  test guards it, gated behind `COGMOD_TEST_SLOW` like the other tests
  that compile a model of their own.

- **[`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  no longer returns the odd 10-40 s response in place of a fast one.**
  The sampler inverts the CDF with Newton’s method, and took a step
  under `1e-10` in log time as convergence. A pass that lands far in the
  tail finds the survival and the density both denormal, and their ratio
  makes the step look tiny while the residual is still hundreds of log
  units off, so the draw was accepted where it stood. It hit about one
  draw in 4,000 to 20,000 at short boundaries or strong drifts
  (`boundary = 0.3`, `bias = 0.3`, `drift = -5` is one such cell;
  `boundary = 2`, `bias = 0.7`, `drift = 6` another), which is rare in
  [`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  and a handful of absurd draws per observation in
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
  where every observation gets thousands. Convergence now also requires
  a small residual, and a test pushes every draw back through the CDF to
  check it lands on its own quantile.

### Performance

- **[`rcogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  is 3x faster, and
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  on a
  [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
  model up to 16x.** The sampler inverts a series whose length is set by
  the fastest response it could be asked for, and used that length for
  every draw: 41 terms at the default start point where the median draw
  needs 6, 205 at a start point of 0.1. It now runs in stages - 16 terms
  settle the bulk of the draws, and only the responses too fast for that
  many go round again with four times as many, until the full series is
  reached. Nothing is approximated: a draw is only accepted from a stage
  whose series is exact at its root, and the draws agree with the CDF to
  `1e-12`. The gain is largest where the parameters vary across draws,
  as they do in
  [`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html),
  because one extreme posterior draw used to set the series length for
  all of them. Converged draws now also drop out of the Newton
  iteration, and a fast response starts from the single-barrier
  small-time approximation rather than from the floor of the bracket.

- **The choice families’ `posterior_predict_*()` methods take a vector
  of observations**, returning the draws stacked with those for `i[1]`
  first.
  [`brms::posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
  calls the method once per observation, and with a few dozen draws per
  call about half of each call is fixed cost and the loop adds as much
  again; predicting in chunks of ~50 observations from a prepared
  `brmsprep` instead runs a posterior predictive check on 2,500 DDM
  trials in about a third of the time. The recipe is in
  [`?posterior_predict_cogmod_ddm`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md).
  The DDM sampler’s own fixed cost per call is also down by about 15%,
  from indexed assignment in place of
  [`ifelse()`](https://rdrr.io/r/base/ifelse.html) and no column copies
  while every draw is still active; the draws are bit-identical.

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

### Documentation

- The performance article is reorganised from the suggestions with no
  downside to the ones that need judgement, and gains four sections.
  **Compiler optimizations**: stanc’s `O1` and CmdStan’s
  `STAN_CPP_OPTIMS` and `STAN_NO_RANGE_CHECKS`, passed through
  `stan_model_args`, and what each one does. **Mass matrix adaptation**
  (`metric = "dense_e"`): why the boundary/ndt and drift/boundary
  trade-offs of evidence accumulation models make the default diagonal
  metric a poor fit, what the dense metric costs as the number of
  parameters grows, and how to pass it through either backend. **Warm
  starts**: reusing the adapted metric and step size that `brms` keeps
  in a fit’s metadata to shorten the warmup of a refit; carrying a pilot
  fit’s metric, step size and posterior means over to the same model on
  more participants, which has more parameters, by mapping the metric
  across by parameter name (about twice the effective draws per second
  of a cold start on a mixed LNR and a mixed DDM, where the pilot’s
  initial values alone bought nothing); and a `cmdstanr`-level pipeline
  that initializes MCMC from Pathfinder draws and their unconstrained
  covariance, then wraps the result back into a `brmsfit`, with the
  reasons never to fix the metric to a variational approximation. The
  approximation section now also covers the **Laplace approximation**
  (`algorithm = "laplace"`) and how it compares with Pathfinder. Each
  section reports what the option bought on the DDM, LBA, LNR and RDM in
  a local benchmark; the scripts behind those numbers live in
  `benchmarks/` (not part of the installed package) and can be rerun on
  any model.

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

CRAN release: 2026-09-12

- CRAN Publication.
