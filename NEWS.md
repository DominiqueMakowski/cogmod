# cogmod 0.2.0

## New features

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

## Bug fixes

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

  `cogmod_rdm()`, `cogmod_lba2()` and `cogmod_ddm()` - the other choice+RT
  families - have not moved yet and still use `tau` + `minrt`.

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

First CRAN release.

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
