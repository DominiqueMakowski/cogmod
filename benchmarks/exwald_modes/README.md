# The ex-Wald's modes on the RT article's data (2026-10-12)

`benchmarks/driftone_priors/` found that `cogmod_exwald()` can split its
chains on the RT article's data. Seed 1 left one chain about 11 log-density
units below the other three, in both prior arms (R-hat 1.21-1.23, minimum
bulk ESS 3-4). This is the follow-up: reproduce it, say what the low chain is,
and find whether the starting values or the likelihood put it there.

The model is the article's own:

```r
bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition, tau ~ Condition,
   family = cogmod_exwald())
```

fitted to `readRDS("vignettes/models/m_wald4.rds")$data`: three participants of
rtdists' `speed_acc`, correct responses only, 2,180 Accuracy and 2,105 Speed
trials. 4 chains, 250 warmup and 250 sampling iterations, `cogmod_inits()`
after `set.seed(seed)`.

**Outcome (2026-10-12).** Documented rather than changed:

- `?rcogmod_exwald` has a section on the two solutions.
- The RT article's ExWald section carries a short version.
- The priors and inits are as they were.
- The drift prior below (P1) was measured but not adopted.

The second mode is a property of the model, not of this data set's size
("More data?"), so no default makes it go away.

## Summary

- **It reproduces.** On this machine, 0.3.4 and seed 1 left chain 3 at a mean
  `lp__` of 2358.0 against 2368.6-2369.2 (R-hat 1.54, minimum bulk ESS 7.5).
- **The posterior has four modes, not two.** Each condition has two solutions.
  In W the Wald's own skew makes the right tail; in E a near-Gaussian Wald and a
  long exponential stage do (an ex-Gaussian with a shift). Every
  combination of the two cells is a local mode: WW, WE, EW, EE (Accuracy first).
- **The chains sample the wrong one.** From `cogmod_inits()` they find WW, or
  EW: that is the low chain. **None of 28 chains found WE or EE**, which is
  where the posterior's mass is. Bridge sampling puts it at:

  | priors | WE | EE | WW | EW |
  | --- | --- | --- | --- | --- |
  | 0.3.4 | 80% | 19% | 0.04% | 3e-8 |
  | 0.3.5 | 36% | 63% | 1.1% | 3e-7 |

  The cached vignette fit (`m_exwald.rds`) is WW. Its four chains agree, so
  nothing flags it.
- **They are real modes, not a slow-mixing region.** Every path from WW to EW
  dips at least 22.8 units below WW and 11.8 below EW. From WW to WE the dip is
  21.0 and 26.9. A chain that lands in one stays there.
- **More data would not remove them.** At ~2,000 trials a condition, data
  simulated from the ex-Wald itself do not reliably reveal which solution is
  true. With 8,000-32,000 they mostly do. But the second mode was present in
  all 36 simulated data sets from 2,000 trials up, and the wall between the
  modes grows with the data (see "More data?").
- **The density is not the cause.** The valley between W and E lies on the
  seam `mu^2 = 2 / tau`, where `cogmod_exwald_decision_lpdf()` switches
  branches. But both branches match brute-force convolution to 2e-10 in log
  density. The log-likelihood and Stan's gradient pass through the seam without
  a step, and Stan's `lp__` equals the R twin's to 1e-11.
- **The starting values pick the mode; the jitter is not what decides.** The
  same `cogmod_inits()` draw went to EW on Artemis and to WW here. In the three
  runs of seed 1, a different chain got stuck each time. Chain 3 here was
  captured at warmup iteration 15, by one 1,023-step trajectory, while still 48
  units below the typical set.
- **W and E predict the same response times.** From the 1st to the 95th
  percentile they agree within 8 ms in Accuracy and 3 ms in Speed. They part
  in the slow tail and in the parameters. In WE, Speed instructions raise the
  drift from 2.5 to 9.6 and the boundary from 0.70 to 1.89. In WW they raise the
  drift from 2.5 to 3.6 and leave the boundary at 0.74. The article's
  Speed - Accuracy figure and its LOO table both read off WW.
- **A drift prior moves the mass back to WW; the tau prior does not.** P1 is
  `normal(3, 1)` on `mu`'s intercept and `normal(0, 1)` on its slope, both on
  the softplus link, on top of 0.3.5:
  - it puts 99.9% of the mass in WW by bridge sampling;
  - it moves the WW estimates by under 0.02 and costs nothing per draw;
  - it leaves the E modes in place, so a chain can still be captured (1 of 8
    from `cogmod_inits()`);
  - its slope half breaks `cogmod_priors()`'s rule that the response's slopes
    are left to brms.

  **Not shipped**: see "A fix" for what it would take.

## The four modes

Modes of `lp__` (Stan's unconstrained scale, Jacobian included), found by
optimising from the 16 `cogmod_inits()` starts of seeds 1-4 and one
hand-placed start per combination (`modes.R`). Natural parameters per
condition (`mu`, `boundary`, `tau`, `ndt`); `results/modes_*.csv` has the rest,
including `k2 = mu^2 - 2 / tau`, which says which branch of the density a cell
is on.

| 0.3.4 | lp | log-lik | Accuracy | Speed | mass: Laplace / bridge |
| --- | --- | --- | --- | --- | --- |
| WE | 2379.47 | 2388.96 | 2.47, 0.70, 0.071, 0.336 | 9.64, 1.89, 0.122, 0.231 | 86% / 80% |
| EE | 2374.73 | 2383.65 | 7.16, 1.53, 0.209, 0.266 | 10.23, 2.17, 0.123, 0.215 | 14% / 19% |
| WW | 2373.57 | 2381.46 | 2.49, 0.73, 0.062, 0.333 | 3.59, 0.74, 0.037, 0.306 | 0.04% / 0.04% |
| EW | 2362.58 | 2374.74 | 6.33, 1.17, 0.207, 0.297 | 3.57, 0.71, 0.042, 0.309 | 4e-8 / 3e-8 |

- **Which cell prefers which solution.** Accuracy prefers W by 6.1 log-lik
  units, Speed prefers E by 7.1. The priors widen both gaps a little.
- **The E modes are ridges.** The smallest Hessian eigenvalue of `-lp__` is
  0.5-1.1 there, against 6.1 at WW. In an E cell the shift `ndt` trades against
  the mean of the near-Gaussian Wald, and only the `ndt` priors hold it:
  Speed's likelihood alone runs E down to `ndt` -> 0. HMC chains in E wander
  along the ridge, down to `ndt` 0.19-0.23 in Speed.
- **The ridge's width is part of why E carries the mass.** Its volume adds to
  its height. Bridge sampling agrees with the Laplace log-mass to 0.03 for WW
  and EW and to 0.4 for EE, the widest ridge.
- **Sampling an E mode costs more.** Per post-warmup draw it took 88 leapfrog
  steps in EW, 91-97 in WE and 237-267 in EE, against 41-46 in WW. Step sizes
  fall from 0.09-0.11 to 0.017.

0.3.5 adds a `boundary` prior, `normal(0.5, 1)` with slope `normal(0, 0.5)`.
WE needs a large boundary slope, so the prior pulls it down to 2375.61. EE
(2373.32) and WW (2373.04) are left nearly level. EW is at 2361.41.

HMC confirms the heights. Chains started inside each mode (`fit.R ...
modes:`) stayed there for all 250 sampling iterations, with these mean
`lp__`:

| priors | WW | WE | EE | EW |
| --- | --- | --- | --- | --- |
| 0.3.4 | 2368.6-2369.2 | 2374.7, 2375.3 | 2369.9, 2370.2 | 2358.0 |
| 0.3.5 | 2368.2 | 2371.0 | 2368.7 | 2357.1 |

## The density is not the cause

Each W/E valley sits on the seam `mu^2 = 2 / tau` (`barrier.R`):

- **Accuracy.** The profile of `lp__` along its `tau` bottoms out at 2350.74,
  at `tau` 0.170, where `k2` = -2.0.
- **Speed.** It bottoms out at 2352.59, crossing the seam between 0.090 and
  0.095.

The seam is where the density switches from the Faddeeva branch to the closed
form, so the switch was the first suspect. `numerics.R` checks it on the 0.3.5
program:

- **Values.** `.dexwald_raw()` against brute-force convolution (adaptive
  quadrature split at the Wald mode, rel.tol 1e-12) differs by at most 1.8e-10
  in log density on the Faddeeva branch, and 2e-15 on the closed form. That
  holds at every cell of every mode and along the straight line from WW to EW.
- **The seam itself.** Step the Accuracy drift through `sqrt(2 / tau)` by
  relative offsets of +-1e-2 down to +-1e-10. The log-likelihood changes only
  as its slope predicts, by 2e-8 between -1e-10 and +1e-10, and Stan's
  `d lp / d Intercept` by 1.2e-7: -30.78345666 against -30.78345678.
- **Stan.** `log_prob` equals the R twin to 1.4e-11 at every point. Its
  gradient agrees with central differences to 1.03e-5 relative or better, the
  same on both sides of the seam, and is finite everywhere. At the four modes
  the gradient norm is 1e-5 to 1e-4.

The valley is on the seam for a reason in the model, not the code (this part
is interpretation). The right tail decays at `min(mu^2 / 2, 1 / tau)`, the
slower of the two stages. Either side of the seam one stage owns the tail and
the other is free to shape the bulk. On the seam both decay at the same rate,
so both stages are wide at once and nothing fits the bulk and the tail
together.

## Why the chains land where they do

`cogmod_inits()` starts every chain far below the typical set: the first
warmup draws sat at `lp__` 545-1,698 against 2,369. It starts `ndt` at half the
first percentile of the RTs (0.14-0.20 s with jitter, against a posterior
0.31-0.33), and `tau` at 0.15 s, a hair below the Accuracy valley at 0.17.

In the first 10-20 warmup iterations the chain descends with the unit metric,
and it can still cross between the basins until its `lp__` passes the saddle's,
about 2,351. That is where it gets caught.

Local chain 3 took this route:

1. At iteration 14 it sat at `tau_A` 0.159 and `lp__` 2,321.
2. At iteration 15 one treedepth-10 trajectory (step size 0.0016) took its
   Accuracy drift from 3.2 to 7.2.
3. It never came back.

Chain 1 started on the E side of the seam (`k2_A` = +3.5) and still ended in WW,
by iteration 22.

What decides the outcome is the descent's dynamics, not the draw:

- **Cross-platform.** `cogmod_inits()`' RNG is platform-independent, so all
  three runs of seed 1 started from identical points. (cmdstanr calls the init
  function once before the per-chain calls, so chain k gets the (k + 1)th
  draw.) Linux and Windows floating point then sent different chains into EW:
  chain 2 on Artemis under 0.3.4, chain 4 under 0.3.5, chain 3 here.
- **Overall rate.** Across 28 chains from `cogmod_inits()` under the shipped
  priors (seeds 1-4, both machines), EW caught 3, WE/EE none, WW the rest.
- **How the chains were labelled.** The 12 local chains are labelled from
  their draws. For the 16 Artemis chains only per-chain `lp__` and pooled
  summaries survive:
  - None sits at WE's height, 6 units above WW.
  - EE's mean `lp__` is close to WW's, so that alone cannot rule it out. But
    on seed 2 every parameter has R-hat 1.00 and per-chain medians within 0.34
    posterior SD of each other (`chain_spread`), which one E chain among four
    would break.
  - The low seed-1 chains sit at EW's `lp__` and pull `b[1]`'s 5% quantile
    to -3.5. EW's drift slope is -2.8 on the link, against WW's +1.2.
- **Where the E modes are.** The Speed E mode sits at drift 8-10 and boundary
  1.4-2.2, far beyond the start's 3 and 0.5 on a ridge the descent never
  climbs.

Inits could steer chains to one mode or another. They cannot make a
multimodal posterior unimodal. Steering them to WW under the current priors
would hide that WW is the minor mode.

## What W and E disagree about

Predicted RT quantiles in seconds at the 0.3.4 modes, against observed
(`shape.R`). Speed is compared WW against WE, Accuracy WW against EW.

| | 1% | 5% | 25% | 50% | 75% | 95% | 99% | 99.9% |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| Accuracy observed | 0.412 | 0.450 | 0.539 | 0.635 | 0.765 | 1.145 | 1.620 | 1.961 |
| W | 0.413 | 0.449 | 0.535 | 0.632 | 0.778 | 1.124 | 1.500 | 2.087 |
| E | 0.412 | 0.450 | 0.536 | 0.634 | 0.781 | 1.116 | 1.448 | 1.918 |
| Speed observed | 0.367 | 0.401 | 0.460 | 0.518 | 0.597 | 0.784 | 1.096 | 1.636 |
| W | 0.373 | 0.399 | 0.457 | 0.519 | 0.608 | 0.802 | 1.000 | 1.299 |
| E | 0.372 | 0.401 | 0.459 | 0.518 | 0.606 | 0.803 | 0.998 | 1.274 |

The log-likelihood difference is spread through the distribution rather than
concentrated in the tail (`results/shape_bins.csv`). In Speed, E gains 23 units
on the 25-50% bin and 7 on the fastest 1%, several of which fall below W's
`ndt` of 0.306 s. It loses 19 on the 1-10% bins and 8 on 75-90%.

The mimicry is close enough that only the parameters tell the two apart. In
the E solutions `mu` and `boundary` stop meaning drift and threshold and only
set a near-Gaussian's mean and spread. A Speed condition that raises the
threshold from 0.70 to 1.89 is the consequence.

## More data?

Is this a shortage of data - too many parameters for 2,000 trials a
condition - rather than a property of the model? `data_size.R` fits one
condition at a time (the four parameters plus `poutlier`, 0.3.5's priors) to:

- **Simulated data.** Drawn from the ex-Wald itself at each of the four cell
  solutions above, 500 to 32,000 trials, three replicates.
- **The real trials.** Subsampled to 25% and 50%.

For each it optimises from a W start and an E start. W - E is the lp gap
between the two optima, positive when W is higher (mean over replicates). The
barrier is from the first replicate, measured from W / from E:

| truth | trials | two optima | W - E | barrier |
| --- | --- | --- | --- | --- |
| Accuracy W | 500 / 2,000 / 8,000 / 32,000 | 1, 3, 3, 3 of 3 | 1.2 / 0.5 / 6.4 / 22.2 | - / 2.1, 1.4 / 25.7, 22.4 / 116, 99 |
| Accuracy E | same | 3, 3, 3, 3 | -0.4 / -0.8 / -4.9 / -16.3 | 0.5, 0.5 / 3.1, 3.6 / 30.3, 34.9 / 116, 131 |
| Speed W | same | 3, 3, 3, 3 | -0.1 / -1.5 / 0.5 / 1.5 | 2.6, 2.0 / 6.4, 7.2 / 8.6, 8.9 / 44, 41 |
| Speed E | same | 1, 3, 3, 3 | -0.1 / -3.6 / -9.2 / -37.3 | - / 1.2, 3.0 / 19.4, 28.9 / 98, 144 |
| real Accuracy | 545 / 1,090 / 2,180 | 3, 3, 1 of 3, 3, 1 | 2.2 / 2.9 / 6.3 | 4.3, 2.1 / 8.2, 4.8 / 18.2, 11.9 |
| real Speed | 526 / 1,052 / 2,105 | same | -0.8 / -3.7 / -6.3 | 1.5, 2.8 / 7.2, 8.9 / 20.1, 26.4 |

The hypothesis is half right:

- **Too few trials to choose.** At this data set's ~2,000 trials a
  condition, data generated by the ex-Wald itself do not reliably say which
  solution is the true one. The gaps are 1-5 units either way, and under the
  Speed W truth the wrong solution scores higher in all three replicates. With
  more data the truth wins: by 6-37 units at 8,000-32,000 trials in three of
  the four cells.
- **Speed W is the exception.** Its E mimic is still within 1.5 units at
  32,000 trials. For some parameter regions no realistic amount of data
  separates them.
- **More data never removes the second mode.** From 2,000 trials up it was
  present in all 36 simulated data sets. At 500 it was missing in 4 of 12: the
  two solutions blur together when data are scarce, not when they are
  plentiful.
- **The wall between the modes grows with the data,** roughly in proportion:
  a few units at 2,000 trials, 20-35 at 8,000, 40-145 at 32,000. The real
  data show the same, from about 2-4 at a quarter of the trials to 12-26 at
  all of them.

So more data sharpens which mode holds the posterior mass. But a chain that
settles in the other mode is held there more firmly, and a cold start's
descent still has to choose. The two solutions are a property of the
ex-Wald's likelihood - two ways of making a right tail - not of having too
few trials.

## A fix

Ruled out:

- **The density.** It is exact; see above.
- **Inits alone.** See above.
- **The tau prior.** It is copied from the ex-Gaussian, centred at 0.2 s,
  which is where E puts `tau` - that looked like the culprit. Moving its
  centre down to 0.079 s and 0.049 s raises WW's share from 1.1% only to 7%
  and 16%. The E modes keep a `tau` of 0.12-0.21 s regardless, and the mass
  moves from EE to WE.

What separates W from E is the drift. Under 0.3.5, `mu` keeps brms's
`student_t(3, 0, 2.5)` intercept and a flat slope. Shares of the mass by
`priors.R`, which re-optimises 0.3.5's four modes under each variant (Laplace;
EW is under 0.001% in every row):

| `mu` prior, softplus link | WW | WE | EE |
| --- | --- | --- | --- |
| 0.3.5 (student_t intercept, flat slope) | 1.1% | 39% | 60% |
| slope normal(0, 2) | 2.6% | 1.6% | 96% |
| slope normal(0, 1) | 4.9% | 0.0004% | 95% |
| intercept normal(3, 1) | 13% | 87% | 0.04% |
| intercept normal(3, 1.5) + slope normal(0, 2) | 46% | 26% | 28% |
| intercept normal(3, 1) + slope normal(0, 2) | 86% | 14% | 0.2% |
| **P1: intercept normal(3, 1) + slope normal(0, 1)** | **99.9%** | 0.004% | 0.1% |

The intercept stops EE, whose two cells both sit at drifts of 7-10. The slope
stops WE, whose Speed drift is 6-7 link units above Accuracy's. Neither alone
is enough.

P1 measured with HMC:

- **Bridge sampling.** Over chains started inside each mode: WW 99.9%, EE
  0.12%, WE 0.004%, EW 3e-6.
- **The E modes stay, lower.** Chains started in them stayed for all 250
  draws. Mean `lp__` was 2358.1 in WE, 2359.5 in EE and 2355.2 in EW, against
  2368.6 in WW.
- **The E ridges are bounded.** EE needed 89 leapfrog steps per draw, against
  237-267 under 0.3.4.
- **The WW estimates barely move.** Drift 2.494 -> 2.505 (Accuracy) and
  3.596 -> 3.583 (Speed), boundary 0.734 -> 0.740 and 0.744 -> 0.739.
- **From `cogmod_inits()`, seeds 1-2.** 7 of 8 chains went to WW. Seed 2's
  chain 3 was captured into EW at warmup iteration 20 and stayed there, 13
  units below. P1 does nothing about the capture.
- **No cost in WW.** P1's WW chains took 9.7-11.4k leapfrog steps per 250
  draws, against 10.2-11.5k under 0.3.4. Bulk ESS per leapfrog step is
  1.00-1.20 times 0.3.4's (median 1.08) over the ten quantities, R-hat
  <= 1.01. That compares the WW chains of each seed-1 run, three against
  four.

Open, and why it is not shipped:

- **The slope half breaks a rule.** `cogmod_priors()` deliberately never puts
  a prior on the response's slopes ("they are the effects the model is being
  fitted to estimate", R/cogmod_priors.R). P1 needs one: the intercept alone
  leaves WE at 87%. That is a design decision, not a measurement.
- **It is tuned on one data set.** `normal(3, 1)` puts 95% of the drift
  between 1.3 and 5.0. A fast task with a short decision stage would legitimately
  sit higher, and the prior would push it towards a longer `tau` or a later
  `ndt`. Before shipping, it needs checking on a second data set; the article's
  DoggoNogo simple-RT data is the obvious one.
- **It does not remove the E modes.** A chain can still be captured in EW, as
  now, and R-hat stays the guard. What P1 changes is that the mode the chains
  find is then also the mode with the mass.
- **The usual checks are not done.** The gradient check is unaffected: priors
  are brms's own `normal_lpdf()`. Gradient cost should be too (two scalar
  terms), but it has not been timed with `gradient_cost.R`.
- **The article needs a refit either way.** Under the priors it was fitted
  with, the posterior is not where `m_exwald.rds` is. Had its chains found WE,
  ExWald's LOO entry would move by roughly the 7.5 units by which WE's
  log-likelihood exceeds WW's, and its Speed - Accuracy panel
  would show a boundary effect of +1.7 on the link instead of 0.

## Files

| file | what |
| --- | --- |
| `common.R` | data, formula, the R twin of `lp__` for 0.3.4 / 0.3.5 / P1, the program for each, W/E labels |
| `modes.R` | the modes of a prior set: multi-start optimisation, Hessian, Laplace |
| `barrier.R` | `lp__` profiled along one cell's `tau` between two modes |
| `numerics.R` | brute-force convolution, the seam, Stan against the R twin |
| `shape.R` | where W and E differ: log-lik by RT bin, predicted quantiles |
| `priors.R` | 0.3.5's four modes under drift and tau prior variants |
| `data_size.R` | one condition, simulated (500-32,000 trials) and subsampled: do both optima survive more data? |
| `fit.R` | one HMC run, from `cogmod_inits()` seeds or from the modes |
| `collect.R` | per-chain and per-run summaries of `fit.R` runs |
| `bridge.R` | bridge-sampled mass of each mode from `fit.R` draws |
| `results/` | the CSVs behind every number above (`chains.csv` has every local chain) |

From the package root (`EXW_SCRATCH` holds the draws):

```
Rscript benchmarks/exwald_modes/modes.R 0.3.4      # also 0.3.5, P1
Rscript benchmarks/exwald_modes/barrier.R 0.3.4
Rscript benchmarks/exwald_modes/shape.R 0.3.4
Rscript benchmarks/exwald_modes/priors.R
Rscript benchmarks/exwald_modes/numerics.R
Rscript benchmarks/exwald_modes/fit.R 0.3.4 seed:1 034_seed1
Rscript benchmarks/exwald_modes/fit.R 0.3.4 modes:benchmarks/exwald_modes/results/modes_034.rds 034_modes
Rscript benchmarks/exwald_modes/collect.R <scratch>/034_seed1 <scratch>/034_modes ...
Rscript benchmarks/exwald_modes/bridge.R 0.3.4 <scratch>/034_seed1 <scratch>/034_modes
```

A 4-chain fit takes 30-50 minutes on an 8-core laptop, the E-mode chains being
the slow ones. The 0.3.4 runs used here started their E chains from WE and EE
only, two each, and the seed 3-4 runs were short (150 warmup, 10 sampling
iterations), long enough to see where each chain settled.
