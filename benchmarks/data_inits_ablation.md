# Data-aware starting values: does moment-matching the start buy anything?

2026-09-22. Companion to the plan that generalised `.ndt_start()` into a
per-family layer reading every natural-scale start off the response
(`.data_start()` in `R/cogmod_inits.R`, `data_init` on the `.SHIFTED`
registry entries, `.DATA_INIT_PLAIN` for the two ex-Gaussians). This note
records what the layer was measured to do, because the answer is "nothing
a fit can see", and the plan said up front that this result should be kept
as carefully as a positive one.

## What was built

Twelve families got a rule: the log mean and log SD for the four
log-location families (each corrected for what its shape, degrees of freedom
or start-point range adds to those moments); the Wald inversion
`drift = sqrt(mean / var)`, `boundary = drift * mean` for the Wald and its
closed-form twin for the Birnbaum-Saunders; two-moment inversions for the
Gamma, inverse Gamma, Weibull and Frechet (the last two by root-finding on
the coefficient of variation); and the skewness-based estimator for the two
ex-Gaussians. `ndt` is placed first by `.ndt_start()` and every rule works
on `RT - ndt`. Every rule falls through to the old constants on fewer than
20 usable responses, a variance that is zero or not finite, or any estimate
outside the parameter's support.

Three things the plan had wrong were caught by deriving each rule against
its registry entry rather than assuming it:

- The log-Weibull's Gumbel location is not the log mean:
  `mu = meanlog - 0.5772 sigma`, `sigma = sdlog * sqrt(6) / pi`.
- The log-Student's log SD is `sigma * sqrt(dof / (dof - 2))`, a factor of
  1.29 at the start's `dof = 5`.
- The LogNormal rule is the exact MLE only at `sigmabias = 0`. The constant
  start of 0.5 adds 0.216 to the log mean and a little to its variance;
  both are taken out, using whatever value the model actually has.

The ex-Gaussian rule needed a cap on the sample skewness at 1.7.
`sigma = sqrt(var - tau^2)` is a difference of nearly equal numbers once the
tail dominates, and a sample skewness a little high starves the Gaussian
stage - the expensive direction. Uncapped, over 300 samples of 200 trials
from an ex-Gaussian of skewness 1.4, the worst start was 1.0 log-likelihood
units per observation below the truth and 5% of samples declined; on real
lexical-decision data (speed_acc, per-participant skewness 1.1-3.0, median
2.1) it declined for 11 of 17 participants. Capped at 1.7 the worst case
was 0.04, nothing declined, and every participant started 0.02-0.22 log
units per observation above the constants.

## Evidence that the starts are closer (no Stan)

Mean log-likelihood per observation, at the constant start (with `ndt` from
`.ndt_start()`) and at the data start, against the truth where known.
1000 simulated trials per cell; 2000 real speed_acc RTs (<= 2 s).

| Data | Constant start, short of truth | Data start, short of truth |
| --- | --- | --- |
| Typical simulated (12 families) | 0.04 - 0.99 | 0.01 - 0.23 |
| Slow participant simulated (ndt 0.45 s, longer decision times) | 2.2 - 5.7 | 0.01 - 0.43 |
| speed_acc, gain of data start over constant | - | +0.004 - +1.2 |

On the slow data the constants sit thousands of log-density units below the
data over a whole data set, which is the situation the RDM brittleness
report traced stuck chains to. The Wald and Birnbaum-Saunders constants were
the furthest off on real data (-0.97 vs +0.26 per observation). Adding 1%
slow lapses (2.5-4 s) to the real data shrank every gain but reversed none
by more than 0.01.

Between-chain dispersion of the starts is identical under both schemes,
because the jitter is additive on the link scale and does not depend on
where the centre is. The plan's worry that data starts would weaken R-hat
by starting all chains in one neighbourhood does not arise: the constant
only ever contributed a common offset, which R-hat never benefited from.

## Evidence that the starts help a fit (Stan)

Same data, priors, seeds and settings; the only difference is the `init`
function. 10 participants from speed_acc, 50 trials per condition (1000
rows), a participant random intercept on `mu` and a Condition effect on the
drift or location; 4 chains, 500 warmup + 500 sampling, 2 seeds. The
"constant" scheme is `cogmod_inits()` as it was before this change.

speed_acc as is (median RT 0.53 s):

| Family | Scheme | Warmup s (seed 1 / 2) | Min bulk ESS (1 / 2) | Max R-hat (1 / 2) | Divergences |
| --- | --- | --- | --- | --- | --- |
| LogNormal (`sigmabias = 0`) | constant | 72 / 65 | 269 / 327 | 1.019 / 1.019 | 0 |
| LogNormal | data | 66 / 60 | 390 / 348 | 1.008 / 1.011 | 0 |
| Wald (`sigmadrift = sigmandt = 0`) | constant | 95 / 92 | 440 / 432 | 1.008 / 1.006 | 0 |
| Wald | data | 97 / 93 | 486 / 506 | 1.007 / 1.010 | 0 |
| Ex-Gaussian | constant | 107 / 107 | 349 / 461 | 1.014 / 1.014 | 1 |
| Ex-Gaussian | data | 108 / 100 | 435 / 420 | 1.007 / 1.012 | 0 |

speed_acc shifted to `0.3 + 1.1 * RT`, the warmstart ablation's slow
population, where the Wald's constant threshold of 0.5 faces a data start
of 2.0:

| Family | Scheme | Warmup s (seed 1 / 2) | Min bulk ESS (1 / 2) | Max R-hat (1 / 2) | Divergences |
| --- | --- | --- | --- | --- | --- |
| LogNormal | constant | 80 / 82 | 434 / 428 | 1.012 / 1.012 | 0 |
| LogNormal | data | 83 / 79 | 400 / 461 | 1.006 / 1.008 | 0 |
| Wald | constant | 107 / 108 | 520 / 510 | 1.010 / 1.008 | 0 |
| Wald | data | 105 / 115 | 388 / 563 | 1.020 / 1.009 | 0 |

Warmup moves by under 10% in either direction. Minimum ESS moves by about
as much as one seed differs from the next, in both directions. R-hat and
divergences are indistinguishable. A start several log units per
observation closer to the data is absorbed by 500 warmup iterations without
trace, even when the constant misses the threshold by a factor of four.

## The cluster ablation (Artemis, 2026-09-22)

The laptop runs above had two seeds and one data size, so the question was
re-asked at scale with `benchmarks/inits_ablation/` (a SLURM array, one cell
per (data, family, warmup, scheme, seed), CmdStan driven directly so that the
WARMUP draws are kept). The metric a better start could move is the
cold-start transient the Illusion Game project measured: the first 100-150
warmup iterations at maximum treedepth before the first metric window. So
each cell records the total leapfrog steps spent in warmup, the transient
length (first warmup iteration at which `lp__` reaches the chain's own
post-warmup 5th percentile, worst chain), the fraction of warmup iterations
at treedepth 10, ESS per second with and without warmup, R-hat and
divergences. Four chains everywhere. The constant arm is `cogmod_inits()` as
it was before the layer, rebuilt from the same internals.

Wall-clock is NOT a usable metric on a shared cluster: identical fits
differed threefold with node load (two 480-participant cells doing 635k and
630k leapfrog steps took 3.1 h and 1.9 h). Leapfrog counts are.

**Experiment A** - real lexical-decision RTs (speed_acc): 10 participants x
100 trials, the same shifted to 0.3 + 1.1 RT, and all 31,175 trials of all 17
participants; LogNormal, Wald, ex-Gaussian, Weibull; warmup 150 and 500; 8
seeds per cell; 384 cells, 8 CPUs each. Ratios data / constant, paired over
seeds (* p < .05, ** p < .01):

| data | family | warmup | warmup leapfrog | transient | min ESS / total s | median ESS / s | max R-hat |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 10 ppts | LogNormal | 150 / 500 | 1.02 / 1.00 | 1.10 / 1.12 | 0.91 / 1.10 | 0.88 / 0.97 | 1.00 / 1.00 |
| 10 ppts | Wald | 150 / 500 | 0.98 / 1.01 | 0.82* / 0.83* | 0.89 / 1.05 | 0.92 / 1.04 | 1.00 / 1.00 |
| 10 ppts | ex-Gaussian | 150 / 500 | 1.01 / 1.02 | 1.11 / 1.11 | 0.80 / 1.03 | 1.01 / 0.95 | 1.00 / 1.00 |
| shifted | LogNormal | 150 / 500 | 1.03 / 1.02 | 0.93 / 0.94 | 0.95 / 1.13* | 0.93 / 1.11 | 1.00 / 1.00 |
| shifted | Wald | 150 / 500 | 1.13** / 1.07** | 0.64** / 0.63** | 0.82* / 0.88 | 1.15 / 0.96 | 1.01 / 1.00 |
| shifted | ex-Gaussian | 150 / 500 | 1.01 / 1.00 | 1.06 / 1.05 | 0.96 / 0.99 | 1.06 / 1.04 | 1.01 / 1.00 |
| all 31k | LogNormal | 150 / 500 | 1.02 / 1.03 | 1.04 / 1.05 | 1.21 / 0.79 | 0.85 / 0.99 | 0.99 / 1.00 |
| all 31k | Wald | 150 / 500 | 1.04 / 1.02 | 0.84* / 0.84* | 0.87 / 0.76 | 1.02 / 0.77 | 1.01 / 1.00 |
| all 31k | ex-Gaussian | 150 / 500 | 0.96 / 0.99 | 1.10 / 1.11 | 1.21 / 0.56 | 1.59 / 0.62 | 0.92 / 0.94 |

Over the 18 sound-family groups the data start spent FEWER leapfrog steps in
warmup in 3, the median ratio being 1.01; ESS per second was a coin flip in
both directions; R-hat and divergences were indistinguishable. The one
consistent effect is the Wald's: its data start reaches the typical set in
16-37% fewer iterations on every data set, significantly, and then spends
the saved iterations on a step size the sampler has to walk back down from
(shifted data: 7-13% MORE warmup leapfrog, also significant). A shorter
transient is not a cheaper warmup.

The Weibull rows are excluded from the table because both arms failed on the
shifted and the full data: R-hat 1.5-2.3 and a minimum ESS of 5-7 on most
seeds, chains at treedepth 10 throughout, and six constant-arm cells still
unfinished at the 2 h wall. The shape and the non-decision time trade off on
real RTs and the chains sit in different modes; that is an identifiability
finding about the family, the same with either start, and worth its own
look.

One asymmetry did appear, and it is about robustness rather than speed. On
the full 31k-row data, of the sound families' 48 constant-start fits, 4 went
wrong - the ex-Gaussian at seeds 4 and 6, under both warmup lengths, so the
jittered start decided it: two stuck at R-hat 1.6 with 120 and 458
divergences, two past the 2 h wall of `short` - and resubmitted with 8 h, the
warmup-150 one finished as another stuck fit (R-hat 1.56, ESS 7, 321
leapfrog steps per draw) while its data-start twin had taken 1.5 h to a clean
one - against 0 of 48 data-start fits. Fisher's exact test puts that around p = 0.06; one family, one data
set. It is the kind of effect a start can plausibly have on a large data set
whose constants sit far off (median RT 0.56 s against a constant `mu` of
0.4 s and `tau` of 0.2 s), and it would want confirming before it carried
any weight.

**Experiment B** - the Illusion Game Muller-Lyer data with the production
formula pattern (a 2-D tensor smooth of difficulty x illusion strength plus
a participant intercept on every distributional parameter, `poutlier` with a
participant intercept), at 120 and 480 participants (15k and 61k rows),
LogNormal / Wald / ex-Gaussian, warmup 300 and 1000, 3 seeds, 4 chains x 4
threads. Stopped by decision after 39 of 72 cells: the 480-participant cells
took 3 h each at warmup 300 and the array was holding the account's whole CPU
allowance. At 120 participants, with the transient at 27-35% of warmup
iterations at treedepth 10 (the regime the design was aimed at), the warmup
leapfrog ratios were 0.99-1.04 for all three families at both warmup lengths;
ESS per second favoured the constant arm slightly (median ratio 0.85, 3
seeds). The three finished 480-participant cells agree: LogNormal warmup 300,
constant 635k against data 630k leapfrog steps.

## Conclusion

The kill criterion the plan set - no change in warmup or ESS on the Tier 1
ablation - is met, on the laptop and on the cluster alike. For the RT-only families the constants were never the
problem: the sampler walks in from them in the first few dozen iterations,
and what the earlier stuck-chain work fixed was a start on a *flat region*
(`ndt` above the data, a shape below 1), which `.ndt_start()` and the
constants already avoid. Starting closer to the mode is a different thing
from starting off a plateau, and only the second one matters.

What this does not rule out is the race models. Their stuck chains came from
cold *drift* starts flung along a plateau in the first transition, which is
the flat-region failure, not the distance one; a per-accumulator rule there
would be tested against that, and the RDM ablation is the test. That is
Tier 3 of the plan, and it should be run on its own terms rather than
inferred from these numbers.

## Decision (2026-09-22)

Parked. The layer was removed from the tree the same day: it added about 400
lines across `R/cogmod_inits.R`, `R/core_shifted.R` and the tests for a
change no fit could see. The implementation is kept as
`inits_ablation/data_aware_inits.patch` (apply with `git apply` from the
repository root against the 0.3.3 dev tree it was cut from), the harness that
measured it in `inits_ablation/`, and this note. The comment above
`.ndt_start()` points here.

Two ideas survive the decision and are worth keeping in mind:

- **The EZ-diffusion inversion** (Wagenmakers, van der Maas & Grasman, 2007)
  was the plan's Tier 2 for `cogmod_ddm()` and was never built. It maps the
  proportion of upper-boundary responses, the mean RT and the RT variance onto
  drift, boundary and non-decision time in closed form, and it remains
  interesting for reasons that have nothing to do with a start: as a
  data-check or a summary in `cogmod_checkdata()` / a `summary()`-style
  helper, as a sanity comparison against a fitted DDM's population-level
  posterior, or as a first pass over many participants before a hierarchical
  fit. dRiftDM's `get_ez_diffusion()` (MIT) has the inversion with its edge
  cases; bmm ships an EZDM *family* (GPL-2, read only) so a family is not the
  gap. The measured null above says only that a better *start* does not help
  the sampler, not that the inversion is not useful.
- **The robustness question.** The one asymmetry seen - 4 of 48 constant-start
  fits of the sound families failing on the 31k-row data against 0 of 48 data
  starts, all four the ex-Gaussian - is not what the plan set out to buy and
  rests on four fits of one family. If anyone returns to this, that is the question to ask directly:
many seeds, the ex-Gaussian and the Wald, 30k+ rows, counting failed fits
rather than timing successful ones. The harness in `inits_ablation/` does it
as written; the design table is the only thing to change.

Recipe for the laptop fit comparison, should it be wanted again: build the model
with `chains = 0`, then `update(m0, chains = 4, init = <scheme>, seed = s,
recompile = FALSE)` for each scheme and seed, reading warmup and sampling
time from `rstan::get_elapsed_time(fit$fit)` and ESS / R-hat from
`posterior::summarise_draws()` over the population-level parameters. The
constant scheme is the pre-change `cogmod_inits()` body: `.init_targets()`,
`.ndt_start()` on `sdata$Y`, `.init_plan()`, `.init_fun()`.
