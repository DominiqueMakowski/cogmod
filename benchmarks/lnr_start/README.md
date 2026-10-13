# Where should the LNR's error accumulator start? (2026-10-13)

`cogmod_inits()` started `cogmod_lnr()` with both accumulators at the same
rate, `nuone = mu = 0.7`: a 50/50 race, 50% errors under the `dec()` coding
every example uses (1 = error, `nuone` the error accumulator). The other two
race families had already moved their error accumulator below the correct
one - `cogmod_rdm()` in 0.3.2, `cogmod_lba2()` in 0.3.4 - and a check of
0.3.5's starts against the Illusion Game posteriors (which sit at `nuone`
-0.8 to -1.2 against `mu` 0.7-1.15) proposed `nuone = 0` for the LNR: with
both sigma at 0.5 that implies 16% errors, close to the LBA2's 11%, and sits
0.47 SD below the centre of `nuone`'s `normal(0.7, 1.5)` prior. The prior does
not move. The check traced no failure to the start, so the question is only
whether the cheaper start is worth having, measured as the LBA2's was.

**Outcome.** `nuone = 0` shipped in 0.3.5, for consistency with the RDM and
the LBA2 rather than for speed. Where errors are rare (3-12%), the closer
start cuts the shortfall of the cold start by 750-1,600 log-likelihood units
per 2,000 trials (step 1), but that head start does not reach the sampler.
Refitted cold, it
changed the median warmup cost by 8% or less, cheaper on one model and dearer
on another, with no divergences, stuck chains or moved posteriors under
either start (step 2).

## Step 1: the shortfall of the cold start (`shortfall.R`)

Log-likelihood at the generating values minus log-likelihood at the whole
cold start - every parameter at its 0.3.5 start, `ndt` at half the data's
first percentile, `poutlier` 0.02 - with `nuone` moved through 0.7 (the
current start), 0.35, 0, -0.35 and -0.7. 2,000 simulated trials, mean of 5
seeds; `mu` generated at 1.0 with two sigma pairs, equal (0.5, 0.5, the
start's own) and wide (0.6, 1.2, the ratio the vignette and the Illusion Game
fits show); `sigmabias` fixed at 0, and free (started at 0.5, generated at
0.2, where the Illusion Game posteriors sit).

Shortfall in log-likelihood units per 2,000 trials, `sigmabias` fixed:

| sigmas | errors | true `nuone` | start -0.7 | -0.35 | 0 | 0.35 | 0.7 (current) | gain of 0 over 0.7 |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| equal | 3% | -0.33 | 307 | 334 | 540 | 1,070 | 2,087 | +1,547 |
| equal | 12% | 0.17 | 608 | 364 | 353 | 706 | 1,575 | +1,222 |
| equal | 26% | 0.55 | 1,482 | 798 | 413 | 462 | 1,079 | +666 |
| equal | 43% | 0.88 | 3,286 | 1,965 | 932 | 400 | 496 | -436 |
| wide | 3% | -1.52 | 203 | 246 | 476 | 1,036 | 2,084 | +1,608 |
| wide | 12% | -0.58 | 535 | 309 | 272 | 594 | 1,435 | +1,163 |
| wide | 26% | 0.14 | 1,391 | 818 | 380 | 324 | 833 | +453 |
| wide | 43% | 0.76 | 3,123 | 2,185 | 1,188 | 483 | 355 | -833 |

With `sigmabias` free the picture is the same (`results/shortfall_sim_summary.csv`):
gains of 750-1,200 at 3-12% errors, 20-220 at 26%, and losses of 880-1,200 at
43%.

So the 50/50 start is the costliest of the five at every error rate up to
26%, and `nuone = 0` takes 1,150-1,600 units per 2,000 trials off the
shortfall at 3% errors, 750-1,220 at 12% and 20-670 at 26%. It is not the
LBA2's picture: there the equal start was the costliest *every* time, because
a truncated-Normal race is thin on the fast side. The LNR's cost is roughly
symmetric about the truth, so at 43% errors - a rate almost no two-choice
data set has - the proposed start is 440-1,200 units worse than the current
one. 0.35 would hedge that (never more than 370 worse, 500-1,300 better
elsewhere) at the cost of half the gain where errors are rare; 0 is kept
because it is the counterpart of the other races' starts and the centre of
where real posteriors sit.

On the decision_making article's data (4,620 trials, 7.3% errors), with the
fitted model's posterior medians as the reference:

| start `nuone` | implied errors | shortfall | gain over 0.7 |
| --- | --- | --- | --- |
| 0.7 (current) | 50% | 2,767 | 0 |
| 0.35 | 31% | 1,610 | +1,157 |
| **0** | **16%** | **1,302** | **+1,465** |
| -0.35 | 7% | 1,482 | +1,285 |
| -0.7 | 2% | 1,836 | +931 |
| -1.0 | 1% | 2,121 | +646 |

## Step 2: does it change a fit? (`prep.R`, `fit_arm.R`, `summarise.R`)

Cold starts, arm A (`nuone` at 0.7) against arm B (`nuone` at 0, patched into
the registry in-session), 4 seeds x 4 chains, 250 warmup + 250 sampling
iterations, `diag_e`, CmdStan driven directly with the warmup draws kept. A
seed fixes the jitter as well as the sampler, and the jitter is additive on
the link scale, so the two arms' chains differ only in the centre of
`Intercept_nuone` (`summarise.R` checks this from the recorded starts).

| model | data | formula |
| --- | --- | --- |
| `vignette` | the decision_making article's (4,620 trials, 7.3% errors) | `RT \| dec(Error) ~ Condition, nuone ~ Condition, ndt ~ Condition, sigmabias = 0` |
| `vignette_bias` | the same | the same with `sigmabias ~ 1` |
| `sim_hier` | simulated: 20 participants x 120 trials, errors 11.5% in Accuracy and 25.0% in Speed | `~ Condition + (1 \| Participant)` on both rates, `ndt ~ 1 + (1 \| Participant)`, `sigmabias = 0` |

Per chain: leapfrog steps in warmup and in sampling, divergences, treedepth
hits, the adapted step size, `lp__` at the start and the cold-start
transient (first warmup iteration at which `lp__` reaches the chain's own
post-warmup 5th percentile), starts rejected by CmdStan, time; per fit:
R-hat and bulk ESS over the population-level parameters, and bulk ESS per
second of wall time.

The 24 cells of record ran on Artemis as SLURM array 11596604 (2026-10-13),
13 on artemis-rtx-00 and 11 on artemis-rtx-02. Each was one fit of 4 chains
on 4 whole cores (`--hint=nomultithread`), with CmdStan 2.39.0 and brms
2.23.1, from a frozen copy of this tree. All 96 chains finished.

A first pass on a 12-core laptop (CmdStan 2.38.0, the two arms side by side)
got through the `vignette` model and one seed of each of the others before
the laptop had to be shut. Those 12 cells are kept under `results/local/`,
were not mixed with the cluster's, and agree with them.

**Two earlier arrays lost chains, and why.** Arrays 11596401 and 11596580
lost 16 and 20 of their 96 chains, every one SIGKILLed (exit status -9) with
nothing on stderr.
- **The tag.** processx 3.8.2, the copy in Artemis's R module, tags each
  chain process with an id made of `sample()` from R's RNG plus the current
  second.
- **The kill.** When `$sample()` returns, cmdstanr cleans up with processx's
  `kill_tree()`, which SIGKILLs every process on the machine that carries
  the tag.
- **The collision.** `fit_cell.R` called `set.seed(seed)` to draw the starts,
  and the two arms of a seed consume the RNG identically. So whenever both
  arms started sampling on the same node in the same second, their chains
  shared tags, and the first arm to finish killed the other's running
  chains. In 11596401 (all 24 cells on one node), 8 of the 9 pairs that
  logged the same start second lost chains, and none of the 3 that did not.
- **The fix.** `fit_cell.R` now calls `set.seed(NULL)` once the starts are
  drawn. CmdStan's seed is passed explicitly, so the draws do not change.
- **The bookkeeping.** cmdstanr drops a chain it did not see finish and
  renumbers the rest, which in 11596401 misaligned the per-chain times and
  hid the loss. `fit_cell.R` now rebuilds each fit from CmdStan's CSVs and
  records every chain's exit status and draw count.

Neither earlier array is used. 11596401's cells are under
`results/cluster_11596401/`; 11596580's, with the killed chains' files, are
on the cluster in `results_11596580/`. Any harness that calls `set.seed()`
before `$sample()` in several tasks on one node, with processx below 3.9, is
exposed in the same way.

### Results

`Rscript benchmarks/lnr_start/summarise.R`, array 11596604. Per arm, over the
4 seeds: divergences after and during warmup (summed), the largest R-hat,
then medians. Leapfrog steps are per draw per chain; the transient is the
worst chain's.

| model | arm | div | div in warmup | max R-hat | warmup leapfrog | sampling leapfrog | step size | transient | min ESS | wall (s) |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| `vignette` | A (0.7) | 0 | 312 | 1.017 | 47.3 | 18.9 | 0.235 | 27.5 | 416 | 27 |
| `vignette` | B (0) | 0 | 291 | 1.020 | 44.4 | 20.6 | 0.223 | 36.0 | 428 | 24 |
| `vignette_bias` | A (0.7) | 0 | 343 | 1.008 | 96.6 | 62.7 | 0.068 | 29.5 | 396 | 282 |
| `vignette_bias` | B (0) | 0 | 336 | 1.008 | 102.2 | 60.4 | 0.067 | 27.5 | 376 | 298 |
| `sim_hier` | A (0.7) | 0 | 322 | 1.011 | 171.0 | 32.4 | 0.155 | 11.5 | 262 | 44 |
| `sim_hier` | B (0) | 0 | 315 | 1.012 | 177.6 | 32.2 | 0.156 | 10.0 | 252 | 41 |

No chain froze and no start was rejected, in either arm.

B / A paired by seed, median (range over the 4 seeds). Below 1 is cheaper or
shorter, except for ESS per second, where above 1 is better:

| model | warmup leapfrog | sampling leapfrog | transient | min ESS per second | wall |
| --- | --- | --- | --- | --- | --- |
| `vignette` | 0.92 (0.89-0.98) | 1.06 (0.95-1.14) | 1.21 (1.03-1.42) | 1.30 (0.85-1.88) | 0.83 (0.75-1.07) |
| `vignette_bias` | 1.08 (0.98-1.13) | 0.97 (0.82-1.08) | 0.93 (0.65-1.07) | 0.88 (0.71-1.00) | 1.06 (0.91-1.16) |
| `sim_hier` | 1.03 (0.99-1.05) | 1.00 (0.97-1.03) | 1.01 (0.75-1.83) | 1.02 (0.56-1.72) | 0.96 (0.91-1.06) |

The posterior medians of the two arms agree within 0.06 posterior SDs on
every parameter of every model (`results/agreement.csv`).

What this says:
- **Neither start is less reliable.** The arms are identical on every
  failure count.
- **The closer start barely changes cost.** With `sigmabias` fixed, its
  warmup is cheaper in all 4 seeds (by 8% at the median), and the sampling
  phase gives most of that back. With `sigmabias` free, its warmup is dearer
  in 3 seeds of 4 (by 8% at the median). On the hierarchical model there is
  no difference.
- **The head start is spent at once.** Every chain of either arm reaches the
  typical set within 37 iterations. Step 1's 1,465-unit head start on this
  data is gone in the first few dozen iterations, and warmup spends its time
  adapting, not travelling.

So there is no performance case either way, and 0 was chosen for
consistency with the other two races (the registry comment in
`R/core_choice.R` has the reasoning).

## Files

| file | what |
| --- | --- |
| `shortfall.R` | step 1; writes `results/shortfall_sim.csv`, `results/shortfall_sim_summary.csv`, `results/shortfall_vignette.csv` |
| `prep.R` | the data, the formulas and the three programs, compiled once into `work/` for a local run, plus `push/` for the cluster (both gitignored) |
| `fit_cell.R` | one cell: compiles, sets the arm's `nuone` start, samples, writes `results/<model>__<arm>__s<seed>__{chains,pars,inits}.csv`; a row per chain whether it finished or not, and `results/failed/` for any chain that did not |
| `fit_arm.R` | one arm's cells on the laptop, in one process, through `fit_cell()` |
| `task.slurm`, `run.sh` | one array task, and `push` / `smoke` / `submit` / `queue` / `log` / `pull` over SSH |
| `summarise.R` | the tables above, after a count of the chains that finished; paired ratios use only seeds where both arms kept all 4; writes `results/fits.csv` and `results/agreement.csv` |

```
Rscript benchmarks/lnr_start/shortfall.R
Rscript benchmarks/lnr_start/prep.R
benchmarks/lnr_start/run.sh push
benchmarks/lnr_start/run.sh submit
benchmarks/lnr_start/run.sh pull
Rscript benchmarks/lnr_start/summarise.R
```

Or on the laptop, `Rscript benchmarks/lnr_start/fit_arm.R A` and `... B` in
two shells instead of the three `run.sh` steps.
