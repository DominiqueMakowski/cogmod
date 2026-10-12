# Priors proposed by the Illusion Game prior check, refitted (2026-10-12)

A check of 18 Illusion Game fits (324,000 trials each) against the priors they
were sampled with proposed six changes to cogmod's defaults:

- **R1:** priors on `boundary` for the DDM and the Wald families.
- **R5:** priors on the DDM's `bias` slopes.
- **R6:** `poutlier` centred at -6 rather than -5.
- **R2-R4:** the race families' `driftone` priors, whose intercept and slopes
  sat in the far tails, centred on the error accumulator and widened.

A prior is a sensitive change (AGENTS.md), so each one was refitted on the
vignettes' data before shipping, current against proposed.

**Outcome.** R1, R5 and R6 shipped in 0.3.5: no posterior moved beyond noise,
and sampling was as cheap or cheaper. R2-R4 did not ship.

- They predicted the same error rates and RTs.
- They cost 9-26% more leapfrog steps per draw on the RDM.
- On the LBA2 they raised the divergent transitions from 1 in 12 chains to 5
  in 12 and 10 in 11.

The tail position they were meant to fix is structural. Below a drift of a few
hundredths the likelihood is flat, so the posterior of a fenced plateau sits
in the lower tail of whatever prior fences it, and where in the plateau it sits
changes no prediction. The documentation now says so instead.

## Design

| arm | `driftone` intercept / slopes | everything else |
| --- | --- | --- |
| A | RDM `normal(3, 2)` / `normal(0, 2)`, LBA2 `normal(1, 2)` / `normal(0, 1.5)` | 0.3.4 |
| B | `normal(-2, 3)` / `normal(0, 4)` (the proposal) | R1, R5, R6 |
| C | `normal(-2, 2)` / `normal(0, 3)` (a tighter fallback) | as B |

The models are the vignettes' own formulas: 500 iterations, 250 of them warmup,
4 chains, and `dense_e` for the LBA2.

- **`rdm`, `lba2`:** the decision-making article's data (two participants of
  rtdists' `speed_acc`, 4,620 trials, errors kept).
- **`rdm_sub`:** the stress case, a 70% subset with 71 Accuracy errors. That
  is the ~75 on which a flat prior once ran the RDM's drift to the floor with
  half the transitions divergent.
- **`ddm`:** the 4-parameter DDM.
- **`ddm_bias`, `ddm_pout`:** the same DDM with `bias ~ Condition` or with
  `poutlier ~ 1` added.
- **`wald4`, `exwald`, `bisa`:** the RT article's fits (three participants,
  correct responses).

Seeds: 4 per arm for the RDMs, 3 for the LBA2 and 2 for the rest. A seed fixes
both the sampler and `cogmod_inits()`'s jitter, so every arm starts from the
same values. The cells ran on Artemis, one fit on 4 CPUs each, from a frozen
0.3.4 tree (5382c79). Arms B and C patch the registries in-session, so the
Stan programs differ only in their priors.

## Results

Efficiency is bulk ESS per 1,000 leapfrog steps, as a ratio to arm A, over
every parameter (`summarise.R` prints them one by one).

| model | arm | divergent (chains) | max R-hat | leapfrog / draw | ESS per step vs A |
| --- | --- | --- | --- | --- | --- |
| rdm | A | 0 (16) | 1.012 | 39.4 | - |
| rdm | B | 0 (16) | 1.018 | 49.5 | 0.68-0.85 |
| rdm | C | 0 (16) | 1.018 | 43.4 | 0.80-0.92 |
| rdm_sub | A | 0 (16) | 1.010 | 38.1 | - |
| rdm_sub | B | 0 (16) | 1.011 | 46.5 | 0.50-0.83 |
| rdm_sub | C | 0 (16) | 1.012 | 41.6 | 0.49-0.93 |
| lba2 | A | 1 (12) | 1.009 | 15.0 | - |
| lba2 | B | 10 (11) | 1.007 | 16.2 | 0.48-0.72 |
| lba2 | C | 5 (12) | 1.016 | 15.5 | 0.63-0.92 |
| ddm | B | 0 (8), A 0 | 1.006 | 12.6 (A 13.7) | 0.99-1.15 |
| ddm_bias | B | 0 (8), A 0 | 1.019 | 13.7 (A 13.4) | 0.89-1.14 |
| ddm_pout | B | 0 (8), A 0 | 1.005 | 13.0 (A 13.1) | 0.83-1.05 |
| wald4 | B | 0 (8), A 0 | 1.010 | 118.8 (A 149.5) | 1.38-1.81 |
| exwald | B | 0 (8), A 0 | 1.210 | 44.3 (A 48.7) | 0.91-1.25 |
| bisa | B | 0 (8), A 0 | 1.014 | 25.1 (A 24.6) | 0.88-1.28 |

**What B and C did to the race families.** Arms B and C moved `driftone` down
the plateau and the condition slope up to pay for it:

- **RDM, Accuracy intercept:** from -3.97 to -5.79 (B) and -5.05 (C); slope
  from 4.58 to 6.40 and 5.68.
- **LBA2:** the intercept went from -4.71 to about -6.1 and -6.6, and
  `sigmaone` climbed the ray from 2.40 to 2.62-2.70 on its link.
- **The priors' tail flags went away:** the centred intercept's prior quantile
  rose from 0.01-0.03 to 0.18-0.47.
- **What the decision process predicts did not change:** error rates per
  condition to within 0.003 and median correct and error RTs to within 5 ms,
  in every arm, against observed 4.6% and 9.8% errors.

**What the other changes did.**

- **DDM:** no decision parameter's posterior median moved by more than 0.006
  on its link.
- **`poutlier` (R6):** on the link it went from -5.52 to -5.58, a fifth of a
  posterior SD. As a natural-scale proportion it moved by 0.0001-0.0002.
- **The Wald with drift variability (`wald4`):** the `boundary` prior trimmed
  the upper end of the (`mu`, `boundary`, `sigmadrift`) ray, the boundary's
  95th percentile going from 2.7 to 2.1. The fit got cheaper.
- **Gradient cost per evaluation** (`gradient_cost.R`, `--reps 21`, against
  0.3.4): bisa 1.00, exwald 1.00, invgaussian 0.99, lba2 0.96, rdm
  1.03.
  - The DDM was not timed. The harness frees all seven of its parameters,
    whose numerical integration makes a gradient tens of seconds at 5,000
    trials, against which three scalar `normal_lpdf()` terms cannot show.
- **Gradient check:** passes on all six touched families, apart from the DDM's
  known `wiener_lpdf()` waiver.

## Two findings for later

- **`cogmod_lba2_decision_lpdf()` can return NaN.** Under arm B's prior the
  LBA2 sampler reached proposals the density evaluated to NaN
  (`log_mix: lambda2 is -nan`). They were rejected, not fatal, and the reruns
  of seeds 1 and 2 were deterministic across two node types (EPYC 7513 and
  9355): identical draws, 7 and 2 divergent transitions.
  - On the first run those two cells, and one chain of seed 3, failed with "No
    chains finished successfully". The messages were suppressed and the CmdStan
    output had been written to scratch and deleted, so the cause is not
    recorded. The same computation succeeded on the same node an hour later.
  - The NaN is a latent bug in the density. Arms A and C were not run with
    messages, so whether they reach it too is unknown. The gradient check's
    grid does not get there.
  - Fixed the same day (NEWS.md, 0.3.5, "Bug fixes"): `phi(z1) - phi(z2)`
    overflowed to 0 * inf once an accumulator's drift should have carried it
    past the threshold well before the response - here an error accumulator
    with a positive drift and a small SD. Stan counts a NaN proposal as a
    divergent transition, so some of arm B's LBA2 divergences may be the bug.
- **`cogmod_exwald()` can split modes on the RT article's data.** Seed 1 left
  one chain about 11 log-density units below the other three in *both* arms
  (R-hat 1.21-1.23, minimum ESS 3-4), so it is not the priors. The cached
  vignette fit happens to agree across chains.

## Files

| file | what |
| --- | --- |
| `prep.R` | writes `push/`: the vignette data as plain data frames, the 57-cell table, the package tree, the two scripts |
| `refit_cell.R` | one cell: patches the registries for the arm, emits the program, samples, writes the chain diagnostics, parameter summaries, posterior predictions and priors |
| `task.slurm` | one array task: modules, libraries, a scratch directory per cell |
| `run.sh` | `push`, `smoke` (one cell with srun), `submit` (the array, 32 cells at a time), `pull` |
| `collect.R` | stacks the pulled per-cell CSVs into `results/` |
| `summarise.R` | the tables above, from `results/` |
| `results/chains.csv` | per chain: leapfrog steps in warmup and sampling, divergences, treedepth hits, step size, mean lp, E-BFMI, time, node |
| `results/pars.csv` | per parameter and cell: median, 5% and 95%, R-hat, bulk and tail ESS, chain spread |
| `results/pp.csv` | per condition and cell: predicted error rate and median correct and error RT, from 100 draws x 2,000 trials |
| `results/priors.csv` | the priors each arm gave each model |

To rerun from the package root:

```
Rscript benchmarks/driftone_priors/prep.R
benchmarks/driftone_priors/run.sh push
benchmarks/driftone_priors/run.sh smoke
benchmarks/driftone_priors/run.sh submit
benchmarks/driftone_priors/run.sh pull
Rscript benchmarks/driftone_priors/collect.R
Rscript benchmarks/driftone_priors/summarise.R
```

Arm A is whatever tree `prep.R` copies. The comparison above used 0.3.4's;
`git archive 5382c79` gives it back. Most cells took 2-15 minutes; the RDM's
took about 40.
