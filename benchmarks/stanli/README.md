# stanli as a backend for cogmod: exact, and faster per gradient and end to end

2026-10-02 to 04 and 2026-10-12, on two platforms:

- **Windows**: a laptop (i7-1265U: 2 performance + 8 efficiency cores),
  CmdStan 2.38, cmdstanr 0.9.0.9000, R 4.5.3.
- **Linux**: Artemis, the Sussex HPC (AMD EPYC 9334, 9355 and 7513 nodes,
  Rocky 9.5), CmdStan 2.39, cmdstanr 0.9.0, R 4.3.2.

Both run brms 2.23.1 from GitHub (`e71e9d7`; CRAN has 2.23.0), so both sample
the same Stan program. stanli 0.18.1 first, then 0.19.0 (released 2026-10-03,
in response to [seantalts/stanli#422](https://github.com/seantalts/stanli/issues/422)),
then 0.19.1 (2026-10-04) together with a hand-written LNR program the
maintainer posted as [cogmod#5](https://github.com/DominiqueMakowski/cogmod/issues/5)
(section "0.19.1, and the rewrite from cogmod#5"; on Linux these ran on EPYC
9355 and 7513 nodes), then 0.20.0 (2026-10-11), which fixed what our 0.19.1
report raised (section "0.20.0").

| file | what |
| --- | --- |
| `bench.R` | the harness: `emit`, `grad` (exactness and cost per gradient), `fit` (sampling, 2 x 2) |
| `bisect.R` | cost of each way of writing the likelihood, stanli only |
| `repro_issue.R` | stanli-only reproducer of #422 |
| `repro_parallel_init.R` | stanli-only reproducer: chains ran serially when given inits (0.19.x; fixed in 0.20.0) |
| `multilevel.R` | seven multilevel decision models (LNR, LNR-6, DDM-4/5/7, LBA, RDM) on simulated data: `emit`, `fit` (cmdstanr against stanli), `grad` (exactness and cost at posterior draws); section "Multilevel decision models" |
| `chains_probe.R` | one chain against four, as threads and as processes, in both engines |
| `probe_bottlenecks.R` | per-gradient cost of DDM-4 and LNR variants that each remove one suspect (Windows) |
| `postprocess.R` | which post-processing tools work on a fit through the `"mock"` route, against a cmdstanr fit (`results/postprocess_0.20.0.csv`) |
| `lnr_select_functions.stan`, `lnr_branchfree_functions.stan` | the rewritten likelihoods |
| `lnr_rewrite.stan` | cogmod#5's whole program, verbatim (variant `rw`) |
| `scalar_dpar.R` | brms's `sigma ~ 1` as a vector against a real, on a built-in family (gaussian); the script behind brms#1945 |
| `brms_scalar_dpar.patch` | the change as submitted in brms#1946, against brms f131ef1 |
| `brms_scalar_dpar/` | its harness: `smoke.R` (generated code per formula), `emit.R` and `compile_check.R` (eleven models, installed brms against the patched checkout) |
| `hpc/` | the same runs on Artemis: `run.sh`, `task.slurm`, `install.R`, `summarise.R` (`summarise_ml.R` and `summarise_chains.R` for the multilevel runs); `sigill_probe.R` is the per-node probe behind the SIGILL diagnosis, `affinity_probe.R` the one behind the slow-CmdStan node of the 0.20.0 fit |
| `results/` | Windows numbers: `check`, `tail`, `time` (and their `_0.19.1` and `_0.20.0` reruns; `_0.20.0_lse` is `orig` against `lse`), `bisect.csv` (0.19.0) and `bisect_0.20.0.csv` from `bisect.R`, `fit_*.csv` from `fit` (0.19.0; `fit_*_0.20.0.csv` the rerun), `repro_issue.txt`; the three `fit_*_smoke.csv` are the seed-11 smoke run behind the ridge paragraph under "Parallel chains and inits" |
| `results/hpc/` | the Linux ones: `ml_<tag>_<model>_<seed>/` and `mlchains_<tag>_<model>_<seed>/` (tags `0.20.0` on two cores per task, `nosmt` on four; `zen3`, `idle`: the probe on 7513s and on an idle 9334); before them `fit_<seed>/`, `grad_<task>/` (0.19.0), `grad_0.19.1_<task>/`, `grad_0.20.0_<task>/`, `fit_0.20.0_<seed>/`, each with a `node.txt`; `fit_0.20.0rtx01_14/` is the CmdStan-only control on the slow node; `fit_0.3.5smt_<seed>/` and `fit_0.3.5nosmt_<seed>/` the same fit on the 0.3.5 program, two cores and four |

## The question

[stanli](https://github.com/seantalts/stanli) interprets Stan programs over
precompiled Stan Math kernels: no C++ toolchain, a model builds in about 2 s.
[paul-buerkner/brms#1911](https://github.com/paul-buerkner/brms/pull/1911)
proposes it as a brms backend, and a comment there gives a "light" route that
works with brms as released:

```r
dummy <- brm(formula, data, prior = prior, stanvars = stanvars, empty = TRUE)
m <- stanli::stanli_model(code = stancode(dummy), data = standata(dummy))
fit <- stanli::sample_model(m, chains = 4, warmup = 500, samples = 500)
dummy$fit <- stanli::as_stanfit(fit)
fit <- brms:::rename_pars(dummy)   # then summary(), loo(), ... as usual
```

Does it work for a cogmod family, and how does it compare with brms's
cmdstanr backend? Measured on the LNR: the decision_making vignette's model
(`nuone ~ Condition`, `sigmazero ~ 1`, `sigmaone ~ 1`, `sigmabias = 0`,
`ndt ~ Condition`) on its data (speed_acc, participants 1-3, RT <= 2 s;
4620 trials), the same set-up as `benchmarks/lnr_vectorize`.

## Short answer (stanli 0.20.0)

- **It runs the program brms writes, unchanged, and exactly.** The log density
  is CmdStan's to within 4 ULP (identical at most points) and the gradient
  to 2e-14, including in the far tail (x = -65) where a naive branch-free
  copy's gradient is not finite. 0.20.0 gives 0.19.1's log density bit for
  bit at every point measured.
- **It is faster per gradient: 0.42x CmdStan on Windows, 0.75-0.77x on
  Linux**, as on 0.19.1 (0.41x; 0.77x on the same node type). CmdStan gains
  far more than stanli from Linux's faster maths library, so the gap
  narrows there. 0.20.0's `fast_math` adds nothing on this model.
- **It samples as efficiently per leapfrog step**: same posterior, same ESS
  per 1000 leapfrog steps (about 22.5 in every arm, 16 seeds). 0.20.0
  reproduces 0.19.0's draws exactly, seed for seed.
- **End to end it is now faster too.** 0.19.x ran the chains one at a time
  whenever `sample_model()` got inits, so a 4-chain fit from
  `cogmod_inits()` took 2.3x CmdStan's wall time on Linux. We reported it;
  0.20.0 runs them in parallel. The same fit now takes **0.74x CmdStan's
  wall time on Linux** (1.20x the ESS per wall second, 16 seeds) **and
  0.43-0.46x on Windows** (2.7-2.8x).
- **The four full-program cliffs we reported are gone** (a data-only
  return, `log_sum_exp()`, a blend, select plus blend: 4-11x on 0.19.1,
  0.8-1.1x now). The package's own code never hit them.
- **The maintainer's hand-written LNR (cogmod#5) is exact and faster, but
  on CmdStan the speed is not in the likelihood.** Its log density is
  brms's bit for bit. All of its CmdStan gain, and half of its stanli gain,
  comes from passing the two intercept-only sigmas as reals - one softplus
  each instead of one per trial - which brms could do for any program. The
  rewritten likelihood itself costs CmdStan 0.92-1.11x and stanli 0.83-0.97x.
- **Hand-written selects (`sel`) are 1.3x faster than the package's code in
  stanli**, but they cost CmdStan 1.36x (Windows) to 1.9x (Linux) and are
  not quite the package's numerics. stanli will not produce them
  automatically (it is held to CmdStan's results). Not worth a second Stan
  source in every family.
- **On multilevel decision models the lead shrinks to about even**
  (section "Multilevel decision models", seven models, 10 participants x
  100 trials, Linux): per gradient 0.65-0.86x CmdStan on five of six, 1.79x
  on the DDM-4 (Stan Math's `wiener_lpdf()` runs slower in stanli), but per
  fit 0.79-1.32x cmdstanr's wall time (DDM-4 1.57x), ESS per wall second
  0.65-1.17x: inside its sampler stanli loses 10-50% of its gradient speed,
  CmdStan 0-10%. What stanli wins every time is the start, 1 s against
  27-43 s of compiling. On Artemis, `--cpus-per-task=4` is two cores with
  two threads each, which hurts stanli more; the 0.20.0 Linux fit below
  ran that way too, and on four whole cores the LNR there takes 0.69x
  cmdstanr's wall time (0.76x on two).
- **brms can take the result through public API** since 0.20.0:
  `brm(..., backend = "mock", mock_fit = as_stanfit(fit))` gives the same
  `brmsfit` as the internal `rename_pars()` route.

## Multilevel decision models: faster per gradient, no faster per fit

Everything below this section is one model: the vignette's LNR, with
population-level effects only, on 4620 trials. Users fit multilevel models
on fewer trials per participant, so on 2026-10-12 `multilevel.R` put seven
of them through both backends, each fitted to data simulated from its own
RNG (outlier-free; 10 participants x 2 conditions x 50 trials, N = 1000;
error rates 7-16% in Accuracy, 11-24% in Speed). The parameter that carries
the choice gets `~ Condition + (Condition | Participant)`, the threshold
(or the second accumulator) and `ndt` a random intercept, the rest are
intercepts:

| model | family | estimated |
| --- | --- | --- |
| `lnr` | `cogmod_lnr()` | `sigmabias = 0` (the vignette's LNR) |
| `lnr6` | `cogmod_lnr()` | all six parameters |
| `ddm4` | `cogmod_ddm()` | `sigmadrift = sigmabias = sigmandt = 0` |
| `ddm5` | `cogmod_ddm()` | `sigmadrift` |
| `ddm7` | `cogmod_ddm()` | `sigmadrift`, `sigmabias`, `sigmandt` |
| `lba` | `cogmod_lba2()` | `sigmazero = 1` |
| `rdm` | `cogmod_rdm()` | all |

Per (model, seed) task: a fit through each backend, 4 chains x (500 + 500)
from the same `cogmod_inits()` and seed, cmdstanr sampling directly (as
`brm()` does) and stanli through the light route, both ending in
`rename_pars()`; then, in a separate process, CmdStan's and stanli's log
density and gradient at 40 posterior draws of the cmdstanr fit, and the
cost of one gradient there, single-threaded, in 15 alternating blocks.
Artemis, CmdStan 2.39, stanli 0.20.0, 8 seeds per model on EPYC 7513 and
9355 nodes, 4 whole cores per task (`--hint=nomultithread`; see "Two
chains per core" below). `hpc/summarise_ml.R --tag nosmt`,
`results/ml_summary_nosmt.csv`; ratios are stanli over cmdstanr within a
seed, median [min, max] over seeds:

| model | per gradient | per leapfrog, in the fit | fit wall time | ESS per wall second | build: CmdStan / stanli |
| --- | ---: | ---: | ---: | ---: | ---: |
| `lnr` | 0.86 [0.82, 0.91] | 1.04 [0.77, 2.08] | 1.06 [0.85, 2.00] | 1.04 [0.61, 1.32] | 31 s / 1.1 s |
| `lnr6` | 0.81 [0.80, 1.14] | 0.97 [0.73, 1.69] | 1.01 [0.83, 2.17] | 1.17 [0.48, 1.92] | 27 s / 0.9 s |
| `ddm4` | 1.79 [1.47, 2.32] | 1.56 [1.24, 2.46] | 1.57 [1.43, 2.25] | 0.65 [0.41, 0.85] | 39 s / 0.7 s |
| `ddm5` | 1.04 [1.02, 1.11] | 1.10 [0.96, 1.58] | 1.17 [1.14, 1.51] | 0.73 [0.56, 1.02] | 43 s / 0.9 s |
| `lba` | 0.86 [0.84, 0.87] | 1.16 [0.86, 1.57] | 1.32 [0.64, 2.09] | 0.97 [0.57, 1.57] | 37 s / 0.9 s |
| `rdm` | 0.65 [0.64, 0.88] | 0.85 [0.54, 1.93] | 0.79 [0.65, 1.40] | 1.06 [0.82, 1.61] | 27 s / 0.9 s |

The fits took 38-210 s (median, cmdstanr). `ddm7` costs about 0.1 s per
gradient at N = 1000 and is still running (4 seeds, `long` partition).

- **Exact, as on the LNR.** At all 1920 posterior points stanli's log
  density is CmdStan's to within 27 ULP (2e-12 absolute; the RDM's 53000
  ULP are at a log density of 0.02, 1.8e-13 absolute) and its gradient to
  5e-12; the posteriors agree (largest mean difference 3.7 joint MCSEs
  over 72-85 quantities x 8 seeds), and the sampler does the same work: ESS
  per 1000 sampling leapfrogs, total leapfrogs and divergences are the
  same in both.
- **Per gradient, stanli is 0.65-0.86x CmdStan on five of six models** on
  Linux - against 0.75-0.77x for the population-level LNR - and **1.79x on
  the DDM-4**. The DDM's cost is Stan Math's `wiener_lpdf()`, which both
  call, so the gap is in how stanli's kernel runs, not in cogmod's code
  (below).
- **End to end it is about even.** A multilevel fit took 0.79x (RDM) to
  1.32x (LBA) of cmdstanr's wall time, 1.57x for the DDM-4, and ESS per
  wall second ranges from 0.65x to 1.17x. Per-seed ratios spread widely
  (0.6-2.2) on shared nodes. What stanli wins clearly is the start: no
  C++ compile, 1 s against 27-43 s, which is 13-81% of these fits'
  sampling time.
- **The lead per gradient does not survive the sampler.** Inside the fit,
  a cmdstanr leapfrog costs 0.96-1.13x its single-threaded gradient; a
  stanli one 1.04-1.51x (RDM 1.51, LBA 1.34, LNR-6 1.21, LNR 1.16). The
  population-level LNR shows no such gap (0.71x CmdStan per leapfrog in
  its fit, 0.75-0.77x per gradient). Across the multilevel fits the gap is
  weakly related to leapfrogs per iteration (Spearman 0.31; LBA 126, RDM
  64, DDM-4 32) and not to the cost of a gradient (-0.04). stanli
  also spends a little more of its time in warmup on five models of six
  (warmup over sampling time 1.86-4.03 against 1.61-3.84) with the same
  sampling-phase leapfrog counts.

**Two chains per core.** The first run (tag `0.20.0`, results kept)
asked Slurm for `--cpus-per-task=4`. Artemis allocates cores
(`CR_CORE_MEMORY`) with two hardware threads each, so that is two physical
cores: a probe job got CPUs 0, 1, 64 and 65, cores 0 and 1 with their
siblings. Four chains then share two cores, and stanli suffers from it
more than CmdStan: per leapfrog in the fit 1.55x for the LNR (now 1.04x),
wall 1.60x (now 1.06x). `--hint=nomultithread` gives four cores, one
thread each, and costs eight CPUs of the account's 550. The same draws
came out in both runs, seed for seed; only the clock moved. The 0.20.0
LNR fit (section "0.20.0") ran the same way. Rerun both ways on the
current code, 16 seeds each (`fit_0.3.5smt_<seed>/`,
`fit_0.3.5nosmt_<seed>/`; 0.3.5's priors changed the program, so `bench.R
emit` now skips `rw` and `rw_v`, whose blocks date from 0.3.4): stanli
`orig` took 0.76 [0.67, 0.88] of cmdstanr's wall time on two cores and
0.69 [0.47, 0.86] on four, with 1.38x and 1.61x its ESS per wall second.
The population-level LNR keeps its lead, a little larger on whole cores;
it is the multilevel models that lose it.

**One chain, four chains, four processes** (`chains_probe.R`, `run.sh
submit mlchains`): the sampling-phase cost of a leapfrog with one chain
alone, four chains as stanli's threads, and four one-chain R processes, in
both engines; 4 models x 4 seeds, 300 + 300 iterations,
`hpc/summarise_chains.R --tag <tag>`. On two cores per task (tags
`0.20.0` on 9355s, `zen3` on busy 7513s, `idle` one task at a time on an
otherwise idle EPYC 9334) four chains cost stanli 1.5-2.7x one chain on
the three race models, against CmdStan's 1.4-1.6x (2.1x for the RDM on the
busy 7513s), and the DDM-4 alike in both; as much with processes as with
threads, so not contention inside stanli's process (allocator, locks) but
a chain sharing its core with another. On whole cores (tag `nosmt`, run
alongside the 48 fits above) the 9355 nodes show little either way (four
chains 1.04-1.14x one, CmdStan 1.03-1.33x) and the 7513s still up to 2.5x
for stanli against 1.2x; too few seeds per node type to say more.

**Where the time goes** (`probe_bottlenecks.R`, Windows, 40 points around
the init, `results/probe_bottlenecks_win.csv`), µs per gradient:

| program | CmdStan | stanli | stanli / CmdStan |
| --- | ---: | ---: | ---: |
| `ddm4` as written | 2000 | 2071 | 1.04 |
| `ddm4`, router call dropped from the branch condition | 1556 | 1750 | 1.11 |
| `ddm4`, straight to the 4-argument `wiener_lpdf()`, no branches | 1528 | 1694 | 1.11 |
| `lnr`, multilevel, N = 1000 | 882 | 441 | 0.50 |
| `lnr`, same data, no random effects | 824 | 309 | 0.38 |
| `lnr`, multilevel, N = 10000 | 8250 | 6625 | 0.80 |

- The DDM's branches cost stanli nothing extra: with none at all it is
  still slower than CmdStan, so the gap is the `wiener_lpdf()` kernel as
  stanli runs it. The router (`cogmod_ddm_log_density_scale()`, whose value
  only decides a branch) is 22% of CmdStan's gradient and 15% of stanli's.
- Random effects cost stanli +132 µs per 1000 trials and CmdStan +58:
  brms's per-trial `r_1_1[J_1[n]] * Z_1_1[n]` terms.
- stanli's cost per trial grows with N, 0.44 to 0.66 µs from 1000 to 10000
  trials, where CmdStan's stays flat (0.88, 0.83): its advantage shrinks
  from 0.50x to 0.80x as data grow, as if its working set left the cache.

**What to suggest upstream**, from the above: (1) the in-fit gap - stanli
loses 10-50% of its gradient speed inside its own sampler where CmdStan
loses 0-10% (`multilevel.R`, `chains_probe.R`); (2) the per-trial cost
that grows with N and with a sibling thread on the core, both pointing at
memory traffic per gradient - computing the gradient of a sum over trials
block by block, as `reduce_sum()` slices it, would bound it and give
within-chain threads too, which stanli does not do today
(`threads_per_chain = 4` leaves `log_prob_grad()` and `sample_model()` on
one core for brms's per-trial loop, LNR at N = 10000 and DDM-7, Windows);
(3) the DDM kernel, 1.1x (Windows) to 1.8x
(Linux) CmdStan with no branches around it: how are stanli's kernels
built (flags, maths library)?; (4) values used only in a branch condition
could be computed without the tape, which CmdStan cannot do: 15% of the
DDM-4's gradient; (5) batched gathers for `r[J[n]]`. Not yet measured: what
stanli's runtime would gain from an x86-64-v3 build (it has no AVX
instructions at all; see "The cluster run") or vectorised `exp`/`log`/`erfc`
under `fast_math`.

`multilevel.R` notes: `grad` is a process of its own because loading
CmdStan's model methods into a process in which rstan has read a fit
segfaults on Linux (the first cluster smoke run died in `dyn.load()` in
every task); `STANLI_ML_MODELS` is colon-separated because `sbatch
--export` splits at commas; and "ms per leapfrog" is the sampling time
over the sampling leapfrogs that `nuts_params()` counts. Until 2026-10-12
evening `multilevel.R` divided warmup plus sampling time by them, which
mixes each arm's luck in warmup into the ratio; `bench.R`'s `fit` still
does, so its "ms per leapfrog" columns above carry warmup too.

```
benchmarks/stanli/hpc/run.sh push
STANLI_TAG=nosmt STANLI_ML_MODELS=lnr:lnr6:ddm4:ddm5:lba:rdm STANLI_MAX=20 \
  benchmarks/stanli/hpc/run.sh submit ml 48 --hint=nomultithread --exclude=artemis-rtx-01
STANLI_TAG=nosmt STANLI_ML_MODELS=ddm7 STANLI_TIME=4-00:00:00 \
  benchmarks/stanli/hpc/run.sh submit ml 4 --partition=long --hint=nomultithread --exclude=artemis-rtx-01
STANLI_TAG=nosmt STANLI_ML_MODELS=lnr:lba:rdm:ddm4 \
  benchmarks/stanli/hpc/run.sh submit mlchains 16 --hint=nomultithread --exclude=artemis-rtx-01
benchmarks/stanli/hpc/run.sh pull
Rscript benchmarks/stanli/hpc/summarise_ml.R --tag nosmt
Rscript benchmarks/stanli/hpc/summarise_chains.R --tag nosmt
```

## 0.20.0: chains in parallel, the cliffs gone

stanli 0.20.0 (2026-10-11) carries
[#432](https://github.com/seantalts/stanli/pull/432), "Follow-ups from
#422": chains run in parallel when given inits, and the four full-program
cliffs from `bisect.R` are fixed, both reported in our 0.19.1 reply. Also
in it: [#452](https://github.com/seantalts/stanli/pull/452) (what a brms
backend needs: completing a `brmsfit` through brms's `"mock"` backend,
partial inits; see "What the light route needs"),
[#446](https://github.com/seantalts/stanli/pull/446) (an opt-in
`fast_math` per model, not held to CmdStan's 10 ULP), and fused densities
(lognormal among them, #439 and #442). Rerun 2026-10-12: Windows from
r-universe's binary, Linux from the release tarball; `grad`, `bisect.R`
and `repro_parallel_init.R` on four tasks on one EPYC 9355 node, the fit
on 16 seeds.

`bench.R` now takes `--fast` (stanli variants also built with
`fast_math = TRUE`, engine `stanli_<v>_fast`; `orig` and `orig_s` by
default) and `fit` arms with a `_fast` suffix; `bisect.R` and
`hpc/summarise.R --kind fit` take `--tag`. The workaround arms for the
serial chains went (section "Parallel chains and inits").

**End to end, stanli is now the faster backend.** The same 2 x 2 fit as
on 0.19.0 (section "A real fit"). Linux, 16 seeds, median [10th, 90th
percentile], ratios within seed (`results/hpc/fit_0.20.0_<seed>/`,
`hpc/summarise.R --kind fit --tag 0.20.0`):

| arm | wall | wall vs CmdStan `orig` | ESS per wall second vs CmdStan `orig` | ms per leapfrog vs CmdStan `orig` |
| --- | ---: | ---: | ---: | ---: |
| cmdstanr `orig` | 86 s [71, 176] | 1 | 1 | 1 |
| cmdstanr `sel` | 135 s [113, 322] | 1.62 [1.30, 2.00] | 0.64 [0.46, 0.82] | 1.78 [1.40, 2.04] |
| stanli `orig` | 63 s [51, 163] | **0.74 [0.59, 1.02]** | **1.20 [1.04, 1.50]** | 0.85 [0.71, 1.04] |
| stanli `sel` | 45 s [40, 133] | 0.56 [0.47, 0.69] | 1.86 [1.26, 2.23] | 0.65 [0.51, 0.71] |

On 0.19.0 the stanli arms took 2.34x and 1.41x CmdStan's wall time.
Windows, seed 11, rounds 1 / 2 (`results/fit_*_0.20.0.csv`):

| arm | wall | wall vs CmdStan `orig` | ESS per wall second vs CmdStan `orig` | ms per leapfrog |
| --- | ---: | ---: | ---: | ---: |
| cmdstanr `orig` | 142 / 128 s | 1 | 1 | 13.8 / 12.6 |
| cmdstanr `sel` | 230 / 237 s | 1.62 / 1.85 | 0.66 / 0.58 | 22.8 / 23.4 |
| stanli `orig` | 61 / 58 s | **0.43 / 0.46** | **2.83 / 2.67** | 5.9 / 5.8 |
| stanli `sel` | 48 / 47 s | 0.34 / 0.37 | 2.71 / 2.49 | 4.9 / 4.7 |

The Linux run was noisy. 13 seeds ran on one EPYC 9355 node, fully
allocated at the time, `grad`'s tasks among the rest; on those alone
stanli `orig` took 0.72 [0.60, 0.88] of CmdStan's wall time and `sel`
0.53 [0.47, 0.61] (ESS per wall second 1.23 and 2.04). Seeds 14-16 first landed on a second 9355 (`rtx-01`), where
both CmdStan arms ran 5x slower than anywhere else (30-35 ms per leapfrog)
while stanli's ran normally. Not stanli: the node was nearly idle, stanli
leaves the process's CPU affinity alone (`hpc/affinity_probe.R`), and
CmdStan alone there, with no stanli model built, took 49 ms per leapfrog
(`STANLI_ARMS=cmdstanr_orig`, `results/hpc/fit_0.20.0rtx01_14/`). Not
chased further. Their rerun on an EPYC 7513 (`--exclude=artemis-rtx-01`)
put two of them on a busy node, 3-4x slow in every arm; those are what the
90th percentiles above carry.

- **The same draws as on 0.19.0.** Every arm, every seed, both platforms:
  the same leapfrog counts and ESS (on Linux even for seeds 14-16, on an
  EPYC 9334 then and a 7513 now), and on Windows the same posterior means
  to the last digit. #432 says the parallel path keeps the serial one's
  draws bitwise; it does, through all the other changes in 0.20.0. So the
  posterior and sampler-efficiency findings of the 0.19.0 fit stand
  unchanged.
- **Only the clock moved.** The chains overlap: summed chain time over wall
  is 3.6-3.9 for both stanli arms, as for cmdstanr (3.1-3.9). On 0.19.x
  it was 1.00. On Windows the stanli `orig` fit went from 195 / 169 s to
  61 / 58 s.
- **Per step stanli now pays for sharing the machine**, as CmdStan always
  did: 5.9 ms per leapfrog on Windows against 4.4-5.1 ms when its chains ran
  alone, and 0.85x CmdStan's on Linux against 0.61x on 0.19.0. That is now
  in line with `grad` (0.42x and 0.75-0.77x per gradient), as it should be.

**Per gradient nothing changed for this program.** At all 46 points of
`grad` (36 interior, 9 edge, the tail) 0.20.0 gives 0.19.1's log density
bit for bit, every engine, both platforms, and the same gradient errors
(to every digit on Windows, within 5e-15 on Linux, where CmdStan's own
entries move as much between builds on different node types), so the
exactness tables below hold as they are. Within-block ratios to CmdStan
`orig` (Windows one run; Linux four tasks, all on an EPYC 9355, against
0.19.1's one task on a 9355):

| program | Windows 0.19.1 | Windows 0.20.0 | Linux 0.19.1 | Linux 0.20.0 |
| --- | ---: | ---: | ---: | ---: |
| CmdStan `orig_s` | 0.70 | 0.68 | 0.85 | 0.82-0.84 |
| CmdStan `rw` | 0.75 | 0.75 | 0.84 | 0.84-0.86 |
| CmdStan `rw_v` | 1.05 | 1.09 | 1.02 | 1.05-1.06 |
| stanli `orig` | **0.41** | **0.42** | **0.77** | **0.75-0.77** |
| stanli `orig_s` | 0.36 | 0.35 | 0.64 | 0.63-0.64 |
| stanli `rw_v` | 0.37 | 0.39 | 0.71 | 0.72-0.73 |
| stanli `rw` | 0.32 | 0.32 | 0.60 | 0.60-0.61 |
| stanli `sel` | 0.31 | 0.31 | 0.60 | 0.58-0.60 |
| stanli `orig`, `fast_math` | | 0.41 | | 0.74-0.76 |
| stanli `orig_s`, `fast_math` | | 0.34 | | 0.62-0.64 |

The split of cogmod#5's rewrite into sigmas and likelihood comes out as
before: on CmdStan the sigmas are 0.68-0.69 (Windows) and 0.79-0.84
(Linux) and the rewritten likelihood 1.09-1.11 and 1.01-1.06; on stanli
0.83-0.86 and 0.92-0.97. Absolute times on the Linux node ran about 10%
above 0.19.1's (CmdStan `orig` 1365-1445 µs against 1275), on a node that
was fully allocated, 13 of the fit's tasks among them; the ratios did not
move.

- **`fast_math` buys nothing here**: 1.00 against the default on Windows,
  0.98-0.99 on Linux, with the same log density and gradient as the
  default at every point. #446 lists what it changes (stanc's partial
  evaluation, merging duplicate operations, fusion over shared
  parameters); the LNR's cost is elsewhere.
- **`log_sum_exp()` in place of `log_mix()`** in the package's program
  (`lse`; measured because in branch-free code `log_sum_exp()` now beats
  `log_mix()`, below; Windows, `results/*_0.20.0_lse.csv`): stanli `lse`
  takes 0.41x CmdStan `orig` against stanli `orig`'s 0.40x, and CmdStan
  `lse` 1.03x. The gain is in branch-free code only; the package's branchy
  lanes are indifferent.
- **The cliffs** are in the `bisect.R` table under "The select pattern":
  all four gone, and in branch-free code `log_sum_exp()` and a blend now
  beat `log_mix()` by 11-22%.
- `repro_issue.R` (Linux): the never-taken parameter branch costs 0.5x the
  branch-free copy (0.19.1: 0.4x).

## 0.19.1, and the rewrite from cogmod#5

stanli 0.19.1 (2026-10-04) carries [#429](https://github.com/seantalts/stanli/pull/429),
the follow-up the maintainer announced on #422: `erfc`, `log1p`, `log1m_exp`
and `inv_square` run natively inside batched loops, and density calls are
made once per batch of observations, both keeping every result bit for bit.
The R package is 0.19.0's apart from the runtime tag; r-universe still had
0.19.0 that morning, so it came from the GitHub release.

The same night the maintainer opened
[cogmod#5](https://github.com/DominiqueMakowski/cogmod/issues/5) with a
hand-written version of this benchmark's program (`lnr_rewrite.stan`,
verbatim): rows split by `dec` in transformed data, vector arithmetic per
group, every parameter branch a 0/1 blend of clamped arms (all three arms of
`cogmod_log_Phi()` kept, `log_mix()` as an elementwise `log_sum_exp()`), the
argument checks as a branch-free -inf, and the per-trial terms summed in the
original order. It claims brms's log density bit for bit on CmdStan, and
1.09 against 1.29 ms on CmdStan and 0.57 against 0.85 ms on stanli `main`
(Apple silicon).

Both claims hold, but the rewrite changes two things at once. Besides the
likelihood, it declares the two sigmas as reals: brms writes `sigmazero ~ 1`
as a 4620-vector of one intercept put through the softplus element by
element, which `benchmarks/lnr_vectorize` had already measured at a quarter
of the gradient on this laptop. So `bench.R` splits them with two more
variants (`orig_s`: brms's program with the sigmas as reals and nothing else
changed, the same log density bit for bit; `rw_v`: the rewrite with brms's
vector sigmas, what a `custom_family(loop = FALSE)` would get), and `grad`
takes `--cmdstan` and `--stanli` lists and a `--tag`
(`results/*_0.19.1.csv`; Linux in `results/hpc/grad_0.19.1_<k>/`, from
`STANLI_TAG=0.19.1 hpc/run.sh submit grad 4`, summarised by
`hpc/summarise.R --tag 0.19.1`).

µs per gradient, and the ratio to CmdStan `orig` within each block (median
over 21 blocks). Windows: one run. Linux: four tasks, one on an EPYC 9355
node, three on an EPYC 7513 node; the range is over tasks.

| program | Windows | Linux, 9355 | Linux, 7513 |
| --- | ---: | ---: | ---: |
| CmdStan `orig` | 4300 (1) | 1275 (1) | 2305-2390 (1) |
| CmdStan `rw_v` | 4450 (1.05) | 1300 (1.02) | 2205-2320 (0.92-1.00) |
| CmdStan `orig_s` | 2950 (0.70) | 1085 (0.85) | 1920-1955 (0.80-0.84) |
| CmdStan `rw` | 3150 (0.75) | 1070 (0.84) | 1800-1865 (0.75-0.80) |
| stanli `orig` | 1750 (0.41) | 980 (0.77) | 1705-1730 (0.71-0.74) |
| stanli `rw_v` | 1550 (0.37) | 915 (0.71) | 1540-1605 (0.64-0.68) |
| stanli `orig_s` | 1600 (0.36) | 820 (0.64) | 1445-1465 (0.61-0.63) |
| stanli `rw` | 1350 (0.32) | 770 (0.60) | 1300-1370 (0.54-0.59) |
| stanli `sel` | 1300 (0.31) | 765 (0.60) | 1335-1365 (0.56-0.58) |

The two changes, each measured both ways round (within-block ratios;
Linux over the four tasks):

| change | engine | Windows | Linux |
| --- | --- | ---: | ---: |
| sigmas as reals | CmdStan | 0.70-0.71 | 0.80-0.85 |
| | stanli | 0.82-0.87 | 0.84-0.85 |
| the rewritten likelihood | CmdStan | **1.05-1.10** | **0.92-1.02** |
| | stanli | 0.83-0.89 | 0.90-0.94 |
| both | CmdStan | 0.75 | 0.75-0.84 |
| | stanli | 0.75 | 0.76-0.79 |

- **On CmdStan the gain is the sigmas.** The rewritten likelihood is a
  little dearer on Windows and at most 8% cheaper on Linux, which is where
  `benchmarks/lnr_vectorize` left vectorising (about 0.98 on the vignette
  model) and below the 10% that README set for reopening it. Not profiled;
  the blends take two logarithms per log Phi where the branch takes one,
  which would weigh most where libm is slow.
- **On stanli both help**, the likelihood by 6-17%: whole-vector operations,
  which the per-observation lanes of #424 do not match. Hand-written `sel`
  reaches the same speed (it has the sigmas as vectors, but a likelihood
  stanli runs even faster).
- **stanli 0.19.1 on the unmodified program**: 0.41x CmdStan on Windows
  against 0.48x on 0.19.0, 0.71-0.77x on Linux against 0.81-0.86x (different
  nodes, so only roughly comparable).

Exactness, against CmdStan `orig` at the init, 5 points near it, 30 further
out (sd 0.5 on the unconstrained scale) and 9 edge points: `poutlier` 1e-13,
exactly 0 and exactly 1; `ndt` 1.5x (more responses below it), the same with
`poutlier` = 0 (log density -inf), `ndt` = 2.5 ms; both sigmas at softplus(-3)
and softplus(3); the two nus at +-3. Plus the far-tail point above.

Interior points, Windows then Linux (where the four tasks ran the same 36
points, so the counts are out of 144 evaluations):

| engine | log density identical | largest ULP apart | gradient, largest relative error |
| --- | ---: | ---: | ---: |
| CmdStan `orig_s` | all | 0 | 1.7e-14 |
| CmdStan `rw`, `rw_v` | all | 0 | 2.7e-13 |
| stanli `orig` | 25/36, 136/144 | 4, 1 | 1.6e-14, 2.1e-15 |
| stanli `rw`, `rw_v` | 21/36, 92/144 | 6, 4 | 1.7e-13 |
| stanli `sel` | 3/36, 8/144 | 22, 28 | 1.1e-13 |

On CmdStan, `rw` gives brms's log density bit for bit at every point on
both platforms, edge and tail points included, as claimed. Its gradient is
not bitwise (summation order), and differs in finiteness at exactly the
places cogmod#5 lists: at `poutlier` = 1 it is non-finite where brms's is
finite, and where the log density is -inf it is finite (brms's is not).

## History: 0.18.1, #422, and the fixes

On 0.18.1 stanli refused the program brms writes over one call: `log_mix()`
inside an `if` on a parameter. Comparisons used as values (`log(y > mu)`)
were refused too. With `log_mix()` written out (`lse` below) the package's
density compiled and was exact but cost 6.6x CmdStan per gradient: every
parameter-dependent `if` became an interpreted "runtime-control region", and
a never-taken one cost 2.7x on its own. Rewritten without control flow
(`bf`, `sel`), the same density cost 0.23-0.30x, and a 4-chain fit was 3.7x
faster per chain than CmdStan.

We reported it as #422 on 2026-10-02. Within a day the maintainer merged
[#423](https://github.com/seantalts/stanli/pull/423) (parameter-dependent
branches in per-observation functions were O(N²) per gradient; `log_mix()`
and comparisons lowered inside regions) and
[#424](https://github.com/seantalts/stanli/pull/424) (a loop whose body
branches on a parameter is compiled once into one register program and run
over 64-observation lanes; #423 alone had made the ex-Gaussian, GEG and LNR
about 10x slower under sampling).
#424 also added all 22 cogmod families, generated with brms, to stanli's own
test corpus (`tests/cogmod`, with a note and our MIT licence). Both shipped
in 0.19.0. The maintainer's own numbers (Apple silicon, one chain of
50 + 50, `notes/performance/2026-10-03-cogmod-sampling.md` in their repo):

| family | stanli 0.18.1 | stanli `main` | CmdStan | `main` / CmdStan |
| --- | ---: | ---: | ---: | ---: |
| gamma | 275 µs | 138 µs | 197 µs | 0.70 |
| weibull | 269 µs | 124 µs | 193 µs | 0.64 |
| exgaussian | 838 µs | 270 µs | 412 µs | 0.66 |
| geg | 1469 µs | 485 µs | 639 µs | 0.76 |
| loggamma | 1864 µs | 313 µs | 454 µs | 0.69 |
| lognormal | 17550 µs | 551 µs | 755 µs | 0.73 |
| exwald | 26629 µs | 693 µs | 1586 µs | 0.44 |
| choco | 2946 µs | 666 µs | 812 µs | 0.82 |
| betagate | 1394 µs | 655 µs | 597 µs | 1.10 |
| lba1 | over 5 min | 707 µs | 906 µs | 0.78 |
| lba2 | over 5 min | 1241 µs | 1423 µs | 0.87 |
| rdm | over 5 min | 1698 µs | 2425 µs | 0.70 |
| lnr_bench (this folder's LNR) | does not compile | 1216 µs | 1265 µs | 0.96 |
| lnr | does not compile | 1977 µs | 2269 µs | 0.87 |

On 0.18.1 the DDM failed on the six-argument `wiener_lpdf()` and the inverse
Gaussian on a guard return; both compile now. `lnr_bench` peaks at 138 MB
resident with lanes.

The maintainer's reply on #422 (2026-10-03) settles the select question.
stanli will not if-convert branches into selects. It is held to CmdStan's
numerics, so that a model gives the same answer under both, and selects do
not. And when they tried it, computing both arms was slower than the
per-observation branch it already takes. `sel` will likely stay faster
because a loop with no parameter branch runs as whole-vector operations,
which lanes cannot match. Their next work keeps every result bit for bit:
computing `erfc`, `log1p` and the like inside the batched loop, and one
density call per batch rather than per observation. They also listed where
`sel` differs from the package's code (see "The select pattern").

## Variants

| name | what |
| --- | --- |
| `orig` | the program `brm()` writes, unchanged |
| `lse` | `orig` with `log_mix(poutlier, lp_out, lp_dec)` written out as `log_sum_exp(log(poutlier) + lp_out, log1m(poutlier) + lp_dec)`; nothing else changes |
| `bf` | `orig`'s functions block replaced by `lnr_branchfree_functions.stan`: `log Phi` as `log(0.5 * erfc())` with no tail series, winner and loser picked by `dec` (data), `fmax(Y - ndt, 1e-300)` for the `t <= ndt` return. `sigmabias = 0` only |
| `sel` | `orig`'s functions block replaced by `lnr_select_functions.stan`: as `bf`, but `log Phi` keeps the series below -25 as a select, and `t <= ndt` is a mask inside `log_mix()`. No checks. `sigmabias = 0` only |
| `orig_s` | `orig` with `sigmazero` and `sigmaone` declared `real` in the model block: one softplus each instead of 4620. The same log density, bit for bit |
| `rw` | cogmod#5's program (`lnr_rewrite.stan`): its own transformed data and model block, sigmas as reals. `sigmabias = 0` and this formula only |
| `rw_v` | `rw` with the sigmas as brms writes them, `vector[N]` through the softplus, sliced per group |

All keep the parameters, priors and generated quantities, so `rename_pars()`
and the rest of brms see the program they expect.

## Exact

At the init and five points around it (`results/check.csv`,
`results/hpc/grad_*/check.csv`): stanli's log density for `orig` equals
CmdStan's within 1e-11 on Windows (values of -645 to -3547) and exactly on
Linux.
Gradients, largest relative difference from CmdStan `orig`'s:

| variant | Windows | Linux |
| --- | ---: | ---: |
| stanli `orig` | 5.2e-16 | 4.2e-16 |
| stanli `sel` | 1.8e-15 | 2.0e-15 |
| stanli `bf` | 9.5e-15 | 9.5e-15 |
| stanli `lse` | 1.4e-14 | 1.4e-14 |
| CmdStan `sel` | 7.7e-16 | 7.7e-16 |

The tail (`results/tail.csv`): both sigmas at softplus(-3) = 0.049 and one
trial replaced by a 12 s response, which puts the loser's survival at
x = -65. `orig` is finite there, by design. In stanli 0.19.0, every variant
gives `orig`'s log density; the gradients:

| variant | Windows | Linux |
| --- | --- | --- |
| `orig` | finite, 1.0e-15 | finite, 2.6e-16 |
| `sel` | finite, 1.4e-15 | finite, 1.3e-15 |
| `lse` | finite, 1.0e-14 | finite, 1.0e-14 |
| `bf` | **not finite** | **not finite** |

`bf` is the 0.3.3 LNR bug: the value survives, the gradient does not.

## Cost per gradient

R loops over 200 points, a different one per call, alternating the programs
block by block (21 blocks), one core each. Windows: one run
(`results/time.csv`). Linux: four tasks, all on an EPYC 9334 node
(`results/hpc/grad_*/time.csv`); ms is the median over tasks, the ratio's
range is over tasks.

| program | Windows 0.18.1 | Windows 0.19.0 | Linux 0.19.0 |
| --- | ---: | ---: | ---: |
| CmdStan `orig` | 4.20 ms (1) | 4.70 ms (1) | 2.07 ms (1) |
| CmdStan `bf` | 4.05 (0.96) | 4.70 (1.00) | 1.93 (0.90-0.94) |
| CmdStan `sel` | 5.70 (1.36) | 6.40 (1.36) | 3.99 (1.89-1.96) |
| stanli `orig` | refused | 2.25 (**0.48**) | 1.72 (**0.81-0.86**) |
| stanli `lse` | 27.8 (6.6) | 2.30 (0.49) | 1.81 (0.86-0.89) |
| stanli `bf` | 0.95 (0.23) | 1.15 (0.24) | 0.84 (0.39-0.41) |
| stanli `sel` | 1.25 (0.30) | 1.40 (0.30) | 1.08 (0.50-0.53) |

On Windows the block-by-block ratio of stanli `orig` to CmdStan `orig`
ranged 0.45-0.52 (10th-90th percentile). Absolute times on the laptop move
10-20% between runs (CmdStan `orig`: 4.20 ms, then 4.70); the ratios within a
run are what to compare.

- **The platform matters more to CmdStan than to stanli.** This laptop's
  maths library is slow (`exp` at 33 ns, `erfc` at 48 ns:
  `benchmarks/lnr_vectorize/README.md`), and CmdStan is bound by it: 2.3x
  faster on Linux. stanli ships its own prebuilt kernels and gains 1.3x. So
  stanli's lead is 2x on Windows and 1.2x on Linux.
- **`sel` costs CmdStan more on Linux** (1.9x against 1.36x): its selects
  evaluate every arm, and where `erfc` is cheap the extra series weighs more.
  Measured, not profiled.
- **`log_sum_exp()` is no longer a cliff inside the branchy program**
  (`lse` = `orig`), but it still is in branch-free code (`bisect.R`, below).

`grad` carries on without CmdStan `bf` if that executable goes missing:
Windows Defender quarantined `lnr_bf.exe` straight after compiling, twice on
2026-10-03, a false positive it also makes on other CmdStan executables.

## The select pattern

A branch `if (x < -25) lo(x) else hi(x)` becomes

```stan
real xl = fmin(x, -25);       // lo's arm, clamped into its domain
real xh = fmax(x, -25);       // hi's arm, likewise
real w = 1 - step(x + 25);    // 1 strictly below -25
return w * lo(xl) + (1 - w) * hi(xh);
```

Both arms are evaluated, but each on an input where it is finite, so the
arm not picked has a finite value and a weight of 0. Its input is the clamp
edge, a constant, so it passes back an exact zero, never `0 * inf`. The
picked arm gets the whole gradient.

Stan's `step(0)` is 1, so the weight has to be written to put the boundary
on the same side as the branch's strict `<`: `step(-25 - x)` would pick the
series at x = -25 itself, where the two arms agree in value but their
derivatives differ by 186 ULP. Likewise the `t <= ndt` mask is
`lp_dec - step(ndt - Y) * 1e300`, which drops the decision component at
Y = ndt as `orig` does; `(step(Y - ndt) - 1) * 1e300` kept it there. Both
were wrong in the first version and fixed on 2026-10-03, after the
maintainer pointed them out; the timings did not move.

Where `sel` still differs from the package's code, all confirmed:

- **Above x = 0**, `log(0.5 * erfc(-x / sqrt2))` in place of
  `log1p(-0.5 * erfc(x / sqrt2))`: the same gradient and an absolute error
  below 1e-16, but not the same value bit for bit: 41-49 ULP at x = 3, and 0
  in place of -1.13e-19 at x = 9. A third arm would fix it at the price of a
  second `erfc`.
- **At `poutlier = 0` with Y below ndt**, `orig` returns -inf; the mask gives
  -1e300 and a NaN gradient (`log_mix()`'s partials take 0 * inf). Either
  way the proposal is rejected.
- **No checks.** `dec` in {0, 1} and Y > 0 are guaranteed only by
  `.cogmod_checkdata()`, which runs from `cogmod_priors()` alone: brms
  declares Y with no lower bound. And a softplus link underflows to exactly
  0 below about -745: during warmup of the cmdstanr `sel` fit,
  `lognormal_lpdf()` rejected a zero scale six times per run, with a message,
  where `orig` returns -inf silently.

None of it matters for fitting this model; it is why a compiler that must
keep a program's meaning cannot write `sel`.

What it took for stanli to keep the select version fast (`bisect.R`, µs per
gradient; every row has the same log density):

| `log Phi` | `t <= ndt` | mixture | checks | 0.18.1 Win | 0.19.0 Win | 0.19.0 Linux | 0.19.1 Linux | 0.20.0 Win | 0.20.0 Linux |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: | ---: | ---: |
| *`orig`, unmodified* | | | | refused | 1700 | 1680 | 1760 | 1400 | 1010 |
| erfc only (`bf`) | clamp to 1e-300 | `log_mix()` | none | 900 | 900 | 825 | 1070 | 900 | 640 |
| select | clamp to 1e-300 | `log_mix()` | none | 1000 | 1100 | 1050 | 1360 | 1000 | 770 |
| select | mask | `log_mix()` | none | 1100 | 1100 | 1060 | 1370 | 1100 | 790 |
| select | mask | `log_mix()` | data only (`dec`, `Y > 0`) | 4800 | 5000 | 5070 | 6370 | **1100** | **790** |
| erfc only | clamp | `log_sum_exp()` + `log1m()` | none | 4700 | 4800 | 4075 | 5770 | **700** | **510** |
| erfc only | blend: `w * log_mix(...) + (1 - w) * (...)` | `log_mix()` | none | 5500 | 6200 | 5060 | 6470 | **800** | **570** |
| select | blend | `log_mix()` | none | 10900 | 11600 | 9580 | 14880 | **1100** | **790** |

0.19.0 Linux is the median over four tasks on an EPYC 9334 node; 0.19.1
Linux over the three on an EPYC 7513 (slower: the fourth task, on a 9355,
read 0.56-0.6x these throughout); 0.20.0 Linux over four tasks on an EPYC
9355. Against 0.19.1's task on the same kind of node (990, 640, 760, 780,
3280, 2620, 3860, 7450) only the four bold rows moved. On 0.19.0 the unmodified program costs
1.6x the select version, and the cliffs in branch-free code are all still
there: a data-only early return, `log_sum_exp()`, or a blend costs 4-6x.
0.19.1 leaves them as they were: against the same row without the
construct, the data check costs 4.7x (0.19.0: 4.8x), `log_sum_exp()` 5.4x
(4.9x), the blend 6.0x (6.1x), select plus blend 10.9x (9.1x).
`orig` has the same data checks and costs 1.7 ms, so it is not the check
itself but the check in an otherwise branch-free function. #424 says a
branch on data alone "keeps the old lowering", which fits. Alone, in a
minimal program, `log_sum_exp()` and the data check cost 0.9-1.2x
(`repro_issue.R`, part 4), so the cliffs need the rest of the program.
`bisect.R` writes each full program to `DIR/bisect_<variant>.stan`.

**0.20.0 removes all four** ([stanli#432](https://github.com/seantalts/stanli/pull/432),
with the two causes: stanc marks a comparison of two data values as
autodiffable, so a data-only return made the loop count as
parameter-controlled; and lanes reading one parameter in more than one
place, as `log(p)` and `log1m(p)` do in every observation, were never fused
into vector kernels). Against the same row without the construct, Windows
then Linux: the data check costs 1.00x and 1.00x, select plus blend 1.10x
and 1.03x, and `log_sum_exp()` and the blend are now *cheaper* than
`log_mix()`: 0.78x and 0.80x, 0.89x and 0.89x. #432 measured 3.4-5.7x
on the same variants against `main`; on our Linux node type it is
4.2x (data check), 5.1x (`log_sum_exp()`), 6.8x (blend) and 9.4x (select
plus blend). The rows without the constructs, `orig` included, are where
0.19.1 left them.

`repro_issue.R` on 0.19.0 (`results/repro_issue.txt`,
`results/hpc/grad_*/repro_issue.txt`): every construct compiles, and the
never-taken parameter branch costs 0.7x (Windows) and 0.9x (Linux) the
branch-free copy, against 2.7x on 0.18.1. On 0.19.1 (Linux,
`results/hpc/grad_0.19.1_*/repro_issue.txt`) it costs 0.4x: the branchy
copy now runs faster than the branch-free one.

## A real fit: {orig, sel} x {cmdstanr, stanli}

On stanli 0.19.0, whose chains ran one at a time (next section); the
0.20.0 rerun, with the same draws and the chains in parallel, is in the
section "0.20.0". 4 chains x (500 warmup + 500 sampling), in parallel:

- `cmdstanr_orig`: `brm(backend = "cmdstanr")`, the standard route;
- `cmdstanr_sel`: `sel` compiled by cmdstanr, read back with
  `brms:::read_csv_as_stanfit()`;
- `stanli_orig`, `stanli_sel`: the light route.

Every arm starts its chains from the same `cogmod_inits()` and seed; stanli
takes them through `unconstrain()` as a chains x parameters matrix.
Compiling and building happen before the clock starts. Arms run in turn,
the first one rotating with the seed. Windows: seed 11, two rounds
(`results/fit_summary.csv`). Linux: 16 seeds, one round each, 12 on an EPYC
9355 node and 4 on an EPYC 9334 (`results/hpc/fit_*/`).

Linux, median [10th, 90th percentile] over the 16 seeds, ratios within
seed:

| arm | slowest chain | ms per leapfrog | per leapfrog vs CmdStan `orig` | ESS per 1000 leapfrog | min bulk ESS | wall | wall vs CmdStan `orig` |
| --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| cmdstanr `orig` | 63 s | 6.08 | 1 | 22.6 | 874 | 65 s [54, 90] | 1 |
| cmdstanr `sel` | 102 s | 9.85 | 1.73 [1.36, 1.96] | 22.7 | 891 | 103 s [92, 155] | 1.71 |
| stanli `orig` | 38 s | 3.60 | 0.61 [0.51, 0.70] | 22.3 | 897 | 140 s [133, 190] | **2.34** |
| stanli `sel` | 23 s | 2.27 | 0.40 [0.31, 0.42] | 23.2 | 886 | 89 s [83, 119] | 1.41 |

Windows, seed 11, rounds 1 and 2 (the draws repeat; only timings differ):

| arm | slowest chain | ms per leapfrog | ESS per 1000 leapfrog | wall |
| --- | ---: | ---: | ---: | ---: |
| cmdstanr `orig` | 127 / 149 s | 12.5 / 14.8 | 21.6 | 130 / 151 s |
| cmdstanr `sel` | 218 / 219 s | 21.7 / 21.9 | 22.4 | 220 / 221 s |
| stanli `orig` | 60 / 44 s | 5.1 / 4.4 | 25.7 | 195 / 169 s |
| stanli `sel` | 36 / 42 s | 3.5 / 4.2 | 19.9 | 130 / 157 s |

- **The same posterior.** Posterior means differ from cmdstanr `orig`'s by at
  most 0.13 posterior SDs over the 16 seeds and 0.06 on Windows; R-hat at
  most 1.010, no divergences, every bulk ESS above 600. `loo_compare()` on
  Windows: elpd differences 0.2 or less (SE 0.1).
- **The same sampler efficiency.** ESS per leapfrog step and leapfrog counts
  match across all four arms, so stanli's NUTS and adaptation are as
  efficient here as CmdStan's. What differs is the cost of a step.
- **Per step, trust `grad` over the fit.** In the fits stanli's chains ran
  one at a time (next section) while CmdStan's four ran at once, so stanli's
  chains had the machine to themselves. That flatters it most on the laptop:
  in a toy test four concurrent chains each ran about 2x slower than one
  alone. Hence 0.30-0.41x per step in the Windows fit against 0.48x per
  gradient, and 0.61x against 0.84x on Linux.
- **The wall time is the sum of stanli's chains.** On Linux, 4 x the mean
  chain over the wall is 1.00 for both stanli arms and 3.7 for both cmdstanr
  arms. That is the serial-chains bug below, not post-processing:
  `as_stanfit()` and `rename_pars()` take under a second. Fixed in 0.20.0.

## Parallel chains and inits (stanli 0.19.0 and 0.19.1; fixed in 0.20.0)

`sample_model()` ran the chains one after another whenever it got `init`,
a shared vector or a chains x parameters matrix, whatever `parallel_chains`
said. Without `init` the same model ran them in parallel.
`repro_parallel_init.R`, a toy with a `log_mix()` likelihood, N = 5000, on
Windows, 0.19.0 (0.19.1 the same, 2026-10-04: wall/sum 0.35, 1.00, 1.00):

```
init none   wall  2.49 s, chains max  2.48, sum  7.59: wall/sum 0.33
init vector wall  8.11 s, chains max  2.59, sum  8.09: wall/sum 1.00
init matrix wall  8.73 s, chains max  3.42, sum  8.70: wall/sum 1.00
```

It was not the model: a branch, a user function, `log_mix()`,
`lognormal_lpdf()`, a data ternary all ran in parallel without inits. The
cause was native, as `sample_model()` suggested by passing the inits and
`parallel_chains` through unchanged: the C API handled supplied inits with
a serial loop of its own
([stanli#432](https://github.com/seantalts/stanli/pull/432)). On 0.20.0
the same script, the same laptop (2026-10-12):

```
init none   wall  2.48 s, chains max  2.47, sum  7.44: wall/sum 0.33
init vector wall  2.06 s, chains max  2.05, sum  6.39: wall/sum 0.32
init matrix wall  2.58 s, chains max  2.57, sum  6.80: wall/sum 0.38
```

#432 says the draws are bitwise those of the serial loop. The 0.18.1 fit
started from stanli's random inits, so whether 0.18.1 had the same
behaviour is unknown.

Two workaround arms were written for it and never produced a result:
one R process per chain (`stanli_orig_proc`, `stanli_sel_proc`) and
stanli's own random inits (`stanli_orig_rand`). On the cluster
(`run.sh submit fitpar 16`, 2026-10-03 22:37) all 16 tasks exited 132
within 13-22 s on EPYC 7513 nodes: the AVX-512 bridge of the cluster
section, not the workaround. On Windows, seed 11, the chain at stanli
seed 1101 from `cogmod_inits()` row 1 went onto a slow-accumulator ridge
(`b_nuone_Intercept` about -2.4 +- 2.4, R-hat 1.6 in a short smoke run)
under both `orig` and `sel`, with the same leapfrog counts in the other
three chains of each, so deterministic, not a hang: ~30000 leapfrog steps
in the smoke run, and the full run was stopped after 15 minutes on that one
chain (2026-10-04). The plain stanli arms at seed 11 (stanli seed 11, four
chains) never met that ridge. Both arms, and the `fitpar` cluster kind,
were removed with the fix (2026-10-12).

## The cluster run

`hpc/run.sh` drives Artemis over SSH (the VPN must be up), modelled on
`benchmarks/lba_screen/`:

```
benchmarks/stanli/hpc/run.sh push          # scripts + minimal package tree + stanli tarball
benchmarks/stanli/hpc/run.sh install       # brms e71e9d7, loo, posterior, stanli + runtime
benchmarks/stanli/hpc/run.sh submit fit 16 # one seed per task, 4 CPUs
benchmarks/stanli/hpc/run.sh submit grad 4 # grad + bisect + repros per task, 2 CPUs
benchmarks/stanli/hpc/run.sh pull          # -> results/hpc/
Rscript benchmarks/stanli/hpc/summarise.R

# a rerun on a new stanli, kept apart from the old results:
STANLI_TAG=0.20.0 benchmarks/stanli/hpc/run.sh submit grad 4
STANLI_TAG=0.20.0 benchmarks/stanli/hpc/run.sh submit fit 16
Rscript benchmarks/stanli/hpc/summarise.R --tag 0.20.0
Rscript benchmarks/stanli/hpc/summarise.R --kind fit --tag 0.20.0

# rerun some indices (a range, not a count), here away from one node:
STANLI_TAG=0.20.0 benchmarks/stanli/hpc/run.sh submit fit 14-16 --exclude=artemis-rtx-01
```

Everything goes into the experiment's own directories and R library
(`/mnt/lustre/users/psych/dmm56/cogmod_stanli`): brms from GitHub at
`e71e9d7` with the loo (2.10.1) and posterior (1.7.0) it needs, and stanli
(`STANLI_VERSION`, 0.20.0 since 2026-10-12), ahead of the production library
on `R_LIBS`. The production library
is never written to. Each task works in its own copy of the tree on scratch.
Wall time on a shared node is noise, so the fit's comparisons are within a
seed, and ESS per leapfrog step needs no clock.

**The SIGILL (2026-10-03), diagnosed 2026-10-04.** `artemis-general-02`, a
VM whose CPU reads "AMD EPYC-Genoa Processor", killed R with SIGILL (exit
132, R's handler: `cause 'illegal operand'`) in `stanli_model()` on all
three `fit` seeds it got, and that evening all 16 `fitpar` tasks died the
same way on the bare-metal EPYC 7513 nodes, while the EPYC 9334 and 9355
nodes never faulted. Seeds 10-12 were rerun with
`--exclude=artemis-general-02`. The maintainer asked for the VM's CPU
flags and `R CMD config CFLAGS`, and the second question was the answer:
the R module's Makeconf has `-march=native`, the library is on Lustre and
shared by every node, and `sacct` shows the 0.19.0 install job ran on
`artemis-rtx-01`, an EPYC 9355 (znver4). stanli's R bridge, `stanli.so`,
came out with one AVX-512F instruction, `vcvttsd2usi %xmm0,%r14d` at file
offset 0x2e2f, the `(uint32_t) Rf_asReal(seed)` cast in
`stanli_r_model_new`; every fault address R printed ends in `e2f`. The VM
and the 7513s have AVX2 but no AVX-512 (x86-64-v3; the 9355 is v4). The
0.19.1 reinstall landed on a 7513, so the current bridge is portable, which
is why both runtimes (the v0.19.0 asset re-downloaded, same bytes) build
and sample the LNR on all three node types today, and why `grad` under
0.19.1 ran on a 7513 without a fault. Rebuilding the bridge on `rtx-01`
and running `hpc/sigill_probe.R` against it on the VM or a 7513 reproduces
exit 132 at the same offset on demand; the same build on the VM has no
`zmm`/EVEX instruction. The runtime `.so` (both releases) has no AVX-512
and no AVX at all, so it was never the cause. `dmesg` on those nodes
carries no `traps:` line despite `debug.exception-trace = 1`, and the
journal is not readable, so R's address was the evidence. Fix in
`run.sh install`: a `Makevars.portable` with `-march=x86-64-v3` through
`R_MAKEVARS_USER`; `--exclude` is no longer needed. stanli's Linux runtime
needs glibc 2.28 at most; the nodes have 2.34.

## What the light route needs

- `brms:::rename_pars()` is internal. After it, `summary()`, `fixef()`,
  `loo()` and `loo_compare()` worked on the result. The comment it came from
  says refit-based functions (`reloo`, moment matching) do not.
- **A public route since stanli 0.20.0**
  ([stanli#452](https://github.com/seantalts/stanli/pull/452), their
  `docs/teaching-support.md`): brms's own `"mock"` backend, which takes a
  finished `stanfit`:

  ```r
  fit <- brm(formula, data, prior = prior, stanvars = stanvars,
             backend = "mock", mock_fit = stanli::as_stanfit(stanli_fit))
  ```

  On the LNR (2026-10-12, 4 x 200 draws) it gives the same `brmsfit` as
  `rename_pars()`: identical `fixef()`, variables and `loo()`, and
  `posterior_predict()` runs on it; 0.33 s against 0.05 s. `bench.R`
  keeps `rename_pars()` so that its wall times stay comparable with the
  0.19.0 runs. Methods that sample again (`update()`, `kfold()`) do not
  work either way; stanli says a `backend = "stanli"` for brms exists as a
  branch, not opened upstream.
- **What works on the result** (`postprocess.R`, 2026-10-12, stanli
  0.20.0, each tool also run on a cmdstanr fit of the same model): every
  tool that reads the draws. `loo()`, `loo_compare()`, `waic()`,
  `log_lik()`, `posterior_predict()`, `posterior_linpred()`,
  `hypothesis()`, `nuts_params()`/`rhat()`, `priorsense`,
  `bayestestR::describe_posterior()`, `emmeans`, and `modelbased` and
  `marginaleffects` (on the LNR on the link scale; on an RT-only lognormal
  also `estimate_contrasts()`, `avg_comparisons()`, `conditional_effects()`,
  `pp_check()` and `bayes_R2()`), with the same numbers as cmdstanr's to
  Monte Carlo error; `loo()` works after `saveRDS()`/`readRDS()` too. The
  LNR's `posterior_epred()`, `conditional_effects()`, `pp_check()` and
  response-scale marginal means fail on both backends: cogmod gives a race
  model no expectation (numerical integration per draw and observation),
  and `posterior_predict()` returns two columns per trial (RT, choice).
  What fails on stanli alone is what needs the model's log density at new
  parameter values: `loo(moment_match = TRUE)` and `bridge_sampler()`
  (both recompile with rstan, then fail), so no Bayes factors, and
  `update()` ("'compiled_model' is not supported in the mock backend"),
  so no `kfold()` or `reloo`. An unofficial wrapper posted on brms#1911,
  [brm_stanli](https://github.com/StaffanBetner/brm_stanli)
  ([comment](https://github.com/paul-buerkner/brms/pull/1911#issuecomment-5949749007)),
  adds moment matching without an rstan recompilation and `update()`; not
  tried here, since an official backend would replace all of this.
- `cogmod_inits()` goes in through `unconstrain()` and
  `sample_model(init = <matrix>)`, or `sample_cstan(init = <list of lists>)`.
  On 0.19.x the first ran the chains serially (above); 0.20.0 runs them in
  parallel. 0.20.0 also accepts partial inits (`unconstrain(partial =
  TRUE)`, partial lists in `sample_cstan()`), as CmdStan does; the
  second route was not tried.
- Reading cmdstanr's output into a stanfit, as the `cmdstanr_sel` arm does,
  needs `brms:::read_csv_as_stanfit()`: `rstan::read_stan_csv()` cannot
  parse CmdStan 2.38's CSV header ("object 'n_kept' not found").
- brms#1911 itself may change all this. Paul Bürkner's reply there plans to
  stop converting fits to `stanfit`, so the post-processing would be per
  backend.

## For cogmod

- **Nothing in the package needs to change.** Its own Stan code runs in
  stanli 0.20.0 exactly and faster than in CmdStan per gradient, 2.4x on
  Windows and 1.3x on Linux, with the same sampler efficiency, and a fit
  from `cogmod_inits()` now finishes in 0.74x CmdStan's wall time on Linux
  and under half of it on Windows.
- **No select source.** `sel` buys another 1.3x in stanli, but costs CmdStan
  1.36-1.9x, so it would have to be a second, stanli-only copy of every
  family's Stan code, with its own gradient checks, and it is not quite the
  package's numerics. The maintainer has ruled out producing it
  automatically. Not worth it.
- **No vectorised LNR either.** cogmod#5's likelihood is exact, but on
  CmdStan it buys nothing on Windows and 0-8% on Linux, and as a
  `custom_family(loop = FALSE)` it would bring every cost
  `benchmarks/lnr_vectorize` listed (a signature per formula, no
  `weights()`/`cens()`/`trunc()`/mixtures, misaligned `dec` under
  threading). Its 6-17% in stanli does not pay for that either.
- **What does pay is the sigmas.** Writing an intercept-only dpar as a real
  instead of a vector of N copies is 15-30% per gradient on CmdStan and
  13-18% on stanli, with the same log density bit for bit. brms writes the
  vector, so it is brms's to change; until then, leaving an intercept-only
  dpar out of `bf()` gets the real (on its natural scale, which
  `cogmod_priors()` and `cogmod_inits()` handle; its sampling efficiency
  against the linked intercept is not measured).
- **What decides it is upstream, and it is now brms's side.** stanli has
  fixed both things we raised (chains in parallel with inits, the cliffs in
  branch-free code), and since 0.20.0 a stanli fit becomes a `brmsfit`
  through brms's own `"mock"` backend. What is left is whether brms merges
  a backend (#1911; stanli keeps a `backend = "stanli"` branch of brms,
  not opened upstream), without which `update()`, `kfold()` and moment
  matching cannot refit.
- Only the LNR was measured here. stanli's own corpus has all 22 families
  compiling and matching CmdStan's gradients, most of them faster.

## Upstream

[seantalts/stanli#422](https://github.com/seantalts/stanli/issues/422)
(2026-10-02), with `repro_issue.R` as its reproducer: the `log_mix()` gap,
comparisons as values, the cost of a never-taken branch, and a request to
if-convert small pure branches. Answered by #423 and #424. Our follow-up
(2026-10-03 morning) thanked them for #423, withdrew a point about per-chain
inits (our misuse: `sample_model()` takes unconstrained values, not lists),
and asked which families failed to compile. Their reply (2026-10-03 evening)
announced 0.19.0, declined if-conversion with reasons, listed where `sel`
differs (all confirmed, two fixed), and asked for the variants behind the
full-program cliffs. Our reply with the 0.19.1 results
([comment](https://github.com/seantalts/stanli/issues/422#issuecomment-5978573943),
2026-10-04): the exactness, cost and sampling figures above, the serial
chains with `init` (reproducer `repro_parallel_init.R`, offered as its own
issue), the `bisect.R` variants they asked for, the SIGILL node, and yes to
keeping the families in their corpus. The maintainer's two follow-ups the
same day asked for the VM's CPU flags, the runtime it loaded, the kernel's
trap line and `R CMD config CFLAGS`, and noted the runtime scans as
baseline x86-64. Our diagnosis
([comment](https://github.com/seantalts/stanli/issues/422#issuecomment-5980192431),
2026-10-04 afternoon; cluster section): the bridge compiled with
`-march=native` on an AVX-512 node, so nothing for stanli to fix.

2026-10-04: 0.19.1 released with #429, announced on #422, and
[cogmod#5](https://github.com/DominiqueMakowski/cogmod/issues/5) opened with
the rewrite (section above). Answered the same morning
([comment](https://github.com/DominiqueMakowski/cogmod/issues/5#issuecomment-5978542829)):
exact as claimed; on CmdStan the gain is the sigmas (0.70-0.85), the
likelihood rewrite 0.92-1.10, the two multiplicative; the sigmas sent to
brms (below), the looped likelihood kept.

2026-10-04, brms: the scalar intercept-only dpars went upstream as
[brms#1945](https://github.com/paul-buerkner/brms/issues/1945) (issue: the
0.41x and 0.70-0.85x figures above, `scalar_dpar.R` as the script) and
[brms#1946](https://github.com/paul-buerkner/brms/pull/1946) (PR, see
Status). Both drafts deleted once posted, as with the others. On
2026-10-09 Paul Bürkner asked why anyone would write a formula for a
parameter they do not predict, and whether the speed was worth the special
case. Our answer
([comment](https://github.com/paul-buerkner/brms/issues/1945#issuecomment-6079343933)):
`sigma ~ 1` and no formula are not the same model (link scale and prior,
so also what a `sigma ~ condition` against `sigma ~ 1` comparison holds
fixed), `~ 1` is a placeholder in iterative work and a default in wrapper
packages, and nothing tells users it costs anything. He will review the PR.

2026-10-11: stanli 0.20.0, with
[#432](https://github.com/seantalts/stanli/pull/432) (merged 2026-10-04,
"Follow-ups from #422"): the chains run in parallel with inits, and the four
`bisect.R` cliffs are fixed, both reported in our 0.19.1 reply. Also
[#452](https://github.com/seantalts/stanli/pull/452) (brms support: the
`"mock"` route, partial inits) and
[#446](https://github.com/seantalts/stanli/pull/446) (`fast_math`). See
the section "0.20.0". Our reply with the rerun
([comment](https://github.com/seantalts/stanli/issues/422#issuecomment-6121453941),
2026-10-12) closed #422.

## Status

Exploratory; nothing outside this folder and `.gitignore` changed. #422
is closed: 0.20.0 fixed both things our 0.19.1 reply reported (rerun
2026-10-12, section "0.20.0"); cogmod#5 is answered and can be closed by its author. The brms suggestion is posted and under review:
[brms#1945](https://github.com/paul-buerkner/brms/issues/1945) (issue) and
[brms#1946](https://github.com/paul-buerkner/brms/pull/1946) (PR, from the
`scalar-intercept-only-dpars` branch of the fork at `../brms`, base
f131ef1). It declares an intercept-only dpar as `real` and drops its `[n]`
in the likelihood, gated off mixtures, skew-normal, logistic-normal,
`rescor`, autocorrelation, `loop = FALSE` custom families and
likelihood-block `stanvars` that mention the parameter (66 lines of R, one
updated and one new test; brms's stancode suite passes).
`brms_scalar_dpar.patch` is the diff as submitted; the PR branch is the
live version if review changes it. `brms_scalar_dpar/` has the harness
(`smoke.R`, `emit.R`, `compile_check.R`): eleven affected models compile
from both versions with the same log density and gradient.
2026-10-12: stanli 0.20.0 fixed the parallel chains with inits, and the
end-to-end rerun is done (section "0.20.0"). Next: run a few other
families through the same harness (the LNR is the only one measured here).
The reply reporting the rerun and closing #422 is posted (2026-10-12; the
draft deleted, as the others were). The performance article now has an experimental section
on the light route (the `"mock"` backend); its figures come from this
folder, and it is to be updated if brms#1946 (scalar intercept-only dpars)
or a stanli backend for brms (#1911) is merged.
2026-10-12 evening: seven multilevel decision models (section "Multilevel
decision models"): per gradient faster on five of six, per fit about
even. `issue_draft_multilevel.md` is a draft for stanli (the gap inside
the sampler, cost per trial growing with N, the DDM kernel), not posted;
delete it once posted. The DDM-7 (tag `nosmt`, 4 seeds on `long`, job
11595009) was still running.
