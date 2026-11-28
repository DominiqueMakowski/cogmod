# nutpie through stanli: ESS per second against cmdstanr

2026-11-15, on Artemis (AMD EPYC 7513 nodes, one task on a 9355; 4 whole
cores per task): stanli 0.26.0, nutpieR 1.9.2 (built with Rust 1.100),
CmdStan 2.39, cmdstanr 0.9.0, brms 2.23.1 (GitHub `e71e9d7`), R 4.3.2.
Five multilevel decision models: three race models (LNR, LBA, RDM) and
two DDMs (DDM-4, DDM-5), 18 seeds each. The LBA and RDM were then run again
over 20 seeds, with nutpie's low-rank mass matrix added.

stanli 0.26.0 (2026-11-11) adds `sample_nutpie()`, opened as
[seantalts/stanli#535](https://github.com/seantalts/stanli/issues/535):
nutpie's sampler (nuts-rs, through [nutpieR](https://github.com/andytimm/nutpieR)
>= 1.9.0) runs on stanli's runtime over its BridgeStan interface. The log
density and gradient are stanli's own, and no C++ toolchain is needed. Only
the sampler changes: nutpie's mass-matrix adaptation (it uses gradients as
well as draws, with a shorter default warmup) and, optionally, its
low-rank-plus-diagonal mass matrix.

stanli's own NUTS already takes 0.53-0.75x cmdstanr's wall time on the
multilevel race models, but 0.97-1.11x on the DDMs (`../stanli/README.md`,
section "0.21.0"). The question here: **does nutpie buy more effective
draws per second on top of stanli's gradient, and how far is each from
cmdstanr?**

## Short answer

ESS is the smallest bulk ESS over the population-level parameters; ESS per
second is over the wall time of sampling. Against cmdstanr, the seeds
without a straggler chain (section "Stragglers"); against stanli, all 18.

| model | nutpie vs cmdstanr | stanli vs cmdstanr | nutpie vs stanli |
| --- | --- | --- | --- |
| LNR | **2.76** [2.26, 3.12] (8 seeds) | 1.61 [1.18, 2.41] | **1.64** [1.30, 2.07] |
| RDM | **3.35** [1.92, 5.80] (13 seeds) | 1.64 [0.98, 2.52] | **1.95** [1.01, 3.42] |
| LBA | **1.60** [1.21, 2.99] (10 seeds) | 1.58 [0.98, 2.18] | **1.23** [0.83, 1.72] |
| DDM-4 | **1.74** [0.86, 2.86] (17 seeds) | 1.12 [0.70, 1.51] | **1.62** [1.15, 2.61] |
| DDM-5 | **1.98** [1.70, 3.15] (10 seeds) | 0.98 [0.66, 1.42] | **2.06** [1.07, 2.87] |

- **nutpie is the fastest of the three on every model.** It takes
  0.22-0.31x cmdstanr's wall time and 0.28-0.47x stanli's. In ESS per
  second against cmdstanr:
  - LNR and RDM: 2.8-3.4x;
  - DDMs: 1.7-2.0x;
  - LBA: 1.60x, no better than stanli's own NUTS (1.58x).
- **On the DDMs, the whole gain is the sampler.** stanli's NUTS alone is
  level with cmdstanr there (1.12x and 0.98x), so nutpie's 1.6-2.1x over
  stanli is all it brings. That is about what it adds over stanli on the
  race models (1.2-2.0x).
- **The gain is fewer gradient evaluations.**
  - nutpie's warmup takes about a third of Stan's leapfrog steps on the
    race models, and a fifth to a quarter on the DDMs.
  - Its sampling takes 0.4-0.65x Stan's steps, at step sizes 2.1-4.5x larger.
  - A step costs what it costs in stanli's sampler (wall time per step of
    the longest chain): 0.44 against 0.41 ms on the LNR, 1.28 against 1.27
    on the RDM, 0.99 against 0.96 on the LBA, 1.61 against 1.54 on the
    DDM-4, and 1.82 against 1.84 on the DDM-5.
- **But nutpie returns fewer effective draws** for the same 4 x 500 draws.
  - Against stanli: 0.82x on the RDM, 0.78x on the LNR, 0.70x on the
    DDM-5, 0.47x on the DDM-4, and **0.41x on the LBA**.
  - Its weak spot is the population-level intercepts. On the LBA every
    dpar's intercept gets a quarter to a half of Stan's bulk ESS, and the
    SDs of the participants' boundaries and non-decision times get about
    60%. The Speed effects, the correlation and the other SDs stay within
    about 20%.
  - On the RDM it is the intercepts of `boundary`, `ndt` and `sigmabias`.
  - On the DDM-4 it is above all `b_ndt_Intercept`, at 0.44x stanli's.
  - Sections "The LBA" and "The DDMs" have the details.
- **The posterior is the same, but nutpie's ESS looks a little optimistic
  on the LBA and RDM.**
  - No parameter's mean is off consistently across seeds, and posterior
    SDs agree within 4%.
  - The standardised difference from cmdstanr is the mean difference over
    the joint MCSE. Its spread should have SD 1.
  - For nutpie that SD is 1.09 on the LBA, 1.06 on the RDM, 1.04 on the
    DDM-5, and 0.99-1.00 on the LNR and DDM-4. For stanli it is 0.95-1.02.
  - If all of the excess is nutpie's, its MCSEs are too small by up to about
    15-20% in variance.
- **nutpie diverges more**, out of 36,000 draws per arm:

  | model | nutpie | cmdstanr | stanli |
  | --- | --- | --- | --- |
  | LBA | 25 in 13 of 18 seeds | 6 in 5 | 7 in 5 |
  | DDM-5 | 6 in 4 | 1 | 1 |
  | LNR | 5 in 3 | 3 | 1 |
  | RDM | 4 in 3 | 0 | 0 |
  | DDM-4 | 1 | 1 | 0 |

  R-hat is at most 1.04 everywhere.
- **With the compile counted**, nutpie gives more ESS per second of build
  plus wall time than cmdstanr:
  - LNR: 4.7x, taking CmdStan's compile as 43 s uncontended;
  - RDM: 4.5x;
  - DDM-5: 2.5x;
  - DDM-4: 2.1x;
  - LBA: 2.0x.

  On the DDMs and the LBA, sampling takes 2-5 minutes, so the compile
  matters less.

- **The low-rank mass matrix (`adaptation = "low_rank"`) is worse**
  (section "The low-rank mass matrix"). On the LBA and RDM it recovers
  part of the intercepts' deficit, but it loses more on the participants'
  SDs. That leaves it at 0.43x (LBA) and 0.34x (RDM) diagonal nutpie's ESS
  per second, better in 1 seed of 20 on each model. It also diverges in
  every LBA seed: 315 divergences, against 28 for the diagonal. This was
  run with 500 warmup draws; nuts-rs's own default for low rank is 800.

## Results

`summarise.R` on the 90 tasks (`results/summary.csv`). Medians [min, max]
over seeds; ratios are against cmdstanr in the same seed. Steps are
leapfrog steps summed over the 4 chains, in thousands.

| model | arm | wall (s) | min bulk ESS | ESS/s | ESS/s vs cmdstanr, 18 seeds | without stragglers | warmup steps | sampling steps | steps per ESS | step size |
| --- | --- | --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LNR | cmdstanr | 62.6 [41.9, 172.3] | 549 [443, 671] | 8.6 | 1 | 1 | 194 | 95 | 517 | 0.112 |
| LNR | stanli | 32.8 [28.8, 42.3] | 599 [476, 679] | 17.6 | 2.30 [1.18, 5.62] | 1.61 [1.18, 2.41] | 193 | 93 | 484 | 0.111 |
| LNR | nutpie | 15.1 [13.5, 23.8] | 463 [408, 561] | 29.3 | 3.14 [2.25, 7.42] | 2.76 [2.26, 3.12] | 71 | 58 | 285 | 0.271 |
| RDM | cmdstanr | 274 [216, 360] | 511 [436, 711] | 1.8 | 1 | 1 | 288 | 129 | 824 | 0.072 |
| RDM | stanli | 139 [125, 219] | 521 [311, 643] | 3.8 | 1.65 [0.98, 2.86] | 1.64 [0.98, 2.52] | 288 | 129 | 813 | 0.071 |
| RDM | nutpie | 58.2 [50.5, 96.0] | 405 [302, 516] | 6.9 | 3.90 [1.92, 5.80] | 3.35 [1.92, 5.80] | 96 | 76 | 418 | 0.218 |
| LBA | cmdstanr | 259 [203, 372] | 494 [388, 609] | 1.8 | 1 | 1 | 421 | 259 | 1378 | 0.038 |
| LBA | stanli | 204 [150, 317] | 568 [493, 695] | 2.9 | 1.58 [0.98, 2.18] | 1.58 [0.98, 2.18] | 420 | 254 | 1183 | 0.040 |
| LBA | nutpie | 73.3 [52.4, 84.7] | 231 [185, 326] | 3.2 | 1.77 [1.21, 2.99] | 1.60 [1.21, 2.99] | 130 | 100 | 994 | 0.179 |
| DDM-4 | cmdstanr | 145 [130, 204] | 674 [492, 763] | 4.1 | 1 | 1 | 238 | 63 | 450 | 0.138 |
| DDM-4 | stanli | 129 [109, 171] | 629 [462, 730] | 4.6 | 1.12 [0.70, 1.51] | 1.12 [0.70, 1.51] | 238 | 63 | 484 | 0.137 |
| DDM-4 | nutpie | 40.4 [31.2, 52.2] | 296 [188, 341] | 8.0 | 1.83 [0.86, 2.86] | 1.74 [0.86, 2.86] | 51 | 34 | 287 | 0.294 |
| DDM-5 | cmdstanr | 220 [164, 387] | 611 [427, 711] | 2.7 | 1 | 1 | 247 | 74 | 515 | 0.123 |
| DDM-5 | stanli | 160 [142, 247] | 553 [444, 752] | 3.4 | 1.13 [0.66, 2.84] | 0.98 [0.66, 1.42] | 250 | 74 | 592 | 0.123 |
| DDM-5 | nutpie | 51.4 [44.9, 73.2] | 375 [257, 425] | 6.6 | 2.25 [1.69, 5.61] | 1.98 [1.70, 3.15] | 60 | 46 | 281 | 0.278 |

The LNR's CmdStan compile took 96 s on the node that ran 13 of its tasks
(all compiling at once), 43-61 s on the other two, and 41 s for the smoke
task alone, so the 7.2x that `summarise.R` reports for its ESS per second of
build + wall is inflated by the crowded compiles; hence the 43 s figure
above. The LBA, RDM and DDM arrays were spread more thinly (compile 30-83 s).

### The LBA

On the LBA, nutpie adapts to a step size 4.5x Stan's (0.18 against 0.04)
and takes 3.2x fewer warmup steps and 2.5x fewer sampling steps per draw,
but its draws are much more autocorrelated on the population intercepts.
Median bulk ESS over seeds, nutpie against cmdstanr: `b_ndt_Intercept` 265
vs 822, `b_boundary_Intercept` 284 vs 995, `b_Intercept` 536 vs 986,
`b_driftone_Intercept` 644 vs 1080, `b_sigmabias_Intercept` 668 vs 1492,
`b_sigmaone_Intercept` 895 vs 1652, and the participants' SDs of
`boundary` 320 vs 497 and of `ndt` 508 vs 849. The Speed effects (1901 vs
1676, 2351 vs 2338), the correlation (1518 vs 1683) and the other SDs
(678-751 vs 704-932) stay close. The `ndt` and `boundary` intercepts set nutpie's minimum in 12
of its 18 seeds, the SD of the participants' boundaries in the other 6 (for
Stan, that SD in 15 or more). The RDM shows the same pattern on fewer
parameters (`b_boundary_Intercept` 512 vs 1362, `b_ndt_Intercept` 498 vs
837, `b_sigmabias_Intercept` 629 vs 1970; the other intercepts within 20%),
and an intercept sets nutpie's minimum in 8 of 18 seeds, a random-effect SD
in the other 10, so it keeps most of its lead.

stanli also skips its O1 optimisation on the LBA program (structural cost
26,516 against a budget of 20,000) and runs the unoptimised one, for the
stanli and nutpie arms alike; `STANLI_NO_O1_FALLBACK=1` forces O1. That is
what a user gets by default, so it was left on; it may be part of why
stanli's lead per step is smallest there (0.96 against cmdstanr's 1.28 ms).

### The DDMs

The DDM-4 fixes the three variabilities at 0 and runs on Stan's 4-argument
`wiener_lpdf()`. The DDM-5 estimates `sigmadrift` and runs on the
5-argument one.

**stanli against cmdstanr.** On these programs stanli's step is a little
cheaper than CmdStan's. In the sampling phase it costs 0.91x on both DDMs
(1.54 against 1.69 ms, 1.84 against 2.04). Its adaptation takes the same
steps as CmdStan's, so its ESS per second ends up near cmdstanr's: 1.12x
and 0.98x. On the DDM-4 this is a change from `../stanli/README.md`'s 0.21.0
run, where a step cost 1.15x CmdStan's. Since then both the program (the
DDM-4 router change, `benchmarks/ddm_router/`) and stanli have changed, so
this run does not say which of the two did it.

**nutpie's weak spot** is the same as on the LBA, but milder. Median bulk
ESS, nutpie against stanli, on the DDM-4:

- `b_ndt_Intercept`: 296 against 666. It sets nutpie's minimum in 17 of 18
  seeds. stanli's minimum is `b_ndt_Intercept` in 7 seeds, the
  participants' drift SD in 7 and the boundary SD in 4.
- The participants' `ndt` SD: 0.49x.
- `b_boundary_Intercept`: 0.52x.
- The participants' `boundary` SD: 0.58x.
- The bias, drift and correlation terms: 0.60-0.67x.
- The Speed effects: 0.82-0.88x.

The DDM-5 is more even. The intercepts of `boundary`, `bias`, `ndt` and
`sigmadrift`, and the participants' `ndt` and `boundary` SDs, get 0.62-0.67x
stanli's ESS, and the Speed effects 0.90-0.98x. nutpie still comes out
1.6x (DDM-4) and 2.1x (DDM-5) ahead of stanli in ESS per second, because
it takes far fewer steps.

### The low-rank mass matrix

The LBA and RDM were rerun with a fourth arm, `sample_nutpie(adaptation =
"low_rank")`, at the same 500 warmup draws. The run is tagged `lowrank` and
summarised in `results/summary_lowrank.csv`: 20 seeds per model, 38 tasks
on EPYC 7513 nodes and 2 on a 9355. Each arm's draws are fixed by its seed,
data and inits. So for seeds 1-18, the other three arms reproduced the first
run's draws exactly, and only their timings are new. Their ratios agree with
the first run: nutpie against stanli gives 1.27 on the LBA and 2.00 on the
RDM, against 1.23 and 1.95.

The table gives medians over seeds. ESS per second against cmdstanr leaves
out the seeds with a straggler: 1 of 20 on the LBA, 5 on the RDM.

| model | arm | min bulk ESS | ESS/s vs cmdstanr | steps per ESS | sampling steps (k) | step size | divergences (seeds) | max R-hat |
| --- | --- | --- | --- | --- | --- | --- | --- | --- |
| LBA | stanli | 565 | 1.56 | 1192 | 254 | 0.040 | 7 (5) | 1.015 |
| LBA | nutpie | 231 | 1.83 | 994 | 100 | 0.177 | 28 (15) | 1.042 |
| LBA | nutpie_lowrank | 105 | 0.88 | 1903 | 59 | 0.234 | 315 (20) | 1.103 |
| RDM | stanli | 521 | 1.56 | 813 | 129 | 0.071 | 0 | 1.017 |
| RDM | nutpie | 410 | 3.21 | 413 | 76 | 0.218 | 5 (4) | 1.021 |
| RDM | nutpie_lowrank | 140 | 1.26 | 1113 | 41 | 0.261 | 25 (13) | 1.056 |

Within each seed, low rank gives 0.43x [0.19, 1.14] the diagonal's ESS
per second on the LBA, and 0.34x [0.21, 1.13] on the RDM. It is also below
stanli's own NUTS on both models (0.53x and 0.83x).

**Population intercepts.** On these the low-rank matrix does what it was
added for. Median bulk ESS, low rank against diagonal (LBA, then RDM):

- `b_Intercept` (LBA): 1024 against 536. cmdstanr has 986.
- `b_ndt_Intercept`: 324 against 265, and 716 against 498.
- `b_boundary_Intercept`: 338 against 284, and 642 against 512.

**Participants' SDs.** On these it loses more than it gained:

- `sd_Participant__boundary_Intercept`: 105 against 332 (LBA), and 147
  against 562 (RDM). This SD sets low rank's minimum ESS in 18 of 20 seeds
  on each model.
- `sd_Participant__ndt_Intercept`: 249 against 496, and 251 against 599.
- The other SDs: 0.5-0.65x the diagonal's.
- The participants' boundary deviations themselves: 0.35-0.6x.

**Steps.** Low rank takes 0.54-0.59x the diagonal's sampling steps, at a
step size 20-30% larger, and its maximum tree depth is one lower. Each draw
moves less, though: low rank needs 2.2x (LBA) and 2.7x (RDM) the leapfrog
steps per effective draw. Its warmup takes 1.08-1.16x the diagonal's steps.

**Divergences.** On the LBA, low rank diverged 315 times: in every seed,
and in 63 of 80 chains (one chain had 57). The diagonal diverged 28 times
in 15 seeds. On the RDM, low rank diverged 25 times in 13 seeds, against
the diagonal's 5 in 4. On the LBA, seeds with more divergences have a
lower boundary-SD ESS (Spearman -0.70). R-hat passes 1.05 for 29
parameter-seed pairs in 6 seeds, all of them participant SDs or
deviations (boundary, ndt, the drifts).

**The posterior** still agrees with cmdstanr's. No parameter is off
consistently across seeds, sd(z) is 1.09 (LBA) and 1.08 (RDM), and
posterior SDs agree within 5%. But a fit with that many divergences is not
one to keep.

A likely reason, not measured, since the draws are not saved: in brms's
non-centred parameterisation, each participant SD and its standardised
deviations lie on a curved ridge, and no fixed linear metric fits a curved
ridge. A metric that straightens the population directions then adapts a
step size too large for the curvature around the SDs.

Not tested: nuts-rs's own warmup for low rank, 800 draws, where this run
used 500. Low rank estimates a covariance from the warmup draws, so it
may need more of them than the diagonal does.

### Stragglers

Many tasks shared a node (13 of the LNR's on one), and in some of them a
chain spent 1.5-5.3x its usual time per leapfrog step for a whole phase,
with the same number of steps as the other chains: the node, not the
sampler. Chains over 1.5x the median of their model and arm:

| model | cmdstanr | stanli |
| --- | --- | --- |
| LNR | 12 chains in 10 seeds, up to 5.3x | 6 in 3, up to 1.6x |
| LBA | 7 in 5, up to 2.4x | 10 in 5, up to 2.6x |
| RDM | 3 in 2, up to 1.7x | 8 in 3, up to 1.7x |
| DDM-4 | 1 in 1, up to 1.6x | none (at most 1.48x) |
| DDM-5 | 6 in 6, up to 3.4x | 5 in 2, up to 1.6x |

`summarise.R` drops any seed with such a chain in either arm from the
ratios against cmdstanr. That is 10 of 18 seeds on the LNR, 8 on the LBA,
5 on the RDM, 1 on the DDM-4 and 8 on the DDM-5. On the LNR the worst cases
took cmdstanr's wall time from about 50 s (42-63 s in the clean seeds) to
100-170 s, which is why its 18-seed ratios run higher and wider.

nutpie reports no phase times, so its own stragglers cannot be seen. The
ratios of nutpie against stanli therefore use all 18 seeds. Dropping the
seeds where stanli straggled leaves their medians at:

- LNR: 1.64;
- LBA: 1.23;
- RDM: 1.86, from 1.95;
- DDM-4: 1.62, unchanged (stanli had no stragglers);
- DDM-5: 1.99, from 2.06.

## Design

`bench.R`, one (model, seed) per call:

- **Data**: 10 participants x 100 trials (50 Accuracy, 50 Speed), simulated
  from each family's RNG without outliers, with each participant's own
  value of the parameter that carries the choice, of the second accumulator
  or the boundary, and of the non-decision time. Models and truths are
  `../stanli/multilevel.R`'s `lnr`, `lba`, `rdm`, `ddm4` and `ddm5`, so the two benchmarks
  fit the same things. The data seed is fixed, so seeds differ only in the
  samplers' randomness.
- **Models** (with `cogmod_priors()`, `cogmod_stanvars()` and the same
  `cogmod_inits()` in every arm):
  - LNR: `nuzero ~ Condition + (Condition | Participant)`,
    `nuone ~ Condition + (1 | Participant)`, `ndt ~ 1 + (1 | Participant)`,
    `sigmazero ~ 1`, `sigmaone ~ 1`, `sigmabias = 0`;
  - LBA (`cogmod_lba2()`): `driftzero ~ Condition + (Condition | Participant)`,
    `driftone ~ Condition + (1 | Participant)`, `boundary ~ 1 + (1 | Participant)`,
    `ndt ~ 1 + (1 | Participant)`, `sigmaone ~ 1`, `sigmabias ~ 1`,
    `sigmazero = 1`;
  - RDM: `vzero ~ Condition + (Condition | Participant)`,
    `driftone ~ Condition + (1 | Participant)`, `boundary ~ 1 + (1 | Participant)`,
    `ndt ~ 1 + (1 | Participant)`, `sigmabias ~ 1`;
  - DDM-4 (`cogmod_ddm()`): `drift ~ Condition + (Condition | Participant)`,
    `boundary ~ Condition + (1 | Participant)`, `bias ~ 1`,
    `ndt ~ 1 + (1 | Participant)`, `sigmadrift = sigmabias = sigmandt = 0`
    (Stan's 4-argument `wiener_lpdf()`);
  - DDM-5: the same with `sigmadrift ~ 1` (its 5-argument `wiener_lpdf()`).
- **Arms**, in an order that rotates with the seed:
  - `cmdstanr`: brms's program compiled by CmdStan;
  - `stanli`: `sample_model()`;
  - `nutpie`: `sample_nutpie()`, diagonal mass matrix;
  - `nutpie_lowrank`: `adaptation = "low_rank"`, only when `--arms` asks for it.

  Each arm runs 4 chains in parallel, 500 warmup + 500 sampling draws.
  nutpie gets the same warmup unless `--nutpie_warmup default` (nuts-rs's
  own: 400 diagonal, 800 low rank).
- **Clock**: from the call that starts sampling to the `brmsfit` (every arm
  goes through `brms:::rename_pars()`, as in `../stanli/multilevel.R`).
  The compile (cmdstanr) or build (stanli, about 1 s) is timed apart, as
  `build_s`.
- **ESS**: the smallest bulk ESS over the population-level parameters
  (tail and median also kept). ESS per second is that over the wall time
  (`ess_bulk_per_wall_s`) or over build + wall (`ess_bulk_per_total_s`);
  `*_vs_cmdstanr` divides each arm's by cmdstanr's in the same seed.
- **Without a clock**: warmup draws are kept, so leapfrog steps are counted
  in warmup and sampling. stanli and nutpie evaluate the same gradient, so
  `leapfrog_per_ess` (all leapfrog steps over the min bulk ESS) compares
  the two adaptations directly.
- **Checks**: R-hat, divergences, and each posterior mean against
  cmdstanr's in units of their joint MCSE: per pair (`mean_z_max`), and in
  `summarise.R` over every parameter and seed (count beyond 3, SD of the
  standardised differences, any parameter off consistently across seeds,
  ratio of posterior SDs).

nutpie reports one time for warmup and sampling together, so its
`warmup_s` is NA and its `ms_per_leapfrog` (sampling phase) cannot be
computed. Its wall time is measured like the others'.

## Files

| file | what |
| --- | --- |
| `bench.R` | simulate, fit every arm, write `results/{fit,chains,pars}_<model>_<seed>[_<tag>].csv` |
| `summarise.R` | one row per (model, arm) over seeds: median [min, max], ratios paired within seed |
| `hpc/run.sh` | `push`, `install`, `smoke`, `submit`, `queue`, `progress`, `log`, `pull` on Artemis |
| `hpc/task.slurm` | one seed per array task, in its own scratch copy of the tree |
| `hpc/install.sh`, `hpc/install.R` | Rust (rustup, into the experiment's directory), stanli 0.26.0 + runtime, nutpieR from source; checks each engine on a toy model |

## Running it

Locally (the nutpie arms need stanli >= 0.26.0 and nutpieR >= 1.9.0; on
Windows nutpieR has an r-universe binary for R 4.5):

```bash
Rscript benchmarks/nutpieR/bench.R --seed 1 --model lnr
```

On Artemis (VPN up). nutpieR has no binary for R 4.3 on Linux, so `install`
builds it from source with a rustup toolchain put in the experiment's
directory, about 0.5 GB, once. `NP_MODELS` (colon-separated, default `lnr`)
says which models an array covers; `submit` defaults to 18 seeds of each,
and `smoke` runs one short task per model:

```bash
benchmarks/nutpieR/hpc/run.sh push
benchmarks/nutpieR/hpc/run.sh install
NP_MODELS=lba:rdm benchmarks/nutpieR/hpc/run.sh smoke
NP_MODELS=lba:rdm benchmarks/nutpieR/hpc/run.sh submit
NP_MODELS=ddm4:ddm5 benchmarks/nutpieR/hpc/run.sh submit
benchmarks/nutpieR/hpc/run.sh pull
Rscript benchmarks/nutpieR/summarise.R
```

The LNR's 18 seeds ran before the model mapping, as `results/hpc/fit_<seed>/`;
the later runs are in `fit_<model>_<seed>/`. `summarise.R` reads both.

The low-rank arm ran as a tagged rerun. With 4 arms rotating over 4 seeds,
that is 20 seeds per model rather than 18, so 40 tasks for 2 models:

```bash
NP_TAG=lowrank NP_MODELS=lba:rdm NP_ARGS="--arms cmdstanr:stanli:nutpie:nutpie_lowrank" benchmarks/nutpieR/hpc/run.sh submit 40
Rscript benchmarks/nutpieR/summarise.R --tag lowrank --out benchmarks/nutpieR/results/summary_lowrank.csv
```

Low rank is not pursued further. Not run yet: diagonal nutpie at nuts-rs's
own warmup of 400 draws, rather than Stan's 500:

```bash
NP_TAG=npwarm NP_ARGS="--arms stanli:nutpie --nutpie_warmup default" benchmarks/nutpieR/hpc/run.sh submit
```

The library is the experiment's own (`cogmod_nutpie/R_libs`: brms 2.23.1
from GitHub as on the laptop, stanli 0.26.0, nutpieR), ahead of production
on `R_LIBS`; nothing else on the cluster is touched.
