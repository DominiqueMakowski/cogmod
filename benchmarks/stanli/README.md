# stanli as a backend for cogmod: exact, faster per gradient, serial chains with inits

2026-10-02 to 04, on two platforms:

- **Windows**: a laptop (i7-1265U: 2 performance + 8 efficiency cores),
  CmdStan 2.38, cmdstanr 0.9.0.9000, R 4.5.3.
- **Linux**: Artemis, the Sussex HPC (AMD EPYC 9334 and 9355 nodes, Rocky
  9.5), CmdStan 2.39, cmdstanr 0.9.0, R 4.3.2.

Both run brms 2.23.1 from GitHub (`e71e9d7`; CRAN has 2.23.0), so both sample
the same Stan program. stanli 0.18.1 first, then 0.19.0 (released 2026-10-03,
in response to [seantalts/stanli#422](https://github.com/seantalts/stanli/issues/422)),
then 0.19.1 (2026-10-04) together with a hand-written LNR program the
maintainer posted as [cogmod#5](https://github.com/DominiqueMakowski/cogmod/issues/5)
(section "0.19.1, and the rewrite from cogmod#5"; on Linux these ran on EPYC
9355 and 7513 nodes).

| file | what |
| --- | --- |
| `bench.R` | the harness: `emit`, `grad` (exactness and cost per gradient), `fit` (sampling, 2 x 2) |
| `bisect.R` | cost of each way of writing the likelihood, stanli only |
| `repro_issue.R` | stanli-only reproducer of #422 |
| `repro_parallel_init.R` | stanli-only reproducer: chains run serially when given inits |
| `lnr_select_functions.stan`, `lnr_branchfree_functions.stan` | the rewritten likelihoods |
| `lnr_rewrite.stan` | cogmod#5's whole program, verbatim (variant `rw`) |
| `scalar_dpar.R` | brms's `sigma ~ 1` as a vector against a real, on a built-in family (gaussian) |
| `reply_draft.md`, `reply_cogmod5_draft.md`, `brms_issue_draft.md` | drafts for stanli#422, cogmod#5 and a brms issue; deleted once posted |
| `hpc/` | the same runs on Artemis: `run.sh`, `task.slurm`, `install.R`, `summarise.R` |
| `results/` | Windows numbers; `results/hpc/` the Linux ones (`fit_<seed>/`, `grad_<task>/`) |

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

## Short answer (stanli 0.19.1)

- **It runs the program brms writes, unchanged, and exactly.** The log density
  is CmdStan's to within 4 ULP (identical at most points) and the gradient
  to 2e-14, including in the far tail (x = -65) where a naive branch-free
  copy's gradient is not finite. At the points measured on 0.19.0 the
  gradient error is the same to every printed digit.
- **It is faster per gradient: 0.41x CmdStan on Windows, 0.71-0.77x on
  Linux** (0.19.0: 0.48x and 0.81-0.86x). CmdStan gains far more than stanli
  from Linux's faster maths library, so the gap narrows there.
- **It samples as efficiently per leapfrog step** (0.19.0, 16 seeds). Same
  posterior, same ESS per 1000 leapfrog steps (about 22.5 in every arm), so
  the difference is all in the cost of a step.
- **End to end it is slower, for now.** stanli 0.19.0 and 0.19.1 run the
  chains one at a time whenever `sample_model()` is given inits, whatever
  `parallel_chains` says (`repro_parallel_init.R`). Started from
  `cogmod_inits()`, a 4-chain fit took 2.3x CmdStan's wall time on Linux.
  Without inits the chains do run in parallel. Not yet reported upstream.
- **The maintainer's hand-written LNR (cogmod#5) is exact and faster, but
  on CmdStan the speed is not in the likelihood.** Its log density is
  brms's bit for bit. All of its CmdStan gain, and half of its stanli gain,
  comes from passing the two intercept-only sigmas as reals - one softplus
  each instead of one per trial - which brms could do for any program. The
  rewritten likelihood itself costs CmdStan 0.92-1.10x and stanli 0.83-0.94x.
- **Hand-written selects (`sel`) are 1.3x faster than the package's code in
  stanli**, but they cost CmdStan 1.36x (Windows) to 1.9x (Linux) and are
  not quite the package's numerics. stanli will not produce them
  automatically (it is held to CmdStan's results). Not worth a second Stan
  source in every family.

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

| `log Phi` | `t <= ndt` | mixture | checks | 0.18.1 Win | 0.19.0 Win | 0.19.0 Linux | 0.19.1 Linux |
| --- | --- | --- | --- | ---: | ---: | ---: | ---: |
| *`orig`, unmodified* | | | | refused | 1700 | 1680 | 1760 |
| erfc only (`bf`) | clamp to 1e-300 | `log_mix()` | none | 900 | 900 | 825 | 1070 |
| select | clamp to 1e-300 | `log_mix()` | none | 1000 | 1100 | 1050 | 1360 |
| select | mask | `log_mix()` | none | 1100 | 1100 | 1060 | 1370 |
| select | mask | `log_mix()` | data only (`dec`, `Y > 0`) | 4800 | 5000 | 5070 | 6370 |
| erfc only | clamp | `log_sum_exp()` + `log1m()` | none | 4700 | 4800 | 4075 | 5770 |
| erfc only | blend: `w * log_mix(...) + (1 - w) * (...)` | `log_mix()` | none | 5500 | 6200 | 5060 | 6470 |
| select | blend | `log_mix()` | none | 10900 | 11600 | 9580 | 14880 |

0.19.0 Linux is the median over four tasks on an EPYC 9334 node; 0.19.1
Linux over the three on an EPYC 7513 (slower: the fourth task, on a 9355,
read 0.56-0.6x these throughout). On 0.19.0 the unmodified program costs
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

`repro_issue.R` on 0.19.0 (`results/repro_issue.txt`,
`results/hpc/grad_*/repro_issue.txt`): every construct compiles, and the
never-taken parameter branch costs 0.7x (Windows) and 0.9x (Linux) the
branch-free copy, against 2.7x on 0.18.1. On 0.19.1 (Linux,
`results/hpc/grad_0.19.1_*/repro_issue.txt`) it costs 0.4x: the branchy
copy now runs faster than the branch-free one.

## A real fit: {orig, sel} x {cmdstanr, stanli}

4 chains x (500 warmup + 500 sampling), in parallel:

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
  `as_stanfit()` and `rename_pars()` take under a second.

## Parallel chains and inits (stanli 0.19.0 and 0.19.1)

`sample_model()` runs the chains one after another whenever it gets `init`,
a shared vector or a chains x parameters matrix, whatever `parallel_chains`
says. Without `init` the same model runs them in parallel.
`repro_parallel_init.R`, a toy with a `log_mix()` likelihood, N = 5000, on
Windows, 0.19.0 (0.19.1 the same, 2026-10-04: wall/sum 0.35, 1.00, 1.00):

```
init none   wall  2.49 s, chains max  2.48, sum  7.59: wall/sum 0.33
init vector wall  8.11 s, chains max  2.59, sum  8.09: wall/sum 1.00
init matrix wall  8.73 s, chains max  3.42, sum  8.70: wall/sum 1.00
```

It is not the model: a branch, a user function, `log_mix()`,
`lognormal_lpdf()`, a data ternary all ran in parallel without inits.
`sample_model()` passes the inits and `parallel_chains` to the native
sampler unchanged, so the switch is there, and nothing in R works around it.
The 0.18.1 fit started from stanli's random inits, so whether 0.18.1 had the
same behaviour is unknown.

Until it is fixed, the light route with `cogmod_inits()` costs
`chains x` the chain time. Random inits (radius 2) were fine for this model
on 0.18.1, but `cogmod_inits()` is there because brms's random inits are a
poor start for some families (`R/cogmod_inits.R` gives the reasons family by
family).

The workaround arms (`fit --arms ...,stanli_orig_proc,stanli_sel_proc`, one
process per chain) have not produced a result yet. On Windows, seed 11, the
chain at stanli seed 1101 from `cogmod_inits()` row 1 goes onto a
slow-accumulator ridge (`b_nuone_Intercept` about -2.4 +- 2.4, R-hat 1.6 in
a short smoke run) under both `orig` and `sel`, with the same leapfrog
counts in the other three chains of each, so it is deterministic, not a
hang: ~30000 leapfrog steps in the smoke run, and the full run was stopped
after 15 minutes on that one chain (2026-10-04). Use another seed, or
several, before reading anything into it.

## The cluster run

`hpc/run.sh` drives Artemis over SSH (the VPN must be up), modelled on
`benchmarks/lba_screen/`:

```
benchmarks/stanli/hpc/run.sh push          # scripts + minimal package tree + stanli tarball
benchmarks/stanli/hpc/run.sh install       # brms e71e9d7, loo, posterior, stanli + runtime
benchmarks/stanli/hpc/run.sh submit fit 16 # one seed per task, 4 CPUs
benchmarks/stanli/hpc/run.sh submit grad 4 # grad + bisect + repro per task, 2 CPUs
benchmarks/stanli/hpc/run.sh pull          # -> results/hpc/
Rscript benchmarks/stanli/hpc/summarise.R

# a rerun on a new stanli, kept apart from the old results:
STANLI_TAG=0.19.1 benchmarks/stanli/hpc/run.sh submit grad 4 --exclude=artemis-general-02
Rscript benchmarks/stanli/hpc/summarise.R --tag 0.19.1
```

Everything goes into the experiment's own directories and R library
(`/mnt/lustre/users/psych/dmm56/cogmod_stanli`): brms from GitHub at
`e71e9d7` with the loo (2.10.1) and posterior (1.7.0) it needs, and stanli
(`STANLI_VERSION`, 0.19.1 since 2026-10-04), ahead of the production library
on `R_LIBS`. The production library
is never written to. Each task works in its own copy of the tree on scratch.
Wall time on a shared node is noise, so the fit's comparisons are within a
seed, and ESS per leapfrog step needs no clock.

`artemis-general-02`, a VM whose CPU reads "AMD EPYC-Genoa Processor",
killed R with SIGILL (exit 132) in `stanli_model()` on all three seeds it
got, while the EPYC 9334 and 9355 nodes ran stanli without a fault. Seeds
10-12 were rerun with `--exclude=artemis-general-02`. stanli's Linux runtime
needs glibc 2.28 at most; the nodes have 2.34.

## What the light route needs

- `brms:::rename_pars()` is internal. After it, `summary()`, `fixef()`,
  `loo()` and `loo_compare()` worked on the result. The comment it came from
  says refit-based functions (`reloo`, moment matching) do not.
- `cogmod_inits()` goes in through `unconstrain()` and
  `sample_model(init = <matrix>)`, or `sample_cstan(init = <list of lists>)`.
  On 0.19.0 the first runs the chains serially (above); the second was not
  tried.
- Reading cmdstanr's output into a stanfit, as the `cmdstanr_sel` arm does,
  needs `brms:::read_csv_as_stanfit()`: `rstan::read_stan_csv()` cannot
  parse CmdStan 2.38's CSV header ("object 'n_kept' not found").
- brms#1911 itself may change all this. Paul Bürkner's reply there plans to
  stop converting fits to `stanfit`, so the post-processing would be per
  backend.

## For cogmod

- **Nothing in the package needs to change.** Its own Stan code runs in
  stanli 0.19.1 exactly and faster than in CmdStan per gradient, 2.4x on
  Windows and 1.3-1.4x on Linux, with the same sampler efficiency (0.19.0).
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
- **What decides it is upstream.** In stanli: chains in parallel with inits
  (a blocker for real use), and the cliffs in branch-free code, which the
  package's code does not hit. In brms: whether #1911 merges a backend.
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
full-program cliffs. Our reply with the 0.19.0 results is drafted in
`reply_draft.md`, not posted, and now out of date: it predates 0.19.1 and
cogmod#5, and its per-chain-process numbers were never obtained.

2026-10-04: 0.19.1 released with #429, announced on #422, and
[cogmod#5](https://github.com/DominiqueMakowski/cogmod/issues/5) opened with
the rewrite (section above). Not yet answered.

## Status

Exploratory; nothing outside this folder and `.gitignore` changed. Open:
answering #422 and cogmod#5 (the serial chains with `init` are still
unreported), and the brms suggestion. The latter is prototyped
(2026-10-04): `brms_scalar_dpar.patch` against brms f131ef1 declares an
intercept-only dpar as `real` and drops its `[n]` in the likelihood, gated
off mixtures, skew-normal, logistic-normal, `rescor`, autocorrelation and
`loop = FALSE` custom families (57 lines of R, one updated and one new
test; brms's stancode suite passes). `brms_scalar_dpar/` has the harness
(`smoke.R`, `emit.R`, `compile_check.R`): eleven affected models compile
from both versions with the same log density and gradient. Issue and PR
text in `brms_issue_draft.md` and `brms_pr_draft.md`; neither posted.
Next, if stanli fixes parallel chains with inits: rerun `fit` locally and on
Artemis (`hpc/run.sh`) for end-to-end wall times, then run a few other
families through the same harness.
