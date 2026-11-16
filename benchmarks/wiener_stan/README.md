# A Stan-language Wiener density for the DDM-4 and DDM-5

Prototyped and applied 2026-11-16 (0.3.6). The candidate below,
`wfpt_unrolled.stan`, is now `.WFPT_STAN` in `R/model_ddm.R`. The R side's
`.ddm_lfpt0()` took the same two refinements, and the series switch is one
constant, `.DDM_FPT_SWITCH`, for both. The prototype measurements below ran
against scratch worktrees of 0.3.5 (6463109) patched by `apply_variant.R`.
That script and `model_ddm_unrolled.patch` apply only to that tree.

The question was whether cogmod can use the speed-ups proposed upstream for
Stan's `wiener_lpdf()` without waiting for a Stan Math release, and on
stanli as well as CmdStan:
[stan-dev/math#3481](https://github.com/stan-dev/math/issues/3481) (speed, 4-
and 5-parameter forms),
[stan-dev/math#3442](https://github.com/stan-dev/math/pull/3442) and
[stan-dev/math#3435](https://github.com/stan-dev/math/pull/3435) (accuracy and
underflow, mostly the 7-parameter form).

stanli cannot load user C++: models are lowered to ops over precompiled
Stan Math kernels, and `wiener_lpdf()` is one of them, run on a nested tape.
So the density is written in the Stan language, in the prelude, where it
runs on both engines (and rstan). It replaces the DDM-4 and DDM-5 call
site of `cogmod_ddm_decision_lpdf()`: the certificate, the router and both
Stan calls there. The 7-parameter call site is untouched.

## Short answer

`wfpt_unrolled.stan`, per gradient, against cogmod today (0.3.5 @ 6463109).
The setup is `program_ab.R` on the multilevel programs (N = 1000), and
`gradient_cost.R` on the intercept-only ones (N = 5000):

| | CmdStan | `gradient_cost.R` | stanli | stanli vs CmdStan today |
| --- | ---: | ---: | ---: | ---: |
| DDM-4 | **0.86** | 0.80 | 1.06 | 0.76 |
| DDM-5 (`sigmadrift` estimated) | **0.29** | 0.34 | 0.67 | **0.24** |

- **Accuracy.** Within 9.4e-12 (relative) of an independent reference at
  70,000 points, which is the reference's own spread. Stan's wiener5, which
  every DDM-5 trial uses today, is off by up to 1.7e-2 relative near a log
  density of -28 (0.47 log units). In the tails it is off by up to 0.11,
  and its gradients by up to 84x.
- **Finite everywhere.** No point gives a non-finite value or partial,
  including boundaries from 1e-160 to 1e80, where the current code gives a
  finite value with a non-finite gradient at 627-2252 of 5000 points.
- **DDM-5 at DDM-4 cost.** A DDM-5 gradient now costs about what a DDM-4
  gradient does today (1964 against 1786 us), and 1.6x the new DDM-4. The
  current code gives 3.7x.

## Files

| file | what |
| --- | --- |
| `wfpt_direct.stan` | the R side's `.ddm_lfpt0()` / `.ddm_lfpt()` transcribed: an `exp()` per term |
| `wfpt_rec.stan` | the same series with the per-term `exp()`/`sin()` turned into products |
| `wfpt_unrolled.stan` | `rec` written out at fixed length with no loop: **the candidate** (stanli, below) |
| `apply_variant.R` | patches a worktree's `R/model_ddm.R` with one variant, and adds a `cogmod_ddm5` program to its `benchmarks/gradient_programs.R` |
| `model_ddm_unrolled.patch` | the candidate as a diff against 6463109 (`git apply`) |
| `grid.R` | function level: values and partials of every variant and of Stan's own forms, against the reference |
| `program_ab.R` | program level: agreement with the base at posterior draws and tail points, cost per gradient in CmdStan and stanli |
| `tail_audit.R` | at the program point where a variant and the base differ most, scores every trial against the reference |
| `ladder.R` | per-gradient cost of each between-trial variability relative to the `sigmadrift`-only model, the ratios `?cogmod_ddm` quotes |
| `results/` | the runs quoted below (`*_run1` = `direct` and `rec`, before the guard); the CSVs are gitignored like every benchmark CSV |

To reproduce, from the repository root:

```
git worktree add --detach <tmp>/base HEAD
git worktree add --detach <tmp>/proto HEAD
Rscript benchmarks/wiener_stan/apply_variant.R --tree <tmp>/base --variant none
Rscript benchmarks/wiener_stan/apply_variant.R --tree <tmp>/proto --variant unrolled
Rscript benchmarks/wiener_stan/grid.R --base <tmp>/base
Rscript benchmarks/wiener_stan/program_ab.R --base <tmp>/base --protos unrolled=<tmp>/proto --model ddm5
(cd <tmp>/proto && Rscript benchmarks/gradient_check.R --pkg . --families cogmod_ddm4,cogmod_ddm5)
```

`program_ab.R` reads `benchmarks/results/stanli_ml/<model>_inputs.rds`
(`benchmarks/stanli/multilevel.R emit`). Remove the worktrees afterwards.

## The density

The density is Navarro & Fuss (2009), as on the R side: the driftless,
unit-boundary density at rescaled time u = t / a^2, with 1e-12 absolute
accuracy on it.
- The small-time series runs over k = -m..m.
- The large-time series runs over k = 1..K.
- Each is summed with its leading term's exponent factored out, so the log
  stays finite and relatively accurate in both tails.
- The drift and its spread enter in closed form, so `sv` adds a factor, not
  a series.

It keeps `wiener_lpdf()`'s convention: the upper boundary, with `w` measured
from the lower one. The caller is unchanged.

How it maps onto the upstream proposals:

- **#3481 item 1 (term counts as doubles).** The switch point and the term
  counts are thresholds on u, computed once in R (`u* = 0.2028`, where the
  two counts cross). Comparisons record nothing on the tape.
- **#3481 item 2 (the large-time bound).** The bound is the paper's, as the
  R side has it. Stan's 4-parameter form, with its constants, keeps the
  small-time series up to u = 0.96 at 7 terms.
- **#3481 items 4, 5 and 7 (one series for all partials).** Autodiff
  through one series gives every partial: the density, d/dt, d/da and d/dw
  are not four separate series as in wiener5. This is where the DDM-5's
  0.29 comes from.
- **#3481 item 6 (summation in linear space).** The sums are already
  linear, after the leading exponent is factored out.
- **#3442 (tolerances relative to the value).** Factoring the leading term
  turns the absolute bound into relative accuracy in both tails, as
  measured below.
- **#3435, and #3442's 7-parameter work.** Not covered: integrating over
  `sigmabias` or `sigmandt` needs quadrature, which the Stan language only
  offers at fixed nodes.

Two things go beyond the R side, so R and Stan now differ below 2e-11:

- The small-time k = 0 and k = -1 terms are summed as
  `-expm1(-2w/u) - w (1 + e1)`. Written as two terms, they cancel as the
  start point nears the far boundary.
- The large-time `sin(k pi wl)` is taken at whichever of `wl` and `w` is
  nearer 0, with alternating signs. Each keeps its relative precision near
  its own boundary.

The wrapper returns -inf as a constant outside 1e-140 < u < 1e140, where
the log density is below -1e100 or so. This came out of `program_ab.R`'s
first run. At a DDM-5 tail point a random effect had put `boundary` at
8e-96, so u = t / a^2 reached 1e190. Every trial's value was finite, but
the division's own partial, u / a^2 = 4e380, overflowed on the tape. 0 * inf
is NaN under any mixture weight, so the whole gradient was NaN. The
comparison is written on products (`t > 1e140 * a2`), whose partials cannot
overflow. Stan's wiener5 gives a finite gradient at that point, but a wrong
one: d/dlog t = -11 where the value is -1.4e191, and the truth is about
the value itself.

## Why `unrolled`: stanli

The first run timed `direct` and `rec`. On CmdStan they were 0.93 / 0.84 on
the DDM-4 and 0.35 / 0.31 on the DDM-5. On stanli they were 23x and 15x
*slower* than the base.

`STANLI_EXECUTION_REPORT=1` shows why. In the base, the trial loop is one
`OP_REGION_MAP` (613 instructions, with `wiener_lpdf()` as an
`OP_WIENER_LPDF` kernel call). With `rec` it stays an `OP_LOOP` with
`control_selection: runtime`, and it is recorded afresh at every new point:
1.25 ms per gradient at a repeated point, 23.5 ms cycling over 40 draws.
That is the same mechanism as the 13x of the Wald quadrature
(`benchmarks/stanli/README.md`, "0.23.2"). The inner loops' trip counts
depend on parameters (through the term-count thresholds), and that keeps
the region map out.

`unrolled` writes both series at their longest: m = 3 and K = 6, the most
either needs at 1e-12 on its side of u*. The extra terms are below 1e-12
of the sum. With no inner loop, stanli maps the trial loop again (1.25 ms
cycling, base 1.0). On CmdStan it costs what `rec` does (0.859 against
0.856). Its values agree with `rec` to 1e-14 at the program's draws.

**Any future edit that puts a loop back into this function brings the
15-24x cliff back on stanli.** The call-site comment would have to say so.

## Results

Windows 11 (i7-1265U), CmdStan 2.38.0, stanli 0.26.0. stanli was upgraded
from 0.23.2 that day: its 0.26.0 runtime was already installed, from the
nutpie work.

### Function level (`grid.R`)

5000 points per batch: 25,000 log-uniform over wide ranges (as
`benchmarks/ddm_router/grid_ab.R`), 10,000 like a DDM-4 posterior's
trials, and 5000 at extreme boundaries, in each of two modes (sv = 0, and
sv > 0).

The reference is in R and independent of the package:
- the small-time series to k = +-30 below u = 0.05;
- the large-time series to K = 200 above it.

Where both converge (u in [0.05, 0.4]) they agree to 9.45e-12. The
package's `.ddm_lfpt()` is within 2.0e-11 of it.

Value error is |value - ref| / max(1, |ref|), with the body at ref > -30, the
tail at -30 to -1000, and far below that. Gradient error is against
central differences of the reference. A partial whose differences at h and
2h disagree by more than 1e-7 is not judged.

| sv | density | body value | body grad | tail value | tail grad | far value | far grad | non-finite |
| --- | --- | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 0 | cogmod today | 6.6e-4 | 2.0e-2 | 1.4e-3 | 32 | 6.3e-4 | 190 | 0 |
| 0 | classic `wiener_lpdf()`, raw | 6.6e-4 | 2.0e-2 | 4.5e-4 | 3.3e-3 | 9.7e-4 | 6.7e-4 | 1727 values, 43 partials |
| 0 | wiener5 at sv = 0, raw | 1.7e-2 | 0.28 | 3.2e-2 | 75 | 6.3e-4 | 190 | 0 |
| 0 | **unrolled** | **9.4e-12** | 4.1e-7 | **1.6e-13** | 3.6e-7 | **8.2e-15** | 6.7e-4 | 0 |
| > 0 | cogmod today (= wiener5) | 1.7e-2 | 0.38 | 0.11 | 84 | 1.0e-3 | 71 | 0 |
| > 0 | **unrolled** | **9.4e-12** | 4.1e-7 | **1.3e-13** | 3.6e-7 | **3.2e-15** | 7.7e-6 | 0 |

`rec` matches `unrolled` in every cell. `direct` is within 2e-11 in value
and has the same gradients. The 4e-7 gradient figure is the finite
differences' floor: every sound variant reaches it, and so does the
classic form wherever that one is sound. On the 10,000 posterior-like
points, every density agrees with the reference to 3.2e-9 or better, so
these differences live off the posterior.

At boundaries from 1e-160 to 1e-3 and from 200 to 1e80 (5000 points),
these are the counts of a finite value with a non-finite partial:

| density | sv = 0 | sv > 0 |
| --- | ---: | ---: |
| cogmod today | 2252 | 627 |
| wiener5 | 625 | 627 |
| classic | 1231 | - |
| `direct` | 85 | 85 |
| `rec`, `unrolled` | **0** | **0** |

On top of that, the classic form and the current code return -inf with
non-finite partials at 2607 and 67 points. `unrolled` returns a constant
-inf at 1800 points (the guard) and never a -inf with partials.

### Gradient check

`gradient_check.R --families cogmod_ddm4,cogmod_ddm5` passes in the base
and in every variant: 16 of 16 and 18 of 18 points, worst relative error
against central differences 5.9e-8 and 1.7e-7. The figures are the same
in all four trees. `cogmod_ddm5` (`sigmadrift` estimated, the other two
variabilities 0) is a program `apply_variant.R` adds. No CI program reaches
the 5-parameter path today.

### Program level (`program_ab.R`)

These are `multilevel.R`'s ddm4 and ddm5 programs: 10 participants x 2
conditions x 50 trials. The points are 40 posterior draws of a short fit of
the base, and 80 points moved from them by N(0, 1.5) and N(0, 3) on the
unconstrained scale. Each figure is the largest gap to the base:

| program | at draws, log density | at draws, gradient (rel.) | at tail points, log density | at tail points, gradient (rel.) |
| --- | ---: | ---: | ---: | ---: |
| DDM-4 | 7.3e-9 | 4.8e-8 | 3.8e-5 | 8.8e-5 |
| DDM-5 | 2.3e-13 | 1.5e-13 | 2.46 | 0.096 |

The same variant in stanli is within 5.7e-13 of CmdStan at the draws, and
within 4.9e-11 at the tail points.

The DDM-5's 2.46 is the base's error, not the new density's.
`tail_audit.R` at that point (tail76) scores all 1000 trials against the
reference:
- Stan's wiener5 is off by up to 1.10 log units on one trial, 13.6 summed.
  The worst trials have a drift of 10 against a boundary of 4.
- The new density is within 5e-4 on every trial, which is the precision of
  the 6-digit inputs that `print()` gives.

The one non-finite DDM-4 tail point is the same in the base and in every
variant: there log_prob itself is -inf.

Cost: one gradient at each of the 40 draws per block, blocks alternating
over variants and engines in random order, 21 rounds, single-threaded.
Each cell is the median microseconds per gradient, with the median ratio
to the same engine's base within round and [min, max] in brackets:

| program | engine | base | `rec` | `unrolled` |
| --- | --- | ---: | ---: | ---: |
| DDM-4 | CmdStan | 1786 | 1222: 0.86 [0.56, 0.95] | 1250: **0.86** [0.56, 1.14] |
| DDM-4 | stanli | 1212 | 27750: 24.4 [21.1, 27.2] | 1333: **1.06** [0.91, 1.40] |
| DDM-5 | CmdStan | 6625 | 1917: 0.28 [0.24, 0.37] | 1964: **0.29** [0.23, 0.36] |
| DDM-5 | stanli | 2417 | 36750: 15.6 [12.7, 17.0] | 1639: **0.67** [0.59, 0.83] |

Against CmdStan's base, stanli runs at 0.72 (DDM-4) and 0.35 (DDM-5) today,
and at 0.76 and 0.24 with `unrolled`.

AGENTS.md's recipe (`gradient_cost.R`, intercept-only, 5000 trials, 21
blocks), from 9.76 / 36.97 ms in the base:

| variant | DDM-4 | DDM-5 |
| --- | ---: | ---: |
| `direct` | 0.81 | 0.30 |
| `rec` | 0.77 | 0.30 |
| `unrolled` | 0.80 | 0.34 |

`direct` and `rec` were run before the guard.

## Applied in 0.3.6

What changed in the package:

- **`R/model_ddm.R`.**
  - `.WFPT_STAN` holds the two functions, next to `.ddm_lfpt0()`, its R
    twin. It is pasted into `.DDM_STAN_PRELUDE` as a template of its own,
    which keeps both under R's 8192-character `sprintf()` limit.
  - The DDM-4/5 call site is the one call
    `cogmod_wfpt_lpdf(t | boundary, w, v, sigmadrift)`.
  - The certificate is gone. The router stays, for the 7-parameter call
    site only.
  - `.DDM_FPT_SWITCH` is the series switch for both sides.
    `.ddm_fpt_ks()` / `.ddm_fpt_kl()` are the term counts, shared by
    `.ddm_lfpt0()` and the tests.
  - `.ddm_lfpt0()`, `.ddm_lfpt()` and `.ddm_lfpt_resp()` take the start
    point's complement exactly, and `.ddm_lfpt0()` has the `expm1()` pairing
    and the `sin()` choice. Against 0.3.5's it moves by at most 1.6e-13
    (relative) over 200,000 points.
- **`?cogmod_ddm`.** The Implementation section and the cost paragraph are
  rewritten.
- **NEWS, AGENTS.md.** Updated.
- **`benchmarks/gradient_programs.R`.** Gains `cogmod_ddm5`, so CI checks
  the 5-parameter path on every PR.
- **`test-model_ddm.R`.**
  - The code-grep expectations now assert `cogmod_wfpt_lpdf` and the
    absence of the two old calls.
  - A new test holds the 4- and 5-parameter Stan density to the R series
    at 1e-10 over 1800 cells: rescaled times from 1e-3 to 8000, and start
    points 1e-4 from either boundary.
  - Another new test checks `.DDM_FPT_SWITCH` against the term counts,
    and that m = 3 and K = 6 cover them.

Checked on the applied tree, Windows, the same day:

- **Tests.** `test-model_ddm.R` passes: no failures, and the usual model-fit
  test skipped.
- **Gradient check.** `gradient_check.R --families cogmod_ddm` passes:
  DDM-4 16/16 (5.9e-8), DDM-5 18/18 (1.7e-7). The 7-parameter program has
  only its two `.GP_KNOWN` waivers.
- **Gradient cost.** `gradient_cost.R` against a worktree of 6463109, 21
  blocks: DDM-4 11.89 to 9.62 ms (**0.81**), DDM-5 35.25 to 12.69 ms
  (**0.36**).
- **stanli.** The multilevel ddm4 and ddm5 programs from the applied tree
  are bit-identical to the prototype's `unrolled` at all 120 points of
  `program_ab.R` (log density and gradient). All of them lower to
  `OP_REGION_MAP`, with no runtime-selected loop.
- **The variability ladder** (`ladder.R`, 300 intercept-only trials, 11
  rounds). These are the multiples `?cogmod_ddm` quotes, per gradient and
  relative to the `sigmadrift`-only model, median [min, max] within round:

  | program | 0.3.5 | 0.3.6 |
  | --- | ---: | ---: |
  | classic (all three at 0) | 0.28 [0.22, 0.33] | 0.86 [0.83, 0.92] |
  | `sigmadrift` + `sigmandt` | 19.7 [15.0, 23.5] | 69.8 [65.7, 73.7] |
  | `sigmadrift` + `sigmabias` | 16.4 [13.4, 18.3] | 59.0 [54.6, 62.9] |
  | all three | 73 [62, 89] | 256 [232, 267] |

  - The 7-parameter rows run the same code in both trees. Their multiples
    grew because the `sigmadrift`-only model became 4.2x cheaper on these
    data.
  - The 0.3.5 column reproduces the help page's "about 3.5" and "about 18"
    but not its "about 30" for both ranges, which came from another data set
    and point. The help page now says 1.2 to 1.6 for `sigmadrift`: this
    ladder's 1.16, 1.32 from `gradient_cost.R`, and 1.57 on the multilevel
    model. For the ranges it says 60-70 and 250.

**When to switch back** to Stan's `wiener_lpdf()`: once a CmdStan release
carries stan-dev/math#3481 and #3442, rerun `grid.R` and `program_ab.R`
against it. `grid.R` already scores Stan's two forms (`classic`, `w5`)
beside the package's. If Stan's are as accurate and no slower on both
engines, the call in `cogmod_ddm_decision_lpdf()` is the only thing to
change. The comment there says so.

Not measured:
- **Linux.** On the cluster CmdStan builds its model methods with
  `-march=native`, and the stanli gap there has differed from Windows
  before.
- **Sampling.** No fits were run. The gradient is cheaper and more
  accurate, and wiener5's tail gradients are off by up to 84x, which can
  matter in warmup, so ESS per second may move by more or less than the
  cost ratio.
- **The DDM-7.** Its call site is untouched.
