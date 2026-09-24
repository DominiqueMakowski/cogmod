# Vectorising `cogmod_lnr()`: feasible, and not worth it

2026-09-24, on a Windows laptop (i7-1265U: 2 performance + 8 efficiency
cores, 12 threads), CmdStan 2.38, cmdstanr 0.9.0, brms 2.23.1. Harness in
`bench.R` (steps described at its top), the prototype likelihoods in
`proto.R`, the libm timing in `libm.R`, raw numbers in `results/`.

## The question

brms can write a custom family's likelihood as one vectorised call instead of
a loop over observations: `custom_family(loop = FALSE)`. Nothing about this is
new, and none of it is cmdstanr's - it has been in brms since 2.16 (#1084),
and cmdstanr's NEWS up to its development version has nothing on it. Is it
worth doing for the LNR, the simplest race in the package and still slow on
real data (the decision-making vignette's 4620 trials take about 75 s per
chain of 500 iterations)?

With `loop = FALSE` brms writes

```stan
target += cogmod_lnr_lpdf(Y | mu, nuone, sigmazero, sigmaone, sigmabias, ndt, poutlier, dec);
```

with a `vector` for every dpar that has a formula and a `real` for every one
that does not (fixed ones like `sigmabias = 0` included), and the whole `dec`
array (`vars = "dec"`, not `"dec[n]"`).

## Variants

| name | loop | what |
| --- | --- | --- |
| `base` | TRUE | the package as it is |
| `v0` | FALSE | the loop moved inside the function, calling the scalar lpdf |
| `v1` | FALSE | one fused loop, no nested user-function calls, `log(t)` shared |
| `v2` | FALSE | vector arithmetic over index sets (winner 0, winner 1, `t <= ndt`), `cogmod_log_Phi()` vectorised by branch, functions of real dpars computed once |
| `v3` | TRUE | not vectorised: the loser's survival as one `lognormal_lccdf()`, guarded by `cogmod_log_Phi()`'s own switch at z = 25 |

`v0`-`v2` fall back to the scalar function whenever `sigmabias > 0`.

Models, all on the vignette's data (speed_acc, participants 1-3, RT <= 2 s):

| name | formula |
| --- | --- |
| `vignette` | the vignette's: `nuone ~ Condition`, `sigmazero ~ 1`, `sigmaone ~ 1`, `sigmabias = 0`, `ndt ~ Condition` |
| `scalars` | the same with both sigmas left out of `bf()`, so brms passes them as reals |
| `hier` | `(Condition \| Participant)` on both nus, `(1 \| Participant)` on `ndt` |
| `startpoint` | `vignette` with `sigmabias ~ 1`: the general density |

## Feasible: yes, exactly

Every variant gives the base program's log density and gradient, from
CmdStan's `log_prob` method at 20 post-warmup draws of each model, on the real
data and on the same data salted with gradient_check.R's tail responses
(fastest possible trial, 3, 6, 12 s, both responses). Worst case over all
models and both data sets: 1.6e-11 on a log density of 1000-1700, and 4e-12
relative on the gradient; everything finite. `v0` is bit-identical.
`results/check.csv`.

What a real implementation would have to deal with:

- **The signature depends on the formula.** Six dpars besides `mu` can each
  be `real` or `vector`: 64 overloads. `cogmod_stanvars(f)` sees the formula
  and could emit the one it needs; `cogmod_lnr_stanvars()`, which takes no
  formula, could not.
- **brms refuses `loop = FALSE` with `weights()`, `cens()`, `trunc()` or a
  mixture.** The looped likelihood would have to stay, behind an argument.
- **Threading breaks silently.** With `threads = threading(k)` and
  `loop = FALSE`, brms slices `Y` and the dpars but passes the whole `dec`,
  so responses are misaligned with RTs and nothing errors. The lpdf would
  need `if (size(dec) != rows(Y)) reject(...)`.
- **stanc's struct-of-arrays optimisation does not reach it.**
  `stanc --O1 --debug-mem-patterns` on the `v2` program reports every
  model-block vector as array-of-structs: anything passed to a user-defined
  function is demoted. That rules out the mechanism that usually makes
  vectorised Stan code fast.
- The R side (`log_lik()`, `posterior_predict()`, priors, inits) is untouched:
  brms calls those per observation either way.

## Faster: no

Milliseconds per gradient (median of alternating blocks), and the ratio to
`base` within the same run. `vignette`/`scalars`: 21 blocks of 30 NUTS
iterations (~560 gradients); `hier`/`startpoint`: 11 blocks of 2 iterations
(~1280 and ~94 gradients; deep trees at the step size a 200-iteration warmup
adapts).

| model | base ms | v0 | v1 | v2 | v3 |
| --- | ---: | ---: | ---: | ---: | ---: |
| vignette | 3.76 | 1.00 | 0.98 | 0.98 | 0.99 |
| scalars | 2.81 | 1.01 | 0.94 | 0.90 | 0.98 |
| hier | 4.64 | 0.99 | - | 1.00 | 1.02 |
| startpoint | 14.2 | 0.96 | - | 1.04 | 0.96 |

The first two rows were run twice more: `v0`-`v2` at 1.01 / 0.99 / 0.98 and
1.01 / 0.95 / 0.90 (raw data not kept), `v3` at 1.01 and 1.02 with bases of
4.28 and 3.03 ms. Absolute times moved 15% between runs; ratios stayed within
0.03. `startpoint` is within noise at 94 gradients per block, and every
variant runs the same scalar fallback there anyway. `results/time_*.csv`.

The one gain is the 10% `v2` finds when the sigmas are reals, and it comes
from computing `log(sigma)` and `log(poutlier)` once per call instead of once
per trial - which the vectorised signature makes possible only for dpars
without a formula. The formulas people write (`sigmazero ~ 1`, random
effects) get nothing.

## Why not

`profile()` blocks around brms's linear predictors and the likelihood
(`results/profile.csv`):

| program | section | share | tape entries / trial | ns / trial |
| --- | --- | ---: | ---: | ---: |
| vignette base | linear predictors | 42% | 14 | 309 |
| vignette base | likelihood | 58% | 25 | 428 |
| scalars base | linear predictors | 19% | 10 | 100 |
| scalars base | likelihood | 81% | 25 | 424 |
| hier base | linear predictors | 48% | 20 | 405 |
| hier base | likelihood | 52% | 25 | 447 |
| vignette v2 | likelihood | | 22 | 431 |
| vignette v3 | likelihood | | 20 | 492 |

90% of the time is the forward pass, and cutting the likelihood's autodiff
tape by 12% (`v2`) or 20% (`v3`) buys nothing. What the forward pass spends
it on is the maths library. Timed with the same toolchain (`libm.R`), one
call costs:

| add | sqrt | log1p | log | exp | erfc |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 0.5 ns | 1.3 ns | 15 ns | 19 ns | 33 ns | 48 ns |

A trial's likelihood needs about ten of those (the winner's two logs, the
loser's log, `erfc` and log, the mixture's `exp` and three logarithms, and
`erfc`'s derivative in the reverse pass), and the linear predictors about
seven more (`exp` for `ndt`, `log1p(exp())` for each softplus sigma, and an
`exp` in each softplus derivative). At 15-48 ns a call, that is most of the
430 + 300 ns. Vectorising in Stan does not reduce the number of calls: each
element still goes through libm once per function. Only hoisting does, and
only for functions of reals.

This is one Windows machine, whose libm is slow (`exp` at 33 ns). On glibc
the transcendental share would be smaller and the tape's larger, so the
vectorised variants might do a little better on Linux. Not measured.

## What does help

- **Within-chain threading**, which already works with the looped family
  (brms writes `dec[nn]`, correctly sliced). Vignette model, grainsize 250,
  11 blocks (`results/threads.csv`): 1.32x faster at 2 threads and 1.72x at
  4, but 0.77x - dearer - at 1, the cost of `reduce_sum`'s slicing. Worth it
  only with idle cores, and on this hybrid CPU the extra threads land partly
  on efficiency cores, so a desktop may scale better.
- **Leaving intercept-only dpars out of `bf()`.** `scalars` against
  `vignette` is the same likelihood with the sigmas as reals: 25-29% cheaper
  per gradient in all three runs (2.81 vs 3.76, 2.76 vs 3.70, 3.03 vs 4.28
  ms). With `sigmazero ~ 1`, brms builds a 4620-vector of one intercept and
  puts every element through the softplus - the linear-predictor profile
  above, 309 against 100 ns per trial. The priors are then natural-scale rows
  (`cogmod_priors()` has them); whether the geometry changes was not
  measured.

## Decision

Decided 2026-09-24: not implemented, for the LNR or any other family. The
package keeps the looped likelihoods; nothing outside this folder changed.
AGENTS.md ("Tried and decided against") points here.

The reasoning: a formula-specific code generator, a second code path for
weights, censoring and mixtures, and a guard against a silent threading bug
buy 0-2% on the models people fit. The other families were not measured, but
the reasons carry over: every family's likelihood is a user-defined function,
which is what keeps it array-of-structs, and every RT and choice family pays
for the same per-trial outlier mixture. `v3` is dropped too: it bought
nothing, and it uses `lognormal_lccdf()`, which AGENTS.md rules out for
normal tails (safe here only because of the z = 25 guard).

Reopen it if a Linux measurement (where libm is cheaper and the tape a larger
share) shows more than the 10% ceiling seen here, or if stanc starts applying
struct-of-arrays to user-defined function arguments.
