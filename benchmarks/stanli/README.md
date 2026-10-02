# stanli as a backend for cogmod: 3-4x faster, with the numerics written as selects

2026-10-02, on a Windows laptop (i7-1265U: 2 performance + 8 efficiency
cores, 12 threads), CmdStan 2.38, cmdstanr 0.9.0, brms 2.23.1, stanli 0.18.1.
Harness in `bench.R` (steps described at its top), the cost bisection in
`bisect.R`, a stanli-only reproducer for upstream in `repro_issue.R`, the rewritten likelihoods in
`lnr_select_functions.stan` and `lnr_branchfree_functions.stan`, raw numbers
in `results/`. Reported upstream as
[seantalts/stanli#422](https://github.com/seantalts/stanli/issues/422).

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

## Short answer

- **As brms writes it, stanli refuses the program** over one call:
  `log_mix()` inside an `if` that depends on a parameter. Written out as
  `log_sum_exp()`, the package's density compiles, is exact, and is **6.5x
  slower** per gradient than CmdStan. stanli compiles an `if` on a parameter
  into an interpreted "runtime-control region", and the LNR's density is
  made of them.
- **Written with no control flow, and the same numerics, it is 3.1x
  faster** per gradient and 3.7x per chain than CmdStan, with the same
  posterior. The same rewritten source in CmdStan is 1.36x slower than the
  original, so it would be a stanli-only source. Each branch becomes a *select*: every arm is evaluated on an
  input clamped into its own domain, so it is finite everywhere, and
  `step()` weights pick one. That keeps `cogmod_log_Phi()`'s far-tail
  series: where a naive branch-free copy's gradient goes non-finite, the
  select's matches CmdStan to 1.4e-15.
- Nothing in the package changed. Shipping this would mean a second,
  stanli-shaped copy of each family's Stan code, for a backend brms has not
  merged; only the LNR at `sigmabias = 0` is done.

## What stanli refuses

From the source, the region compiler appeared to take only arithmetic,
`exp`, `log` and a few more, which would have ruled out `erfc`, `log1p`,
`log1m_exp`, `lognormal_lpdf()` and user-function calls. Small probe programs
say otherwise (`repro_issue.R`, `results/repro_issue.txt`): all of those
compile inside a parameter-dependent `if`. What fails:

| construct | result |
| --- | --- |
| `log_mix()` inside `if` on a parameter | `runtime-control region: function log_mix` |
| `log_mix()` with no branch | compiles |
| `log_sum_exp()` written out, inside the same `if` | compiles |
| comparison as a value, `log(y > mu)`, branch or not | `unsupported function Greater__` |
| the same comparison in a ternary or an `if`; `step()` | compiles |

Compiling is not the end of it: a parameter-dependent `if` is slow even when
never taken. `repro_issue.R` wraps a `log Phi` in one at N = 5000: same
value, same gradient, 2.7x the cost.

## Variants

| name | what |
| --- | --- |
| `orig` | the program `brm()` writes, unchanged (CmdStan only) |
| `lse` | `orig` with `log_mix(poutlier, lp_out, lp_dec)` written out as `log_sum_exp(log(poutlier) + lp_out, log1m(poutlier) + lp_dec)`; nothing else changes |
| `bf` | `orig`'s functions block replaced by `lnr_branchfree_functions.stan`: `log Phi` as `log(0.5 * erfc())` with no tail series, winner and loser picked by `dec` (data), `fmax(Y - ndt, 1e-300)` for the `t <= ndt` return. `sigmabias = 0` only |
| `sel` | `orig`'s functions block replaced by `lnr_select_functions.stan`: as `bf`, but `log Phi` keeps the series below -25 as a select, and `t <= ndt` is a mask inside `log_mix()`. No checks. `sigmabias = 0` only |

All three keep the parameters, priors and generated quantities, so
`rename_pars()` and the rest of brms see the program they expect.

## The select pattern

A branch `if (x < -25) lo(x) else hi(x)` becomes

```stan
real xl = fmin(x, -25);       // lo's arm, clamped into its domain
real xh = fmax(x, -25);       // hi's arm, likewise
real w = step(-25 - x);       // 1 below -25
return w * lo(xl) + (1 - w) * hi(xh);
```

Both arms are evaluated, but each on an input where it is finite, so the
arm not picked has a finite value and a weight of 0. Its input is the clamp
edge, a constant, so it passes back an exact zero, never `0 * inf`. The
picked arm gets the whole gradient. Where the two branches meet continuously,
as `cogmod_log_Phi()`'s do at -25 (to the last bit), the select reproduces
the branch exactly.

What it took for stanli to keep it fast, measured on the LNR (`bisect.R`,
`results/bisect.csv`: median µs per gradient over 11 alternating blocks;
every row has the same log density):

| `log Phi` | `t <= ndt` | mixture | checks | µs |
| --- | --- | --- | --- | ---: |
| erfc only (`bf`) | clamp to 1e-300 | `log_mix()` | none | 900 |
| select | clamp to 1e-300 | `log_mix()` | none | 1000 |
| select | mask: `lp_dec + (step(Y - ndt) - 1) * 1e300` | `log_mix()` | none | **1100** (`sel`) |
| select | mask | `log_mix()` | data only (`dec`, `Y > 0`) | 4800 |
| erfc only | clamp | `log_sum_exp()` + `log1m()` | none | 4700 |
| erfc only | blend: `w * log_mix(...) + (1 - w) * (...)` | `log_mix()` | none | 5500 |
| select | blend | `log_mix()` | none | 10900 |

So the select itself is nearly free. The `t <= ndt` switch has to live
inside the single `log_mix()` call, as a mask that pushes the decision
component to -1e300 where `log_mix()` gives it an exact zero weight and
adjoint, not around it. Even early returns on data alone cost 4x, although
stanli has the data when it builds the model. Alone, in a minimal program,
`log_sum_exp()` costs 1.2x and the data check nothing
(`results/repro_issue.txt`), so these cliffs depend on the rest of the
program.

The mask relies on the decision component being finite at the clamped time,
which the select `log Phi` guarantees for any parameter values. `bf`'s clamp
alone relies on it being below about -1e3 at t = 1e-300 s, which holds for
any sigma below ~20 on the log scale. `sel` drops the checks: the parameter
ones are guaranteed by brms's links and bounds, the data ones by the data
block and `.cogmod_checkdata()`. One small difference remains: above
x = 0 `sel` takes `log(0.5 * erfc(-x / sqrt2))` rather than the package's
`log1p(-0.5 * erfc(x / sqrt2))`. That loses relative accuracy only in a
value within 1e-16 of 0, and gives the same gradient.

## Exact: yes

At the init and five points around it (`results/check.csv`), the log
density of every variant, in CmdStan and in stanli, equals `orig`'s to every
printed digit (10 significant). Gradients, relative to `orig`'s: `sel` in
stanli 2e-15 at most, `bf` 1e-14, `lse` 1.4e-14.

The tail (`results/tail.csv`): both sigmas at softplus(-3) = 0.049 and one
trial replaced by a 12 s response, which puts the loser's survival at
x = -65. `orig` is finite there, by design. In stanli:

| variant | lp | gradient |
| --- | --- | --- |
| `lse` | = `orig` | finite, 1.0e-14 from `orig` |
| `bf` | = `orig` | **not finite** |
| `sel` | = `orig` | finite, 1.4e-15 from `orig` |

`bf` is the 0.3.3 LNR bug: the value survives, the gradient does not.

## Cost per gradient

R loops over 200 points, a different one per call, 21 blocks alternating the
six programs (`results/time.csv`); one core each (CPU time / elapsed = 1.00).

| program | median ms | IQR | / CmdStan `orig` |
| --- | ---: | ---: | ---: |
| CmdStan `orig` | 4.20 | 3.45-4.45 | 1 |
| CmdStan `bf` | 4.05 | 3.40-4.50 | 0.96 |
| CmdStan `sel` | 5.70 | 4.45-6.15 | **1.36** |
| stanli `lse` | 27.8 | 24.7-28.2 | 6.6 |
| stanli `bf` | 0.95 | 0.85-1.15 | 0.23 |
| stanli `sel` | 1.25 | 1.10-1.40 | **0.30** |

Earlier runs of the same harness gave the same ratios within 0.02
(`bf` 0.23-0.25, `lse` 6.5-6.6, `sel` 0.30-0.32); absolute times moved 20%
between runs. The block-by-block ratio of CmdStan `sel` to `orig` ranged
1.27-1.45 (10th-90th percentile), median 1.36.

**`sel` is not free in CmdStan.** Its selects evaluate every arm: each
`log Phi` computes the asymptotic series (`inv_square`, a polynomial, two
logarithms) as well as `erfc`, where `orig` computes one or the other. On a
maths library as slow as this machine's, that costs 36% more per gradient,
above `gradient_cost.R`'s 1.3 threshold. `bf`, which does no extra work, is
free (0.96). So one Stan source cannot serve both backends at full speed:
the branches are what CmdStan wants and what stanli cannot afford.

Why stanli's graph path beats CmdStan by 3-4x was not investigated. CmdStan
on this machine is bound by the maths library (`exp` at 33 ns, `erfc` at
48 ns: `benchmarks/lnr_vectorize/README.md`). stanli ships its own prebuilt
kernels, which may simply be cheaper on Windows. On Linux, where glibc's
libm is faster, the gap is probably smaller. Not measured.

## A real fit

4 chains x (500 warmup + 500 sampling), in parallel: `brm(backend =
"cmdstanr")` on `orig`, inits from `cogmod_inits()`, against the light route on
`sel`. stanli's `sample_model()` takes no list of per-chain inits ("'list'
object cannot be coerced to type 'double'"), so its chains start from its
own random inits (radius 2). `results/fit_chains.csv`,
`results/fit_pars.csv`.

| | CmdStan `orig` | stanli `sel` |
| --- | ---: | ---: |
| seconds per chain, warmup + sampling (mean) | 112.3 | **30.4** |
| leapfrog steps, all chains | 37672 | 40184 |
| ms per leapfrog step | 11.9 | 3.0 |
| bulk ESS, min-max over 9 parameters | 812-2119 | 968-1680 |
| R-hat, max | 1.004 | 1.006 |
| divergences | 0 | 0 |
| wall time, `brm()` to a brmsfit | 138 s | 37 s |

Posterior means differ by at most 0.05 posterior SDs. `loo_compare()`: elpd
difference 0.0 (SE 0.1). The earlier `bf` fit, same harness, came out the
same (29.9 against 126.7 s per chain, 37456 leapfrog steps), before the
harness moved to `sel`.

A gradient costs about 3x more inside a 4-chain fit than in the
single-process loop above, for both backends alike (11.9 against 4.2 ms,
3.0 against 1.35). That is four chains on a hybrid CPU whose extra cores are
efficiency cores, not anything stanli-specific. The CmdStan wall time
excludes compiling, which cmdstanr had cached; stanli built in about 2 s.

## What the light route needs

- `brms:::rename_pars()` is internal. After it, `summary()`, `fixef()`,
  `loo()` and `loo_compare()` worked on the result. The comment it came from
  says refit-based functions (`reloo`, moment matching) do not.
- brms#1911 itself may change all this. Paul Bürkner's reply there plans to
  stop converting fits to `stanfit`, so the post-processing would be per
  backend.
- `cogmod_inits()` cannot be passed in (above). Random inits were fine here.

## For cogmod

`sel` shows that the package's numerics survive as selects at stanli's fast
speed, so the blocker is maintenance, not maths. To ship it:

- **A second Stan source per family**, emitted when the backend is stanli
  (`cogmod_stanvars()` could switch on it). Not one shared source: `sel`
  costs CmdStan 1.36x per gradient (above). Every branch in every family
  would need the select form: the `sigmabias == 0` dispatch and the `A > 0`
  path of the start-point LNR and LogNormal (several regimes each), the
  DDM, the LBAs. Each select evaluates all of its arms, so a function with
  many regimes pays for all of them.
- **A second set of gradient checks.** `benchmarks/gradient_check.R` runs
  through CmdStan. The select sources compile in CmdStan too (`bf`'s did,
  with the same cost as `orig`), so they could be checked there with no
  extra machinery.
- **brms has to ship the backend**, or users take the light route by hand.

Worth doing if brms merges a stanli backend, or if fitting without a C++
toolchain becomes something users ask for. Before then the cost is a second
copy of every family's Stan code for an unmerged backend.

## Upstream

[seantalts/stanli#422](https://github.com/seantalts/stanli/issues/422)
(2026-10-02), with `repro_issue.R` (stanli only, no cogmod) as its
reproducer, reports the `log_mix()` gap, comparisons as values and the cost
of a never-taken branch, with the LNR numbers as the motivating case. Its
main ask is for small pure branches to be if-converted to a select, as above,
so that code written the ordinary way would not need `sel`'s rewriting; then
`log_mix()` in regions. It offers testing against this harness and a PR for
the `log_mix()` part (PR #223, which added ops to parameter-dependent
regions, is the pattern to follow). Related: #374 (brms sampling
performance).

## Status

Exploratory; nothing outside this folder and `.gitignore` changed. Reopen
when brms ships a stanli backend or stanli makes ordinary branches cheap;
then rerun `bench.R` and extend `sel` to the other families, or drop it if
stanli no longer needs it.
