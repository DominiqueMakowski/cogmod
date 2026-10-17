# A cheaper routing test for the DDM's classic Wiener path

Proposed 2026-10-16, implemented and measured 2026-10-17 (0.3.5). With all
three variabilities fixed at 0 in `bf()` ("DDM-4"), a gradient now costs
0.77x what it did in CmdStan and 0.83x in stanli, and every log density and
every partial derivative is what it was, bit for bit.

| file | what |
| --- | --- |
| `certificate_check.R` | the proof below, checked numerically in pure R against a transcription of the router |
| `grid_ab.R` | function-level A/B: `cogmod_ddm_decision_lpdf()` from a base tree against the current one, value and partials per trial at 200,000 points |
| `program_ab.R` | program-level A/B on `multilevel.R`'s ddm4 program: exactness at posterior draws and tail points, then the cost per gradient in CmdStan and stanli |
| `ddm4_harness.R` | runs `benchmarks/gradient_check.R` or `gradient_cost.R` on a DDM-4 program, which the CI programs never build |

## Where the cost was

`cogmod_ddm_decision_lpdf()` (R/model_ddm.R, in `.DDM_STAN_PRELUDE`)
chooses between Stan's implementations of the Wiener density on every
trial. On the DDM-4 it takes the classic 4-parameter `wiener_lpdf()` only
where that one is sound, and the 5-parameter one ("wiener5") elsewhere:

```stan
if (sigmadrift == 0 && t > 0.01 * square(boundary)
    && cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600) {
  return wiener_lpdf(y | boundary, tau0, w, v);              // classic
}
return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift);    // wiener5, no guard needed
```

The router, `cogmod_ddm_log_density_scale()`, estimates the log density so
that the classic form is kept out of the region where it underflows. There
it returns `-inf` with NaN partials, and one such trial makes the whole
gradient NaN, even under a zero mixture weight (0 * NaN). The comment above
the function explains why.

The router's value is only compared with -600, but Stan records it on the
autodiff tape like any other expression: about 25 nodes per trial,
including `log()` x 4, `sin()`, `sqrt()` and a `pow()`. Its adjoint is
always zero, so it changes nothing in the gradient and costs the full
price. `benchmarks/stanli/probe_bottlenecks.R` (Windows, stanli 0.20.0,
DDM-4 with participants, N = 1000, 40 points near the inits) put it at up
to 22% of CmdStan's gradient and 15% of stanli's: 2000 against 1556 us and
2071 against 1750 us with the router dropped from the condition.

Only the DDM-4 pays this. With `sigmadrift` estimated the `&&`
short-circuits before the router. The 7-parameter call site also calls the
router, but next to that density's 100-274 us per trial it is noise, and it
is left alone.

## The change: a certificate in front of the router

For nearly every trial the router's value sits hundreds of log units above
-600, which a few multiplications and comparisons can prove. That is tested
first, and the router runs only when the cheap test cannot vouch for the
trial:

```stan
if (sigmadrift == 0) {
  real a2 = square(boundary);
  if (t > 0.01 * a2
      && ((boundary < 100 && w > 1e-3 && w < 1 - 1e-3 && t < 50 * a2
           && v * (v * t + 2 * boundary * w) < 600)
          || cogmod_ddm_log_density_scale(t, v, boundary, w, 0, 0, 0) > -600)) {
    return wiener_lpdf(y | boundary, tau0, w, v);
  }
}
return wiener_lpdf(y | boundary, tau0, w, v, sigmadrift);
```

- Comparisons record nothing on the tape; only arithmetic on parameters
  does. The certificate records about 6 cheap nodes (one `square`, a few
  products and sums) instead of the router's 25 with seven transcendental
  functions. `||` short-circuits, so the router runs only for trials the
  certificate rejects.
- `0.01 * a2` is bit for bit the `0.01 * square(boundary)` of before, so
  the small-time test is unchanged.
- The same code serves both boundaries: the caller passes (`-mu`,
  `1 - bias`) for the lower one, and the certificate is symmetric in w.

**Every trial takes the branch it took before.** Where the certificate
holds, the router's value is provably above -600 (next section), so the
old code took the classic branch too. Where it fails, the router runs as
before. The branch is the same, so the log density is the same bit for
bit, and so is the gradient, because the nodes removed only ever received
zero adjoints. That is what makes the change safe without any new argument
about the density itself.

The prelude is a `sprintf()` template, and this change took it past R's
limit of 8192 characters for a format string (7581 before, 8377 after), so
it is now one template per function, joined by `paste0()`. The generated
Stan code is unchanged outside the routing block.

## The proof

At sv = sw = st0 = 0 the router is

    est = -2 log a + min(S, L) + D,   tau = t / a^2,
    S = log w - w^2 / (2 tau) - 0.5 log(2 pi) - 1.5 log tau
    L = log pi - pi^2 tau / 2 + log sin(pi w)
    D = -v^2 t / 2 - v a w = -v (v t + 2 a w) / 2

Inside the certificate w is in (1e-3, 1 - 1e-3), so the router's clamps of
w to [1e-12, 1 - 1e-12] and of its interval ends do nothing, and t > tau0 =
1e-10 makes its `fmax(t, 1e-12)` the identity. Then, for tau in (0.01, T),
a < A, w in (w0, 1 - w0) and D > -K:

- S >= log w0 - 1/(2 * 0.01) - 0.5 log(2 pi) - 1.5 log T, using w <= 1
  and tau > 0.01;
- L >= log pi - pi^2 T / 2 + log sin(pi w0);
- -2 log a > -2 log A.

With T = 50, A = 100, w0 = 1e-3 and K = 300: S >= -63.7, L >= -251.4, so
est > -9.2 - 251.4 - 300 = **-560.6 > -600**, with 39 log units to spare
for rounding. The knobs trade against each other through
`-2 log A + min(S_lb(w0, T), L_lb(w0, T)) - K > -600`, and the large-time
term dominates (`-pi^2 T / 2` is -247 at T = 50). If a bound is changed,
the constraint has to be recomputed, and `certificate_check.R` updated.

`certificate_check.R` checks this numerically. Over 2e6 log-uniform
points, 605002 are certified, none has a router value at or below -600, and
the lowest is -532. Every trial like those of a DDM-4 posterior that
reaches the router is certified, so in the sampling phase the router should
essentially never run.

## Verified

2026-10-17, Windows (i7-1265U), CmdStan 2.38.0, stanli 0.21.0, against a
worktree of the base (`git worktree add --detach <tmp>/base HEAD`).

- **Function level** (`grid_ab.R --base <tmp>/base`): 200,000 points in 40
  batches of 5000, 30 batches log-uniform over wide ranges and 10 like a
  DDM-4 posterior's trials. Each trial's t, v, boundary and w are
  parameters of their own, so one gradient gives every trial's partials,
  including the one on t (which is `Y - ndt` in a model). 95152 points are
  certified; of the rest, 39537 go to the classic form through the router,
  12742 to wiener5, and 52569 fall below the small-time cut. Log density,
  gradient and each trial's value are identical to the base in all 40
  batches.
- **Program level** (`program_ab.R --base <tmp>/base`): the ddm4 program
  of `benchmarks/stanli/multilevel.R` (10 participants x 2 conditions x 50
  trials), at 40 posterior draws of a short CmdStan fit and 80 points moved
  from them by N(0, 1.5) and N(0, 3) on the unconstrained scale. Log
  density and gradient identical at all 120, in CmdStan and in stanli.
- **Gradient check** (`ddm4_harness.R benchmarks/gradient_check.R
  --families cogmod_ddm`): 16 of 16 points ok, worst relative error against
  central differences 5.9e-8. (Its note "drop the `.GP_KNOWN` entry" is
  about the 7-parameter model, which this program does not build.)
- `testthat::test_local(filter = 'model_ddm')`: 710 passed, 1 skipped (the
  slow fit).

## Cost

Same machine and versions. `program_ab.R` times one gradient at each of the
40 posterior draws per block, 21 rounds of blocks alternating over variants
and engines in random order, single-threaded. It ran twice: first with the
certificate spliced into the base by text (`cert`), then with the change
as implemented (`new`). Ratio to base within round, median [min, max]:

| variant | CmdStan, run 1 | CmdStan, run 2 | stanli, run 1 | stanli, run 2 |
| --- | ---: | ---: | ---: | ---: |
| base (us per gradient) | 1813 | 2350 | 1386 | 1857 |
| certificate | **0.77** [0.73, 1.00] | **0.77** [0.64, 0.91] | **0.83** [0.81, 0.95] | **0.83** [0.69, 0.99] |
| router dropped (unsafe; the bound) | 0.75 [0.70, 0.83] | 0.76 [0.62, 0.92] | 0.78 [0.75, 0.82] | 0.77 [0.64, 0.93] |
| wiener5 for every trial | 2.71 [2.62, 2.86] | 2.80 [2.23, 3.50] | 1.20 [1.18, 1.26] | 1.17 [1.04, 1.38] |

The laptop ran warmer in run 2 (every absolute time is about 30% higher),
and the ratios did not move. The certificate takes almost all of what
dropping the router would, without its cost: with the router dropped, 23 of
the 120 points have a non-finite gradient (2 in the base) and the log
density moves by up to 355.

AGENTS.md's recipe, on a DDM-4 program (`ddm4_harness.R` running
`gradient_cost.R`: intercept-only, 5000 simulated trials, median of 21
alternating blocks): 14.61 ms per gradient in the base, 11.25 ms with the
change, **ratio 0.77**.

Not measured: Linux, where the 0.21.0 rerun put the DDM-4 at 1950 us
(CmdStan) and 2361 us (stanli) per gradient
(`benchmarks/stanli/README.md`, section "0.21.0"). `multilevel.R grad
--model ddm4 --build matched` on the cluster would give it. On this laptop
stanli is already the faster engine on this program (0.76-0.79x CmdStan's
gradient in the base, 0.82-0.85x with the change), where the Linux nodes
gave 1.22x.

Since the classic model got cheaper, `?cogmod_ddm`'s cost of estimating
`sigmadrift` relative to it went from "about 2.8 times" to "about 3.5
times": wiener5 against the certificate is 3.5x and 3.6x in the two runs
above, where wiener5 against the base is 2.7x and 2.8x, which is what the
old figure measured.

## Tried and decided against: wiener5 for every trial

The call-site comment justifies the classic form with 2.0 against 5.4 us
per trial (2026-09-18), and the 0.21.0 rerun hinted that wiener5 had
caught up (the DDM-5 cost CmdStan about what the DDM-4 did on Linux). On
the same program, routing every DDM-4 trial through the 5-argument
`wiener_lpdf()` at sv = 0 costs 2.7-2.8x the base in CmdStan and 1.2x in
stanli (table above), and it is not the same density in the tails (next
section). The classic form stays.

## Side finding: the two Stan forms against the R series

`grid_ab.R`'s 150,000 wide-range points, scored by both Stan forms and by
the package's R series `.ddm_lfpt()` (log space, 1e-12 absolute on the
rescaled density; Stan's `wiener_lpdf()` is the upper boundary, so the
reference is `.ddm_lfpt(t, -v, a, 1 - w)`):

- **The classic form**, where the router sends trials to it (84690
  points): within 7e-12 of the R series below tau = 0.1 and 4e-13 above
  tau = 1. In between, 762 points are off by more than 1e-6, at most 4.7e-3,
  all at tau 0.46-0.96 and a start point within 0.022 of the far boundary
  (most within 1e-3).
- **wiener5** disagrees with the R series by more than 1e-3 at 3370 points,
  up to 1.33 log units, at tau 0.10-0.56, and only where the log density is
  below -26. That fits an error controlled in absolute terms on the density.
  Near the small-time cut (tau 0.01-0.1) both forms agree with the R series
  to 1e-11, so there is no step where the two paths meet there. On trials
  like a posterior's the two Stan forms agree to 1.6e-8 in value and 2e-7
  in the gradient, not the 1e-13 the call-site comment says.

This touches the DDM-4 only in the far tail it already sends to wiener5,
where the outlier component dominates the mixture. Models that estimate
`sigmadrift` use wiener5 everywhere, so a trial below -26 there can be off by
up to a log unit before the mixture. Not acted on; a possible Stan Math
report.

## Related, not this change

- **The CI programs still miss this path.**
  `benchmarks/gradient_programs.R` estimates every DDM dpar (`.GP_FIXED`
  fixes nothing for `cogmod_ddm`), so `gradient_check.R` and
  `gradient_cost.R` reach only the 7-parameter call site, never the classic
  one, and the most common DDM configuration has no gradient check in CI.
  `ddm4_harness.R` runs both scripts on it locally; adding a DDM-4 program to
  the CI harness is a separate decision.
- **On stanli's side**, values used only in a branch condition could be
  evaluated without the tape, which CmdStan cannot do from Stan code. It is
  on the checklist of
  [seantalts/stanli#455](https://github.com/seantalts/stanli/issues/455). With
  this change the DDM-4 no longer depends on it.
