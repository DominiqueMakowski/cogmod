# Where should `cogmod_betagate()` start? (2026-10-14)

`cogmod_inits()` covered the other two rating-scale families but refused
`cogmod_betagate()` ("nothing to offer"), so a betagate fit started wherever
`init = 0` or Stan's random `U(-2, 2)` put it. `init = 0` puts `phi`
(`softplus`) at 0.69, a U-shaped Beta (shapes `mu * phi * 2` = 0.69 at
`mu = 0.5`), and `pex` (`logit`) at 0.5, half of all responses on an
endpoint.

**Outcome.** `mu = 0.5, phi = 1, pex = 0.2, bex = 0.5` shipped in 0.3.5.
`phi = 1` is the Uniform at `mu = 0.5`, the max-entropy start
`cogmod_betadiscrete()` also uses. The start is closer to every style below
than `init = 0` is in the worst case and on average (step 1), but as with
`cogmod_lnr()`'s start the head start does not reach the sampler (step 2):
the change is for consistency, so that `cogmod_inits()` covers every family
`cogmod_priors()` does, not for speed.

## Step 1: the shortfall of the start (`shortfall.R`)

Mean log-likelihood per observation at each data set's generating values,
minus the same at the start, on 4,000 simulated trials per style: `default`
(`?rcogmod_betagate`'s defaults), `flat` (`phi = 1`), `confident`
(`phi = 8`), `high` and `low` (`mu` 0.75 and 0.25), `extreme` (`pex = 0.4`),
`ushape` (`phi = 0.5`) and `ceiling` (`mu = 0.85, pex = 0.3, bex = 0.85`,
two thirds of responses at 1). `mu` and `bex` start at 0.5 throughout.

| `phi` | `pex` | worst | mean | |
| --- | --- | --- | --- | --- |
| 1.5 | 0.5 | 0.98 | 0.68 | best worst case, but `pex` at the origin |
| 1 | 0.3 | 1.01 | 0.51 | |
| 1 | **0.2** | **1.23** | **0.50** | shipped: lowest mean of the grid |
| 0.69 | 0.5 | 1.36 | 0.71 | `init = 0` |
| 1 | 0.1 | 1.65 | 0.56 | the ceiling style sets the worst case |
| 2 | 0.1 | 1.84 | 0.67 | `cogmod_choco()`'s start |
| random | | 2.54 | | Stan's `U(-2, 2)`, median over 2,000 draws |

`pex` at 0.1, `cogmod_choco()`'s start, is the better start on all but the
two heavy-extreme styles, and costs most on the ceiling (1.65 per
observation). 0.2 gives up a little on the ordinary styles for that one.
The full grid (`phi` 0.69-6 x `pex` 0.05-0.5) is in `results/shortfall.csv`.

## Step 2: cold fits (`fits.R`)

Intercept-only model with all four parameters in `bf()` and
`cogmod_priors()`, 2,000 trials, one chain of 1,000 + 1,000 per run, 4 seeds
per arm, on three of the styles. Medians:

| style | arm | warmup leapfrog | sampling leapfrog | divergent | min ESS |
| --- | --- | --- | --- | --- | --- |
| default | `cogmod_inits()` | 6,088 | 5,076 | 0 | 1,000 |
| default | `init = 0` | 5,890 | 5,317 | 0 | 861 |
| default | random | 5,920 | 5,474 | 0 | 808 |
| ushape | `cogmod_inits()` | 5,365 | 5,407 | 0 | 816 |
| ushape | `init = 0` | 5,360 | 5,489 | 0 | 921 |
| ushape | random | 5,467 | 5,485 | 0 | 802 |
| ceiling | `cogmod_inits()` | 8,459 | 7,525 | 0 | 603 |
| ceiling | `init = 0` | 8,892 | 6,825 | 0 | 639 |
| ceiling | random | 8,694 | 6,729 | 0 | 729 |

Warmup within 5% across arms, in no consistent direction; no post-warmup
divergence in any of the 36 runs (7-19 during warmup in every arm alike);
posterior medians the same to about 0.01 on the link scale.
