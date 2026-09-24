# Probit or logit for `cogmod_gaussbit()`?

Run 2026-09-24 on a 16-core Windows laptop, CmdStan 2.38, brms 2.23.1, with
`bench.R` at its defaults: 3 simulated datasets per setting, 4 chains x
(1000 + 1000) iterations, brms' default random inits for every model. Raw
numbers in `results/fits.csv` (one row per fit), `results/gradients.csv` and
`results/summary.csv` (medians over the three datasets).

## The models

| name | what |
| --- | --- |
| `gp_rho` | `cogmod_gaussbit()`, `rho ~ 1` |
| `gl_rho` | its logit twin: same Gaussian copula, logistic choice marginal, `q = Phi^-1(logistic(mudec))` per trial (defined in `bench.R` only) |
| `gp_0` | `cogmod_gaussbit()`, `rho = 0` |
| `mv_probit` | native `gaussian()` + `bernoulli("probit")`; the same model as `gp_0` |
| `mv_logit` | native `gaussian()` + `bernoulli("logit")`: the default analysis |

Data: 3600 trials, two conditions, simulated from `cogmod_gaussbit()` at
`rho = 0.3`, about 10% errors. `mixed` has 30 participants x 120 trials with
random intercepts and condition slopes on the RT mean and on the choice, and
random intercepts on `sigma`; `fixed` has the same trials and no participant
structure.

## Results (medians over 3 datasets)

| setting | model | leapfrog / iter | leapfrog / min-ESS | us / gradient | CPU ms / min-ESS |
| --- | --- | ---: | ---: | ---: | ---: |
| fixed | gp_rho | 6.4 | 5.6 | 1815 | 10.2 |
| fixed | gl_rho | 6.4 | 6.0 | 4623 | 27.8 |
| fixed | gp_0 | 6.0 | 6.3 | 1564 | 9.8 |
| fixed | mv_probit | 6.0 | 6.3 | 897 | 5.6 |
| fixed | mv_logit | 6.0 | 6.9 | 628 | 4.3 |
| mixed | gp_rho | 32.0 | 303 | 1963 | 595 |
| mixed | gl_rho | 31.8 | 349 | 4673 | 1633 |
| mixed | gp_0 | 34.1 | 373 | 1782 | 664 |
| mixed | mv_probit | 31.7 | 245 | 1056 | 259 |
| mixed | mv_logit | 31.4 | 304 | 895 | 272 |

No divergences in any of the 30 fits; max Rhat 1.03.

## What it says

- **The link does not change the geometry.** Leapfrogs per iteration are the
  same for every model within a setting (about 6 without random effects, 31-34
  with them), and leapfrogs per effective sample differ by no more than
  models that are *identical* do: `gp_0` and `mv_probit` are the same
  likelihood with the same priors and read 373 and 245 in the mixed setting.
  With three datasets that is the noise floor, and probit against logit sits
  inside it in both settings.
- **The difference is the cost of a gradient.** The logit twin is 2.4-2.6x
  the probit per gradient in both settings, from the `std_normal_log_qf()`
  that maps the logistic marginal onto the copula's normal scale, and so
  2.7x the CPU time per effective sample. Probit is the cheaper
  implementation of the correlated model, for that reason alone.
- **Among the native brms models it goes the other way:** logit is 15-30%
  cheaper per gradient than probit (`bernoulli_logit` against
  `bernoulli(Phi())`) - ahead per effective sample without random effects
  (4.3 against 5.6 ms), level with them (272 against 259).
- **The custom family costs about 1.7x the native model per gradient**
  (`gp_0` against `mv_probit`, the same model): the per-observation loop brms
  writes for a custom family. `custom_family(loop = FALSE)` with a vectorised
  lpdf is the lever, untried.
- **In the mixed model the bottleneck is the RT half, not the choice.** The
  worst-mixing parameter is the population RT intercept in 14 of the 15 mixed
  fits (the choice's slope SD in the other), for every model alike: with 120 trials
  per participant the likelihood pins each participant's mean RT tightly, which
  is the regime where brms' non-centred random effects mix worst. No choice of
  link touches that.
