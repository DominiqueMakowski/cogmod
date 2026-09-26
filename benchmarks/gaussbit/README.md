# `cogmod_gaussbit()`: built, measured, dropped

2026-09-24, on `dev`, before any release. Nothing of it is left in the
package; this directory is the record.

## What it was

A descriptive joint model of choice and RT, meant as the baseline the
sequential sampling models would be compared against: a **Gaussian** RT and a
**probit** choice, correlated `tanh(rho)` at the trial level (Fisher-z `rho`,
identity link). Whatever `rho` is, the RT is exactly `Normal(mu, sigma)` and
`P(dec = 1)` exactly `pnorm(mudec)`; `rho` sets how the choice depends on the
RT (slow errors against fast ones). At `rho = 0` it is `gaussian()` +
`bernoulli("probit")`, likelihood and priors included. It took the same
`rt | dec(response)` data as the LNR, RDM, DDM and LBA, so `loo_compare()` could
put it next to them. It had no `ndt` and no outlier component, so it was a
standalone family outside `.CHOICE`.

The design choices that were settled while it existed are kept too. Probit
rather than logit, because a logit marginal on a Gaussian copula needs a
`std_normal_log_qf()` per trial, 2.4-2.6x the gradient for no change in
geometry (`links/`). The choice dpar was called `mudec`, and `rho` was on
Fisher z because brms has no link onto (-1, 1).

## Why it was dropped

Everything below is on the `decision_making` vignette's data (speed_acc,
participants 1-3, RT <= 2 s, 4620 trials), with every parameter
`~ Condition`. Maximum likelihood, from `decomposition.R`
(`results/decomposition.csv`):

| model | logLik | k |
| --- | ---: | ---: |
| gaussbit, `rho = 0` (= `gaussian()` + `bernoulli("probit")`) | 34.94 | 6 |
| gaussbit, `rho` free | 35.88 | 8 |
| native brms: `gaussian()` + probit `Error ~ Condition * RT` | **35.88** | 8 |
| log-RT + probit, `rho` free | 905.17 | 8 |
| log-RT + probit + `ndt ~ Condition` + `poutlier` ("lognorbit") | 1318.29 | 11 |
| LNR, the vignette's formula | 1325.74 | 9 |

1. **With `rho` free, the likelihood is one brms already fits natively.** The
   probit of the choice given the RT is linear in the RT, so the model is
   `Error ~ Condition * RT` under another parameterisation, and the two agree
   to two decimals above. The native model is 1.7x cheaper per gradient (the
   custom family's per-observation loop; `links/`). It also enters
   `loo_compare()` against the `dec()` families, because the loo yhash of
   `bf(RT ~ .) + bf(Error ~ .)` equals that of `RT | dec(Error)` when RT is
   listed first. That leaves only two things unique to the family:
   - `mudec` stays on the marginal scale, which differs from the native
     model's conditional coefficients only once there are random effects;
   - `posterior_predict()` draws RT and choice jointly rather than the choice
     from the observed RT.
2. **`rho` bought almost nothing:** 0.9 log-likelihood units for two
   parameters (Fisher z 0.056 and 0.006). A probit that is linear in the RT
   can describe only a monotone conditional accuracy function.
3. **As a baseline, the Gaussian mostly measures skew.** In the vignette's
   Bayesian fits gaussbit's `elpd_loo` was 23 against 1316 for the LNR and
   1210-1218 for the RDM and DDM. Of the roughly 1290-unit log-likelihood gap
   above, taking the log of the RT closes about 870, and a shift plus the outlier
   component (which the fastest speed trial, 0.071 s, makes indispensable)
   close about 410. That leaves about 7 for the race itself, nearly all of it
   in the speed condition's errors. Observed error-minus-correct mean RTs are
   +34 ms and +2 ms. The LNR implies +17 and -15 ms, and lognorbit +8 and
   -11 ms.
4. **The informative version is barely cheaper than the LNR.** `lognorbit_cost.R`
   fits a prototype of the shifted lognormal + probit + outliers model and the
   vignette's LNR, with the same data and the same `ndt` / `poutlier` priors,
   4 chains x (500 + 500), three seeds each in alternating order, and then
   times `grad_log_prob()` in 21 alternating blocks. It was run twice:

   | | LNR | lognorbit | LNR / lognorbit |
   | --- | ---: | ---: | ---: |
   | us per gradient, run 1 | 4937 | 4292 | 1.15 |
   | us per gradient, run 2 (`results/lognorbit_gradients.csv`) | 3564 | 2887 | 1.23 |
   | CPU ms per min-ESS, median of 3, run 1 | 287 | 283 | 1.01 |
   | CPU ms per min-ESS, median of 3, run 2 (`results/lognorbit_fits.csv`) | 252 | 229 | 1.10 |
   | wall s per fit, run 2 | 142-145 | 124-128 | |

   There were no divergences, max Rhat was 1.007, and the worst-mixing
   parameter was `ndt` in 5 of the 6 fits. At `sigmabias = 0` one LNR trial
   is a lognormal density for the winner and a normal tail for the loser,
   while a lognorbit trial is a lognormal density and a normal tail plus a
   `cosh` and a `sinh`: the same arithmetic. So the LNR, with two fewer
   parameters and a slightly better fit, costs 0-10% more per effective
   sample. (The "LNR = 4.6x gaussbit" figure that circulated while the family
   existed came from the intercept-only programs of `gradient_check.R`; it
   does not hold for a real model.)
5. **Upkeep.** Being outside the registries, the family's name was hardcoded at
   about 17 sites: `.checkdata_class()`, the `cens()` refusal,
   `cogmod_stanvars()`, `.PRIORS_PLAIN` (with a `free_slopes` option in
   `.priors_dpars()` that only it used), `.INIT_PLAIN`, `helper-stan.R`,
   `_pkgdown.yml`. It also carried a vignette section, a panel of the
   posterior predictive figure and a cached model to refit.

**For a "default analysis" baseline in a vignette**, fit the native
multivariate model: `bf(RT ~ Condition) + bf(Error ~ Condition, family =
bernoulli("probit")) + set_rescor(FALSE)`, adding `+ RT` to the choice formula
for the speed-accuracy coupling. A family would be worth rebuilding only for
someone who specifically wants marginal-scale choice coefficients under random
effects; if so, build the lognorbit form on the shifted machinery, and guard
`u * sinh(rho)`, whose overflow gave the prototype 3 rejected warmup
proposals (a NaN in `log_mix()`).

## Files

| file | what |
| --- | --- |
| `gaussbit.patch` | The whole family as it was at 22d45d8: R file, tests, docs, the priors / inits / checkdata / stanvars hooks, NEWS, pkgdown, and the vignette section. `git apply benchmarks/gaussbit/gaussbit.patch` restores it (checked clean against the tree it was dropped from), then `roxygen2::roxygenise(".")`. Not included: the cached `vignettes/models/m_gaussbit.rds` and the regenerated `man/figures/decision_making1.png`, which were binary; refit and re-render them. |
| `decomposition.R` | The maximum-likelihood table above, by-cell differences and error-RT means. About a minute. Self-contained: the gaussbit density is written out in the script. |
| `lognorbit_cost.R` | The fit and gradient timing above, with the lognorbit prototype written in the script. About 15 minutes. `show` prints the generated program. |
| `results/` | `decomposition.csv`, `lognorbit_fits.csv`, `lognorbit_gradients.csv` (run 2). |
| `links/` | The probit-vs-logit benchmark that chose the link; `bench.R` needs the patch applied. |
