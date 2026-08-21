# Deprecated family and function names

The model families were renamed to a single `cogmod_*` scheme, so that
every family reads as one of this package's and none of them can collide
with a `brms` or `stats` name. The old names still work and are exact
synonyms, but they are no longer documented and will be removed in a
future release.

|  |  |
|----|----|
| old | new |
| `rt_lognormal()` | [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md) |
| `rt_loggamma()` | [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md) |
| `rt_invgaussian()` | [`cogmod_invgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgaussian.md) |
| `rt_gamma()` | [`cogmod_gamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_gamma.md) |
| `rt_invgamma()` | [`cogmod_invgamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invgamma.md) |
| `rt_weibull()` | [`cogmod_weibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_weibull.md) |
| `rt_invweibull()` | [`cogmod_invweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_invweibull.md) |
| `rt_logweibull()` | [`cogmod_logweibull()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_logweibull.md) |
| `rt_exgaussian()` | [`cogmod_exgaussian()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_exgaussian.md) |
| `rt_lba()` | [`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md) |
| `lba()` | [`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md) |
| `lnr()` | [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md) |
| `rdm()` | [`cogmod_rdm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_rdm.md) |
| `ddm()` | [`cogmod_ddm()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md) |
| `choco()` | [`cogmod_choco()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md) |
| `betagate()` | [`cogmod_betagate()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betagate.md) |
| `betadiscrete()` | [`cogmod_betadiscrete()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_betadiscrete.md) |

Every derived function follows its family: `rrt_lognormal()` is
[`rcogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
`rt_lognormal_stanvars()` is
[`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
and so on for the density, the `*_stanvars()`, the `*_lpdf_expose()` and
the `brms` post-processing hooks.

## Usage

``` r
drt_lognormal(x, mu = -0.7, sigma = 0.5, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_lognormal(i, prep)

posterior_epred_rt_lognormal(prep, predict_outliers = NULL)

posterior_predict_rt_lognormal(i, prep, predict_outliers = NULL, ...)

rrt_lognormal(n, mu = -0.7, sigma = 0.5, ndt = 0.2, poutlier = 0)

rt_lognormal(
  link_mu = "identity",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_lognormal_lpdf_expose()

rt_lognormal_stanvars()

drt_loggamma(
  x,
  mu = -0.7,
  sigma = 0.5,
  shape = 0,
  ndt = 0.2,
  poutlier = 0,
  log = FALSE
)

log_lik_rt_loggamma(i, prep)

posterior_epred_rt_loggamma(prep, predict_outliers = NULL)

posterior_predict_rt_loggamma(i, prep, predict_outliers = NULL, ...)

rrt_loggamma(n, mu = -0.7, sigma = 0.5, shape = 0, ndt = 0.2, poutlier = 0)

rt_loggamma(
  link_mu = "identity",
  link_sigma = "softplus",
  link_shape = "identity",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_loggamma_lpdf_expose()

rt_loggamma_stanvars()

drt_invgaussian(
  x,
  drift = 3,
  boundary = 0.5,
  ndt = 0.2,
  sigmadrift = 0,
  poutlier = 0,
  log = FALSE
)

log_lik_rt_invgaussian(i, prep)

posterior_epred_rt_invgaussian(prep, predict_outliers = NULL)

posterior_predict_rt_invgaussian(i, prep, predict_outliers = NULL, ...)

prt_invgaussian(
  q,
  drift = 3,
  boundary = 0.5,
  ndt = 0.2,
  sigmadrift = 0,
  poutlier = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rrt_invgaussian(
  n,
  drift = 3,
  boundary = 0.5,
  ndt = 0.2,
  sigmadrift = 0,
  poutlier = 0
)

rt_invgaussian(
  link_mu = "softplus",
  link_boundary = "softplus",
  link_sigmadrift = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_invgaussian_lpdf_expose()

rt_invgaussian_stanvars()

drt_gamma(x, mu = 3, sigma = 0.15, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_gamma(i, prep)

posterior_epred_rt_gamma(prep, predict_outliers = NULL)

posterior_predict_rt_gamma(i, prep, predict_outliers = NULL, ...)

rrt_gamma(n, mu = 3, sigma = 0.15, ndt = 0.2, poutlier = 0)

rt_gamma(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_gamma_lpdf_expose()

rt_gamma_stanvars()

drt_invgamma(x, mu = 4, sigma = 1.5, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_invgamma(i, prep)

posterior_epred_rt_invgamma(prep, predict_outliers = NULL)

posterior_predict_rt_invgamma(i, prep, predict_outliers = NULL, ...)

rrt_invgamma(n, mu = 4, sigma = 1.5, ndt = 0.2, poutlier = 0)

rt_invgamma(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_invgamma_lpdf_expose()

rt_invgamma_stanvars()

drt_weibull(x, mu = 2, sigma = 0.5, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_weibull(i, prep)

posterior_epred_rt_weibull(prep, predict_outliers = NULL)

posterior_predict_rt_weibull(i, prep, predict_outliers = NULL, ...)

rrt_weibull(n, mu = 2, sigma = 0.5, ndt = 0.2, poutlier = 0)

rt_weibull(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_weibull_lpdf_expose()

rt_weibull_stanvars()

drt_invweibull(x, mu = 3, sigma = 0.4, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_invweibull(i, prep)

posterior_epred_rt_invweibull(prep, predict_outliers = NULL)

posterior_predict_rt_invweibull(i, prep, predict_outliers = NULL, ...)

rrt_invweibull(n, mu = 3, sigma = 0.4, ndt = 0.2, poutlier = 0)

rt_invweibull(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_invweibull_lpdf_expose()

rt_invweibull_stanvars()

drt_logweibull(x, mu = -0.8, sigma = 0.3, ndt = 0.2, poutlier = 0, log = FALSE)

log_lik_rt_logweibull(i, prep)

posterior_epred_rt_logweibull(prep, predict_outliers = NULL)

posterior_predict_rt_logweibull(i, prep, predict_outliers = NULL, ...)

rrt_logweibull(n, mu = -0.8, sigma = 0.3, ndt = 0.2, poutlier = 0)

rt_logweibull(
  link_mu = "identity",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_logweibull_lpdf_expose()

rt_logweibull_stanvars()

drt_exgaussian(x, mu = 0.5, sigma = 0.1, tau = 0.2, log = FALSE)

log_lik_rt_exgaussian(i, prep)

posterior_epred_rt_exgaussian(prep)

posterior_predict_rt_exgaussian(i, prep, ...)

rrt_exgaussian(n, mu = 0.5, sigma = 0.1, tau = 0.2)

rt_exgaussian(
  link_mu = "identity",
  link_sigma = "softplus",
  link_tau = "softplus"
)

rt_exgaussian_lpdf_expose()

rt_exgaussian_stanvars()

drt_lba(
  x,
  drift = 3,
  sigma = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.3,
  poutlier = 0,
  log = FALSE
)

log_lik_rt_lba(i, prep)

posterior_epred_rt_lba(prep, predict_outliers = NULL)

posterior_predict_rt_lba(i, prep, predict_outliers = NULL, ...)

rrt_lba(
  n,
  drift = 3,
  sigma = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.3,
  poutlier = 0
)

rt_lba(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_sigmabias = "softplus",
  link_boundary = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rt_lba_lpdf_expose()

rt_lba_stanvars()

dlba(
  x,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.2,
  response,
  poutlier = 0,
  log = FALSE
)

lba(
  link_mu = "identity",
  link_driftone = "identity",
  link_sigmazero = "softplus",
  link_sigmaone = "softplus",
  link_sigmabias = "softplus",
  link_boundary = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

lba_lpdf_expose()

lba_stanvars()

log_lik_lba(i, prep)

posterior_epred_lba(prep)

posterior_predict_lba(i, prep, predict_outliers = NULL, ...)

rlba(
  n,
  driftzero = 3,
  driftone = 3,
  sigmazero = 1,
  sigmaone = 1,
  sigmabias = 0.5,
  boundary = 0.5,
  ndt = 0.2,
  poutlier = 0
)

dlnr(
  x,
  nuzero = 0,
  nuone = 0,
  sigmazero = 1,
  sigmaone = 1,
  ndt = 0.2,
  response,
  poutlier = 0,
  log = FALSE
)

lnr(
  link_mu = "identity",
  link_nuone = "identity",
  link_sigmazero = "softplus",
  link_sigmaone = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

lnr_lpdf_expose()

lnr_stanvars()

log_lik_lnr(i, prep)

posterior_epred_lnr(prep)

posterior_predict_lnr(i, prep, predict_outliers = NULL, ...)

rlnr(
  n,
  nuzero = 0,
  nuone = 0,
  sigmazero = 1,
  sigmaone = 1,
  ndt = 0.2,
  poutlier = 0
)

drdm(
  x,
  vzero = 3,
  vone = 2,
  boundary = 0.5,
  bias = 0.2,
  ndt = 0.2,
  response = NULL,
  poutlier = 0,
  log = FALSE
)

log_lik_rdm(i, prep)

posterior_epred_rdm(prep)

posterior_predict_rdm(i, prep, predict_outliers = NULL, ...)

prdm(
  q,
  vzero = 3,
  vone = 2,
  boundary = 0.5,
  bias = 0.2,
  ndt = 0.2,
  poutlier = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rdm(
  link_mu = "softplus",
  link_driftone = "softplus",
  link_sigmabias = "softplus",
  link_boundary = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

rdm_lpdf_expose()

rdm_stanvars()

rrdm(
  n,
  vzero = 3,
  vone = 2,
  boundary = 0.5,
  bias = 0.2,
  ndt = 0.2,
  poutlier = 0
)

dddm(
  x,
  drift = 0,
  boundary = 1,
  bias = 0.5,
  ndt = 0.2,
  response,
  sigmadrift = 0,
  sigmabias = 0,
  sigmandt = 0,
  poutlier = 0,
  log = FALSE
)

ddm(
  link_mu = "identity",
  link_boundary = "softplus",
  link_bias = "logit",
  link_sigmadrift = "softplus",
  link_sigmabias = "logit",
  link_sigmandt = "log",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
)

ddm_lpdf_expose()

ddm_stanvars()

log_lik_ddm(i, prep)

posterior_epred_ddm(prep, predict_outliers = NULL)

posterior_predict_ddm(i, prep, predict_outliers = NULL, ...)

rddm(
  n,
  drift = 0,
  boundary = 1,
  bias = 0.5,
  ndt = 0.2,
  sigmadrift = 0,
  sigmabias = 0,
  sigmandt = 0,
  poutlier = 0
)

choco(
  link_mu = "logit",
  link_confright = "logit",
  link_precright = "softplus",
  link_confleft = "logit",
  link_precleft = "softplus",
  link_pex = "logit",
  link_bex = "logit",
  link_pmid = "logit"
)

choco_lpdf_expose()

choco_stanvars()

dchoco(
  x,
  p = 0.5,
  confright = 0.5,
  precright = 4,
  confleft = 0.5,
  precleft = 4,
  pex = 0.1,
  bex = 0.5,
  pmid = 0,
  mid = 0.5,
  log = FALSE
)

log_lik_choco(i, prep)

posterior_epred_choco(prep)

posterior_predict_choco(i, prep, ...)

rchoco(
  n,
  p = 0.5,
  confright = 0.5,
  precright = 4,
  confleft = 0.5,
  precleft = 4,
  pex = 0.1,
  bex = 0.5,
  pmid = 0,
  mid = 0.5
)

betagate(
  link_mu = "logit",
  link_phi = "softplus",
  link_pex = "logit",
  link_bex = "logit"
)

betagate_lpdf_expose()

betagate_stanvars()

dbetagate(x, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5, log = FALSE)

log_lik_betagate(i, prep)

posterior_epred_betagate(prep)

posterior_predict_betagate(i, prep, ...)

rbetagate(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5)

betadiscrete(link_mu = "logit", link_phi = "log", link_pzero = "logit")

betadiscrete_lpmf_expose()

betadiscrete_stanvars()

dbetadiscrete(x, mu = 0.5, phi = 3, k = 5, pzero = 0, log = FALSE)

log_lik_betadiscrete(i, prep)

pbetadiscrete(
  q,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

posterior_epred_betadiscrete(prep)

posterior_predict_betadiscrete(i, prep, ...)

qbetadiscrete(
  p,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
)

rbetadiscrete(n, mu = 0.5, phi = 3, k = 5, pzero = 0)
```
