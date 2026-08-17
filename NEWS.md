# cogmod 0.1.0

First CRAN release.

## Models for subjective scales

* Beta-Gate (`betagate()`, `rbetagate()`, `dbetagate()`), a reparametrised
  ordered beta model.
* Discrete Beta (`betadiscrete()` and friends) for Likert-type responses.
* Choice-Confidence, CHOCO (`choco()`, `rchoco()`, `dchoco()`) for bipolar
  scales.
* Signal detection with confidence ratings (`conf_sdt()`).

## Models for decision making

* Lognormal race (`lnr()`), linear ballistic accumulator (`lba()`), drift
  diffusion with optional across-trial variability (`ddm()`), and the racing
  diffusion model (`rdm()`).

## Models for reaction times alone

* A consistently parametrised set of shifted, right-skewed response
  distributions: `rt_lognormal()`, `rt_invgaussian()`, `rt_gamma()`,
  `rt_invgamma()`, `rt_weibull()`, `rt_logweibull()`, `rt_invweibull()`,
  `rt_exgaussian()` and `rt_lba()`.

## Data

* `wagenmakers2008`, lexical decision data from Wagenmakers et al. (2008).
* `badlm`, a simulated dataset in which two conditions share a mean reaction
  time while differing in shift, spread and tail weight.
