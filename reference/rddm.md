# Drift Diffusion Model (DDM)

The Drift Diffusion Model (DDM) is a widely used model for
decision-making tasks. It assumes that evidence accumulates over time
until it reaches one of two decision boundaries.

The 'ddm' family in this package is not a full re-implementation of the
DDM, but uses of the `brms` wiener family. It is simply a reparametrized
version for consistency with the other race models in the package.

## Usage

``` r
rddm(
  n,
  drift,
  bs,
  bias,
  ndt,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL,
  ...
)

dddm(
  x,
  drift,
  bs,
  bias,
  ndt,
  response,
  log = FALSE,
  sigmadrift = 0,
  sigmabias = 0,
  sigmatau = 0,
  minrt = NULL,
  ...
)

ddm_lpdf_expose()

ddm_stanvars()

ddm(
  link_mu = "identity",
  link_bs = "softplus",
  link_bias = "logit",
  link_tau = "logit",
  link_minrt = "identity",
  link_sigmadrift = "softplus",
  link_sigmabias = "logit",
  link_sigmatau = "logit"
)

log_lik_ddm(i, prep, ...)

posterior_epred_ddm(prep)

posterior_predict_ddm(i, prep, ...)
```

## Arguments

- n:

  Number of simulated trials. Must be a positive integer.

- drift:

  Drift rate. Can take any real value.

- bs:

  Decision threshold (boundary separation). Must be positive.

- bias:

  Starting point bias, as a proportion of the boundary separation
  *measured from the boundary coded `0`*. Must be in (0, 1). See
  Details.

- ndt:

  Non-decision time. Must be non-negative.

- sigmadrift:

  Inter-trial variability in drift rate (`sv` in the usual notation).
  Must be non-negative. Default `0` (no variability, i.e. the classic
  4-parameter DDM).

- sigmabias:

  Inter-trial variability in the starting point, expressed as a fraction
  (in `[0, 1)`) of the maximum allowed range, i.e.
  `sw = sigmabias * min(2*bias, 2*(1-bias))` (`sw`/`sz` in the usual
  notation). Default `0`.

- sigmatau:

  Inter-trial variability in non-decision time, expressed as a fraction
  of `minrt`, i.e. `st0 = sigmatau * minrt` (`st0` in the usual
  notation), with `ndt` the lower bound of the resulting Uniform
  distribution. Default `0`.

- minrt:

  Minimum reaction time. Only required when `sigmatau > 0` (used to
  scale `sigmatau` into `st0`).

- ...:

  Other arguments to be passed to
  [`brms::rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
  or
  [`brms::dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html),

- x:

  The observed reaction time (RT). Must be greater than `ndt`.

- response:

  The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.

- log:

  Logical; if TRUE, returns the log-density. Default: FALSE.

- link_mu, link_bs, link_bias, link_tau, link_minrt:

  Link functions for the parameters.

- link_sigmadrift, link_sigmabias, link_sigmatau:

  Link functions for the extra inter-trial variability parameters (fix
  them to 0, e.g. `sigmadrift = 0`, in the
  [`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  formula to recover the classic 4-parameter DDM).

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

## Response coding

The response coded `1` (`response` here, `dec()` in a `brms` formula) is
the **upper** boundary and the response coded `0` is the **lower** one,
following `brms`'s own
[`wiener()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family. Two consequences are worth keeping in mind when reading a fitted
model:

- `bias` is measured *from the lower boundary*, so `bias > 0.5` places
  the starting point closer to the response coded `1`, and `bias < 0.5`
  closer to the response coded `0`. Since first-passage times are
  shorter for the nearer boundary, `bias > 0.5` makes responses coded
  `1` *faster* than responses coded `0`, and `bias < 0.5` makes them
  slower. At exactly `bias = 0.5` the two conditional RT distributions
  are identical, whatever the drift rate.

- A positive `drift` pushes the accumulator towards the response coded
  `1`. When `0` codes correct responses and `1` codes errors (a common
  choice when the task has no natural stimulus-to-boundary mapping),
  good performance therefore corresponds to a *negative* drift rate.

## Implementation

The underlying 4-parameter process is simulated and evaluated with
[`brms::rwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)/[`brms::dwiener()`](https://paulbuerkner.com/brms/reference/Wiener.html)
(which require the `RWiener` package). The full 7-parameter model is
built on top of these rather than delegated to another package:
between-trial variability is simulated by drawing the per-trial
parameters, and evaluated by combining a closed-form drift correction
with Gauss-Legendre quadrature over the starting point and non-decision
time.

## Examples

``` r
# Simulate data
# data <- rddm(1000, drift = 0.2, bs = 1, bias = 0.5, ndt = 0.3)
# hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
# hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))

# Compute density
# dddm(x = c(0.5, 0.7), drift = 0.2, bs = 1, bias = 0.5, resp = c(0, 1), ndt = 0.3)
```
