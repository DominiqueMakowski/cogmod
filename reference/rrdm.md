# Simulate from the Two-Accumulator Racing Diffusion Model (RDM)

Simulates choice reaction times from a two-accumulator Racing Diffusion
Model (RDM). This is a specialized version where exactly two
accumulators race towards a common threshold. The model assumes
variability in the starting point of the diffusion process, drawn from a
uniform distribution. This version is optimized for performance using
vectorized operations, and allows drift rates to be zero (a zero-drift
accumulator is slow but still finishes; see Details).

## Usage

``` r
rrdm(n, vzero, vone, bs, bias, ndt)

drdm(x, response = NULL, vzero, vone, bs, bias, ndt, log = FALSE)

prdm(q, vzero, vone, bs, bias, ndt, lower.tail = TRUE, log.p = FALSE)

rdm_lpdf_expose()

rdm_stanvars()

rdm(
  link_mu = "softplus",
  link_driftone = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

log_lik_rdm(i, prep, ...)

posterior_predict_rdm(i, prep, ...)

posterior_epred_rdm(prep)
```

## Arguments

- n:

  Number of trials to simulate. Must be a positive integer.

- vzero:

  Drift rate for the first accumulator. Must be a single non-negative
  number.

- vone:

  Drift rate for the second accumulator. Must be a single non-negative
  number.

- bs:

  Threshold parameter, defined as `bs = b - bias`, where `b` is the
  decision threshold and `bias` is the maximum starting point. Must be a
  single positive number.

- bias:

  Maximum starting point parameter. The starting point for each
  accumulator on each trial is drawn from `Uniform(0, bias)`. Must be a
  single positive number.

- ndt:

  Non-decision time (encoding and motor time offset). Must be a single
  non-negative number.

- x:

  Vector of quantiles (observed reaction times).

- response:

  Accumulator whose finishing time is being scored: `0` for the `vzero`
  accumulator, `1` for the `vone` accumulator. This gives the
  *defective* density `f_response(x) * S_other(x)`, which is what a race
  likelihood needs. The default `NULL` instead returns the *marginal*
  density of the winning RT, ignoring which accumulator won.

- log:

  Logical; if TRUE, probabilities p are given as log(p).

- q:

  Vector of quantiles (reaction times).

- lower.tail:

  If `TRUE` (default) return `P(RT <= q)`, otherwise the survival
  `P(RT > q)`.

- log.p:

  If `TRUE`, probabilities are returned on the log scale.

- link_mu, link_driftone, link_sigmabias, link_bs, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Value

A data frame with `n` rows and two columns:

- rt:

  The simulated reaction time (minimum finishing time across the two
  accumulators).

- response:

  The winning accumulator, coded `0` for `vzero` and `1` for `vone`,
  matching the `dec()` coding used by the `brms` families.

## Details

The RDM implemented here follows the formulation where the two
accumulators have drift rates `vzero` and `vone`. The diffusion process
starts at a point `z` drawn from `Uniform(0, bias)`. The process
terminates when either accumulator reaches a threshold `b`. The
parameter `bs` is defined as `bs = b - bias`, representing the distance
from the maximum starting point `bias` to the threshold `b`. Therefore,
the effective distance to threshold for a given trial is
`bs = b - z = bs + bias - z`.

The finishing time for a single accumulator, given its drift rate `v`,
`bs`, `bias`, and `ndt`, is simulated by:

1.  Sampling a starting point `z ~ Uniform(0, bias)`.

2.  Calculating the distance `bs = bs + bias - z`.

3.  If `v > 0`, simulating the time to reach `bs` from an Inverse
    Gaussian distribution with mean `bs / v` and shape `bs^2`. This
    simulation uses an internal implementation based on Michael et al.
    (1976).

4.  If `v = 0`, drawing the driftless first passage time, which is Levy
    distributed: `bs^2 / Z^2` with `Z ~ Normal(0, 1)`. A zero-drift
    accumulator is slow, but it still finishes with probability one, so
    it can win the race.

5.  Adding the non-decision time `ndt` to the finishing times.

The function simulates this process for both accumulators using
vectorized operations. The accumulator that finishes first determines
the response (0 for `vzero`, 1 for `vone`) and the reaction time (RT)
for that trial.

This implementation is based on the description and parameters used in:
Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling
models without random between-trial variability: The racing diffusion
model of speeded decision making. *Psychonomic Bulletin & Review*, *27*,
911-936.
[doi:10.3758/s13423-020-01738-8](https://doi.org/10.3758/s13423-020-01738-8)
(specifically matching the `WaldA` component used within their RDM
simulation).

`prdm()` describes the RT of the race as a whole (whichever accumulator
wins), since `P(min(T0, T1) > q) = S0(q) * S1(q)`. There is no
comparably simple closed form for the per-response defective CDF, so
`prdm()` takes no `response` argument.

The `brms` family names the drift of the first accumulator `mu` (as
`brms` requires) and that of the second `driftone`, and calls the
start-point range `sigmabias` to match
[`lba()`](https://github.com/DominiqueMakowski/cogmod/reference/rlba.md),
where it denotes the same quantity. Note that this is *not* the same
thing as `bias` in
[`ddm()`](https://github.com/DominiqueMakowski/cogmod/reference/rddm.md),
which is a relative starting point in `[0, 1]`. Both drifts use a
softplus link with a lower bound of zero, following
[`rt_invgaussian()`](https://github.com/DominiqueMakowski/cogmod/reference/rrt_invgaussian.md):
a Wald drift must be non-negative for the accumulator to be a proper
first passage time.

Note that `sigmabias` and `bs` are only weakly identified from each
other, because they enter the threshold as the sum `b = bs + sigmabias`
and trade off almost freely: on simulated data with 4000 trials the
profile log-likelihood varies by only about 3 units as `sigmabias`
ranges from 0 to half the threshold, while `bs` slides to compensate.
With flat priors the sampler tends to wander down the `sigmabias -> 0`
ridge (the plain Wald race) and produce divergent transitions. A weakly
informative prior on `sigmabias` fixes this, for example
`brms::prior(normal(-1, 1), class = "Intercept", dpar = "sigmabias")`,
which under the softplus link is centred near `0.31`. The sum
`bs + sigmabias` is well identified either way, so it is the more
trustworthy quantity to interpret and to compare across conditions. The
same caveat applies to
[`lba()`](https://github.com/DominiqueMakowski/cogmod/reference/rlba.md),
which shares this parameterisation.

## References

- Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating
  Random Variates Using Transformations with Multiple Roots. *The
  American Statistician*, *30*(2), 88–90.
  [doi:10.2307/2683801](https://doi.org/10.2307/2683801)

- Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling
  models without random between-trial variability: The racing diffusion
  model of speeded decision making. *Psychonomic Bulletin & Review*,
  *27*, 911-936.
  [doi:10.3758/s13423-020-01738-8](https://doi.org/10.3758/s13423-020-01738-8)

- Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian
  distribution and its statistical application—a review. *Journal of the
  Royal Statistical Society Series B: Statistical Methodology*, *40*(3),
  263-275.

## See also

`rrt_invgaussian`

## Examples

``` r
rdm_pos <- rrdm(n = 1000, vzero = 0.8, vone = 0.6, bs = 0.5, bias = 0.2, ndt = 0.15)

# You can expose the lpdf function as follows:
# rdm_lpdf <- rdm_lpdf_expose()
# rdm_lpdf(0.5, 2, 1.5, 0.2, 0.5, 0.5, 0.2, 0)
```
