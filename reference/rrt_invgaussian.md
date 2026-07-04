# Inverse Gaussian Model (Shifted Wald)

Density, distribution function, and random generation for the Shifted
Wald distribution, also known as the Shifted Inverse Gaussian
distribution. This distribution is commonly used in modeling reaction
times in cognitive tasks. It is characterized by a drift rate (`drift`),
a decision threshold (`bs`), and a non-decision time (`ndt`).

Functions:

- `rrt_invgaussian()`: Simulates random draws from the Shifted Wald
  model.

- `drt_invgaussian()`: Computes the density (likelihood) of the Shifted
  Wald distribution.

- `prt_invgaussian()`: Computes the cumulative distribution function
  (CDF).

## Usage

``` r
rrt_invgaussian(n, drift = 3, bs = 0.5, ndt = 0.2)

drt_invgaussian(x, drift = 3, bs = 0.5, ndt = 0.2, log = FALSE)

prt_invgaussian(
  q,
  drift = 3,
  bs = 0.5,
  ndt = 0.2,
  lower.tail = TRUE,
  log.p = FALSE
)

rt_invgaussian_lpdf_expose()

rt_invgaussian_stanvars()

rt_invgaussian(
  link_mu = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

log_lik_rt_invgaussian(i, prep, ...)

posterior_predict_rt_invgaussian(i, prep, ...)

posterior_epred_rt_invgaussian(prep)
```

## Arguments

- n:

  Number of observations. If `length(n) > 1`, the length is taken to be
  the number required.

- drift:

  Drift rate. Must be positive. Represents the average speed of evidence
  accumulation. Range: (0, Inf).

- bs:

  Decision threshold (boundary separation). Must be positive. Represents
  the amount of evidence needed to make a decision. Range: (0, Inf).

- ndt:

  Non-decision time (shift parameter). Must be non-negative. Represents
  time for processes like stimulus encoding and response execution.
  Range: \[0, Inf).

- x:

  Vector of quantiles (observed reaction times).

- log:

  Logical; if TRUE, probabilities p are given as log(p).

- q:

  Vector of quantiles (observed reaction times).

- lower.tail:

  Logical; if TRUE (default), probabilities are `P[X <= x]`, otherwise,
  `P[X > x]`.

- log.p:

  Logical; if TRUE, probabilities p are given as log(p). Defaults to
  FALSE.

- link_mu, link_bs, link_tau, link_minrt:

  Link functions for the parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Details

The Shifted Wald distribution describes the time it takes for a Wiener
diffusion process starting at 0 to reach a threshold `bs` \> 0, given a
positive drift rate `drift` \> 0. The resulting time is then shifted by
a non-decision time `ndt` \>= 0.

The distribution is mathematically equivalent to shifting an Inverse
Gaussian distribution with mean `mu = bs / drift` and shape parameter
`lambda = bs^2`. That is,
`ShiftedWald(drift, bs, ndt) = InverseGaussian(mean = bs/drift, shape = bs^2) + ndt`.

The random generation algorithm implemented here is based on the method
described by Michael, Schucany, and Haas (1976), as used in the
`statmod` package.

## References

- Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating
  Random Variates Using Transformations with Multiple Roots. *The
  American Statistician*, *30*(2), 88–90.
  [doi:10.2307/2683801](https://doi.org/10.2307/2683801)

- Anders, R., Alario, F., & Van Maanen, L. (2016). The shifted Wald
  distribution for response time data analysis. *Psychological Methods*,
  *21*(3), 309–327.
  [doi:10.1037/met0000063](https://doi.org/10.1037/met0000063)

- Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation
  of the ex-Gaussian and shifted Wald parameters: A diffusion model
  analysis. *Psychonomic Bulletin & Review*, *16*(5), 798–817.
  [doi:10.3758/PBR.16.5.798](https://doi.org/10.3758/PBR.16.5.798)

- Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian
  distribution and its statistical application—a review. *Journal of the
  Royal Statistical Society Series B: Statistical Methodology*, *40*(3),
  263-275.

## Examples

``` r
# Simulate 1000 RTs
rts <- rrt_invgaussian(1000, drift = 3, bs = 0.5, ndt = 0.2)
# hist(rts, breaks = 50, main = "Simulated Shifted Wald RTs", xlab = "Reaction Time")
```
