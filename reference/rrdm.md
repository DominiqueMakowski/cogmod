# Simulate from the Two-Accumulator Racing Diffusion Model (RDM)

Simulates choice reaction times from a two-accumulator Racing Diffusion
Model (RDM). This is a specialized version where exactly two
accumulators race towards a common threshold. The model assumes
variability in the starting point of the diffusion process, drawn from a
uniform distribution. This version is optimized for performance using
vectorized operations and allows one (but not both) drift rate to be
zero.

## Usage

``` r
rrdm(n, vzero, vone, bs, bias, ndt)

drdm(x, vzero, vone, bs, bias, ndt, log = FALSE)
```

## Arguments

- n:

  Number of trials to simulate. Must be a positive integer.

- vzero:

  Drift rate for the first accumulator. Must be a single non-negative
  number.

- vone:

  Drift rate for the second accumulator. Must be a single non-negative
  number. At least one of `vzero` or `vone` must be positive.

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

- log:

  Logical; if TRUE, probabilities p are given as log(p).

## Value

bias data frame with `n` rows and two columns:

- rt:

  The simulated reaction time (minimum finishing time across the two
  accumulators).

- choice:

  The index of the winning accumulator (1 for `vzero`, 2 for `vone`).

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

4.  If `v = 0`, the finishing time is considered infinite (`Inf`).

5.  Adding the non-decision time `ndt` to finite finishing times.

The function simulates this process for both accumulators using
vectorized operations. The accumulator that finishes first determines
the choice (1 for `vzero`, 2 for `vone`) and the reaction time (RT) for
that trial. If one drift rate is zero, the accumulator with the positive
drift rate will always win.

This implementation is based on the description and parameters used in:
Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling
models without random between-trial variability: The racing diffusion
model of speeded decision making. *Psychonomic Bulletin & Review*, *27*,
911-936.
[doi:10.3758/s13423-020-01738-8](https://doi.org/10.3758/s13423-020-01738-8)
(specifically matching the `WaldA` component used within their RDM
simulation).

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
```
