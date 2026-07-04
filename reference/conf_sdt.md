# Confidence Signal Detection Theory Model (EXPERIMENTAL)

Functions for the Signal Detection Theory (SDT) model of confidence.
This model assumes that a single sample of sensory evidence is used to
generate both a discrimination response and a confidence rating.

Functions:

- `rconf_sdt()`: Simulates random draws from the SDT model.

- `dconf_sdt()`: Computes the likelihood/density of observed responses
  and ratings.

## Usage

``` r
rconf_sdt(n, dprime, c, thetazero, thetaone, truth = NULL)

dconf_sdt(
  truth,
  response,
  confidence,
  dprime,
  c,
  thetazero,
  thetaone,
  log = FALSE
)

conf_sdt_stanvars()

conf_sdt_custom_family(
  link_mu = "identity",
  link_c = "identity",
  link_tzeroone = "log",
  link_tzerotwo = "log",
  link_toneone = "log",
  link_tonetwo = "log"
)

log_lik_conf_sdt(i, prep)

posterior_predict_conf_sdt(i, prep, ...)

posterior_epred_conf_sdt(prep)
```

## Arguments

- n:

  Number of simulated trials.

- dprime:

  Sensitivity parameter(s) (d'). Can be a single value or a vector.

- c:

  Response bias parameter. Must be a single value.

- thetazero:

  A numeric vector of confidence criteria for response=0. Must be sorted
  in decreasing order.

- thetaone:

  A numeric vector of confidence criteria for response=1. Must be sorted
  in increasing order.

- truth:

  Stimulus value (0 or 1). For simulation, can be a vector of length
  `n`. If NULL, it's sampled. For density, it's the observed stimulus.

- response:

  Observed response (0 or 1). For density calculation.

- confidence:

  Observed confidence rating (integer \>= 0). For density calculation.

- log:

  Logical; if TRUE, returns the log-density.

- link_mu:

  Link function for the dprime parameter.

- link_c:

  Link function for the c parameter.

- link_tzeroone, link_tzerotwo, link_toneone, link_tonetwo:

  Link functions for the criteria offset parameters.

- i, prep:

  For brms' functions to run: index of the observation and a `brms`
  preparation object.

- ...:

  Additional arguments.

## Details

The SDT model assumes that for a given stimulus `truth` (0 or 1),
sensory evidence `x` is drawn from a normal distribution. When `truth`
is 1, the evidence is drawn from `N(dprime / 2, 1)`; when `truth` is 0,
it's drawn from `N(-dprime / 2, 1)`. The `response` is 1 if `x > -c` and
0 otherwise. The decision criterion is thus `-c`. The confidence rating
is determined by which region the evidence `x` falls into, defined by
the confidence criteria `thetaone` and `thetazero`. A rating of 0
corresponds to the lowest confidence (i.e., evidence closest to the
decision criterion).

**Parameter Interpretation:**

- `dprime`: The sensitivity (d'). It reflects the observer's ability to
  distinguish between the two stimulus categories. Higher values
  indicate better performance.

- `c`: The response bias. It reflects the observer's tendency to favor
  one response over the other, independent of the stimulus. A positive
  value indicates a bias towards response 1, and a negative value
  indicates a bias towards response 0. A value of 0 is unbiased.

- `thetazero` & `thetaone`: The confidence criteria. These thresholds
  partition the evidence space into different confidence levels for each
  response.

  - For `response = 1`, `thetaone` are the criteria. They must be
    greater than the decision criterion (`-c`). For a given response,
    higher confidence ratings correspond to values of `x` further away
    from the decision criterion. For example, `confidence = 0` if
    `-c < x < thetaone[1]`, `confidence = 1` if
    `thetaone[1] < x < thetaone[2]`, and so on.

  - For `response = 0`, `thetazero` are the criteria. They must be less
    than the decision criterion (`-c`). For example, `confidence = 0` if
    `thetazero[1] < x < -c`, `confidence = 1` if
    `thetazero[2] < x < thetazero[1]`, and so on.

## References

- Green, D. M., & Swets, J. A. (1966). Signal detection theory and
  psychophysics. Wiley.

## Examples

``` r
# Simulate data with 2 criteria per response (3 confidence levels: 0, 1, 2)
sim_data <- rconf_sdt(
  n = 1000, dprime = 1, c = 0.2,
  thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
)
table(sim_data$response, sim_data$confidence)
#>    
#>       0   1   2
#>   0 106 236  94
#>   1 261 212  91

# Calculate density for confidence=2 (highest) for response=1
dconf_sdt(
  truth = 1, response = 1, confidence = 2, dprime = 1, c = 0.2,
  thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
)
#> [1] 0.1586553
```
