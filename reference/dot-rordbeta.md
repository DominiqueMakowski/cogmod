# Generate Ordered Beta Variates

This function will generate ordered beta random variates given values
for the mean (`mu`), dispersion (`phi`) and cutpoints governing the
ratio of degenerate (discrete) to continuous responses.

## Usage

``` r
.rordbeta(n = 100, mu = 0.5, phi = 1, cutzero = -1, cutone = 1)
```

## Arguments

- n:

  Number of variates to generate.

- mu:

  Value of the mean of the distribution. Should be in the \\0,1\\
  interval (cannot be strictly equal to 0 or 1). If length is greater
  than 1, should be of length `n`.

- phi:

  Value of the dispersion parameter. Should be strictly greater than 0.
  If length is greater than 1, should be of length `n`.

- cutzero, cutone:

  Two numeric values for the cutpoints on the logit scale. Second value
  should be strictly greater than the first value.

## Value

A vector of length `n` of variates from the ordered beta distribution.
