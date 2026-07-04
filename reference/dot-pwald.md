# Analytical CDF for Wald with Uniform Start Point Variability

Calculates the cumulative distribution function (CDF) for the time it
takes a diffusion process with drift `drift`, starting point
`z ~ U(0, bias)`, to reach threshold `b = bs + bias`, shifted by `ndt`.
Based on Tillman et al. (2020). This function calculates the CDF for the
*unadjusted* time `x`.

## Usage

``` r
.pwald(x, drift, bs, bias, ndt, lower.tail = TRUE, log.p = FALSE)
```
