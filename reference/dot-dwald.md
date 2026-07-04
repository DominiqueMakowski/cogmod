# Analytical PDF for Wald with Uniform Start Point Variability

Calculates the probability density function (PDF) for the time it takes
a diffusion process with drift `drift`, starting point `z ~ U(0, bias)`,
to reach threshold `b = bs + bias`, shifted by `ndt`. Based on Tillman
et al. (2020). This function calculates the density for the *unadjusted*
time `x`.

## Usage

``` r
.dwald(x, drift, bs, bias, ndt, log = FALSE)
```
