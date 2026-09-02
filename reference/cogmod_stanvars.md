# The Stan code a cogmod family needs, read off the model

Returns the `stanvars` argument for
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html), for
whichever `cogmod` family the model uses. It is a front end to the
per-family `<family>_stanvars()` functions -
[`cogmod_lognormal_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
[`cogmod_choco_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_choco.md),
[`cogmod_ddm_stanvars()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_ddm.md)
and the rest - which remain available and unchanged.

## Usage

``` r
cogmod_stanvars(formula, ...)
```

## Arguments

- formula:

  A
  [`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
  formula carrying the family, a `cogmod` family object, or a fitted
  `brmsfit`.

- ...:

  Passed to the family's own `<family>_stanvars()` function.

## Value

A `stanvars` object, to pass to `brms::brm(stanvars = )`.

## Details

`brms` needs the custom likelihood injected into the generated Stan
program, and every `cogmod` family ships one. Calling this instead of
the family's own function means the family is named once, in
[`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html),
rather than twice:

    f <- brms::bf(RT ~ Condition, ndt ~ Condition,
                  family = cogmod_lognormal())

    brms::brm(f, data = df,
              prior    = cogmod_priors(f, df),
              init     = cogmod_inits(f, df),
              stanvars = cogmod_stanvars(f))

## A warning it may emit

[`cogmod_lba1()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba1.md)
and
[`cogmod_lba2()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lba2.md)
have a likelihood that is *exactly* constant along the ray that
multiplies the drift rates, their SDs, the start-point range and the
threshold offset by a common factor. If the formula pins none of them to
a constant, this warns: the RT distribution is still identified, but the
individual parameters are not, and the fit will converge to whatever the
priors say about that direction rather than fail. Fixing any one member
in [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html) -
conventionally `sigmazero = 1` - silences it. Leaving a parameter *out*
of [`bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
does not count: `brms` estimates it anyway.

## What it accepts

A
[`brms::bf()`](https://paulbuerkner.com/brms/reference/brmsformula.html)
formula carrying the family, the family object itself, or a fitted
`brmsfit` (useful for recompiling or for
[`update()`](https://rdrr.io/r/stats/update.html)).

## See also

[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md),
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)

## Examples

``` r
f <- brms::bf(RT ~ 1, ndt ~ 1, family = cogmod_lognormal())
cogmod_stanvars(f)
#> [[1]]
#> [[1]]$name
#> [1] ""
#> 
#> [[1]]$sdata
#> NULL
#> 
#> [[1]]$scode
#> [1] "\n// Log-likelihood for one observation from the shifted LogNormal model.\n// Y: observed reaction time.\n// mu: mean of the decision time on the log scale (meanlog).\n// sigma: SD of the decision time on the log scale (> 0).\n// ndt: non-decision time, same unit as Y (> 0).\n// poutlier: proportion of responses from the outlier process, in [0, 1].\n//\n// The outlier component is a half Normal with scale 0.2 s. It keeps the density\n// strictly positive below `ndt`, where the shifted decision component has none.\n// That is what removes the hard min-RT boundary and lets `ndt` be estimated\n// directly rather than as a fraction of an observed minimum. The scale is a\n// constant in SECONDS. This family expects reaction times in seconds; give it\n// another unit and the component contributes nothing anywhere in the data,\n// which silently reinstates the min-RT boundary it exists to remove.\n//\n// It is written out rather than called as normal_lpdf() because both of its\n// parameters are constant: written that way, Stan recomputes the normalising\n// constant for every observation on every leapfrog step.\nreal cogmod_lognormal_lpdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    // Parameter checks\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n\n    // The leading constant includes the log(2) that folds the symmetric\n    // Normal onto [0, Inf).\n    real lp_out = 1.3836465597893728 - 12.5 * square(Y);\n    real t_adj  = Y - ndt;\n\n    // Faster than the non-decision time: only the outlier component can have\n    // produced this response.\n    if (t_adj <= 0) return log(poutlier) + lp_out;\n\n    return log_mix(poutlier, lp_out, lognormal_lpdf(t_adj | mu, sigma));\n}\n\n// Log CDF and log survival of the same mixture, for brms's cens() addition\n// term. See ?rcogmod_invgaussian for what censoring a reaction time means and\n// when it is the right model.\nreal cogmod_lognormal_lcdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n    // Outlier CDF: 2 Phi(Y / s) - 1 = erf(Y / (s sqrt(2)))\n    real lF_out = log(erf(Y * 3.5355339059327369));\n    real t_adj  = Y - ndt;\n    if (t_adj <= 0) return log(poutlier) + lF_out;\n    return log_mix(poutlier, lF_out, lognormal_lcdf(t_adj | mu, sigma));\n}\n\nreal cogmod_lognormal_lccdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return 0;\n    // Outlier survival: 2 Phi(-Y / s), through the lower tail (see above)\n    real lS_out = 0.69314718055994529 + std_normal_lcdf(-Y * 5);\n    real t_adj  = Y - ndt;\n    // Not yet past the non-decision time: the decision process cannot have\n    // finished, so its survival is exactly 1.\n    if (t_adj <= 0) return log_mix(poutlier, lS_out, 0);\n    return log_mix(poutlier, lS_out, lognormal_lccdf(t_adj | mu, sigma));\n}\n"
#> 
#> [[1]]$block
#> [1] "functions"
#> 
#> [[1]]$position
#> [1] "start"
#> 
#> [[1]]$pll_args
#> character(0)
#> 
#> 
#> attr(,"class")
#> [1] "stanvars"

# Equivalent to naming the family a second time:
cogmod_lognormal_stanvars()
#> [[1]]
#> [[1]]$name
#> [1] ""
#> 
#> [[1]]$sdata
#> NULL
#> 
#> [[1]]$scode
#> [1] "\n// Log-likelihood for one observation from the shifted LogNormal model.\n// Y: observed reaction time.\n// mu: mean of the decision time on the log scale (meanlog).\n// sigma: SD of the decision time on the log scale (> 0).\n// ndt: non-decision time, same unit as Y (> 0).\n// poutlier: proportion of responses from the outlier process, in [0, 1].\n//\n// The outlier component is a half Normal with scale 0.2 s. It keeps the density\n// strictly positive below `ndt`, where the shifted decision component has none.\n// That is what removes the hard min-RT boundary and lets `ndt` be estimated\n// directly rather than as a fraction of an observed minimum. The scale is a\n// constant in SECONDS. This family expects reaction times in seconds; give it\n// another unit and the component contributes nothing anywhere in the data,\n// which silently reinstates the min-RT boundary it exists to remove.\n//\n// It is written out rather than called as normal_lpdf() because both of its\n// parameters are constant: written that way, Stan recomputes the normalising\n// constant for every observation on every leapfrog step.\nreal cogmod_lognormal_lpdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    // Parameter checks\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n\n    // The leading constant includes the log(2) that folds the symmetric\n    // Normal onto [0, Inf).\n    real lp_out = 1.3836465597893728 - 12.5 * square(Y);\n    real t_adj  = Y - ndt;\n\n    // Faster than the non-decision time: only the outlier component can have\n    // produced this response.\n    if (t_adj <= 0) return log(poutlier) + lp_out;\n\n    return log_mix(poutlier, lp_out, lognormal_lpdf(t_adj | mu, sigma));\n}\n\n// Log CDF and log survival of the same mixture, for brms's cens() addition\n// term. See ?rcogmod_invgaussian for what censoring a reaction time means and\n// when it is the right model.\nreal cogmod_lognormal_lcdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n    // Outlier CDF: 2 Phi(Y / s) - 1 = erf(Y / (s sqrt(2)))\n    real lF_out = log(erf(Y * 3.5355339059327369));\n    real t_adj  = Y - ndt;\n    if (t_adj <= 0) return log(poutlier) + lF_out;\n    return log_mix(poutlier, lF_out, lognormal_lcdf(t_adj | mu, sigma));\n}\n\nreal cogmod_lognormal_lccdf(real Y, real mu, real sigma, real ndt, real poutlier) {\n    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return 0;\n    // Outlier survival: 2 Phi(-Y / s), through the lower tail (see above)\n    real lS_out = 0.69314718055994529 + std_normal_lcdf(-Y * 5);\n    real t_adj  = Y - ndt;\n    // Not yet past the non-decision time: the decision process cannot have\n    // finished, so its survival is exactly 1.\n    if (t_adj <= 0) return log_mix(poutlier, lS_out, 0);\n    return log_mix(poutlier, lS_out, lognormal_lccdf(t_adj | mu, sigma));\n}\n"
#> 
#> [[1]]$block
#> [1] "functions"
#> 
#> [[1]]$position
#> [1] "start"
#> 
#> [[1]]$pll_args
#> character(0)
#> 
#> 
#> attr(,"class")
#> [1] "stanvars"
```
