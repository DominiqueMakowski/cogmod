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
#> [1] "\n// log(Phi(x)), the one piece of arithmetic every normal tail in this package is\n// built from. Neither of Stan's two routes to it is good enough for both jobs\n// it has here, which is to be right in the far tail *and* to hand back a\n// usable derivative there. All three claims below were measured against\n// central differences of the log probability over 20000 responses.\n//\n// The erfc route - std_normal_lcdf() is not it, but lognormal_lcdf(),\n// lognormal_lccdf() and the log(u1) + log1m(u2 / u1) the LogNormal used to\n// write are - has good partials, to about 4e-6 on a summed gradient of order\n// 1e3. But erfc underflows near x = -38, and then the value is log(0) and the\n// partials are inf or 0/0. That is not a harmless -inf in a mixture:\n// log_mix() in the lpdf stays finite when the decision component is -inf, but\n// reverse mode multiplies the (zero) adjoint into the stored partial, and\n// 0 * inf is NaN, so a single response turns the gradient of the whole model\n// to NaN - 'Gradient evaluated at the initial value is not finite' at the\n// start of a fit, divergent transitions afterwards.\n//\n// std_normal_lcdf() has the range: its value is exact against R's\n// pnorm(log.p = TRUE) as far as x = -1e7. Its partials are not - they sat\n// 1.7e-3 from central differences where the erfc route sat 4e-6 on the LNR,\n// and 2e-4 to 7e-2 on the RDM, which took every tail through it - so it is not\n// a drop-in for the tails of a race, where those partials are the gradient of\n// the drifts and the scales.\n//\n// So: erfc in the body of the distribution, and below x = -25 the asymptotic\n// expansion of the tail,\n//\n//   Phi(x) = phi(x) / (-x) * (1 - 1/x^2 + 3/x^4 - 15/x^6 + 105/x^8 - 945/x^10)\n//\n// whose leading term is the exponent itself. Nothing underflows, the result\n// stays finite and differentiable as far as x = -1e150, and the six terms\n// agree with pnorm(log.p = TRUE) to 4e-16 relative from x = -25 down - the\n// last bit of a double - so the two branches meet with no step in the density.\nreal cogmod_log_Phi(real x) {\n  if (x < -25) {\n    real z = inv_square(x);\n    real series = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 - 945 * z))));\n    return -0.5 * square(x) - log(-x) - 0.91893853320467274 + log(series);\n  }\n  if (x > 0) return log1p(-0.5 * erfc(x * 0.7071067811865476));\n  return log(0.5 * erfc(-x * 0.7071067811865476));\n}\n\n// log(Phi(y + c) - Phi(y)) for c > 0, from whichever tail keeps the two terms\n// from cancelling: the upper one when y > 0, where both CDFs sit near 1, the\n// lower one otherwise. Taken as the larger tail plus log(1 - ratio) in log\n// space, the way .lognormal_ldiff_pnorm() does it - the quotient u2 / u1 this\n// replaces divides two minute numbers and lost its own accuracy long before\n// either underflowed (8e-4 against central differences, against 1e-6 here).\nreal cogmod_lognormal_ldiff_Phi(real y, real c) {\n  real hi;\n  real lo;\n  if (y > 0) {\n    hi = cogmod_log_Phi(-y);\n    lo = cogmod_log_Phi(-(y + c));\n  } else {\n    hi = cogmod_log_Phi(y + c);\n    lo = cogmod_log_Phi(y);\n  }\n  return lo < hi ? hi + log1m_exp(lo - hi) : negative_infinity();\n}\n\n// Log density of the accumulator's finishing time with start-point range A.\n// At A = 0 this is the LogNormal itself, at the LogNormal's cost.\nreal cogmod_lognormal_acc_ldens(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return lognormal_lpdf(t | meanlog, sigma);\n  real a = (meanlog - log(t)) / sigma;\n  real c = log1p(A) / sigma;\n  real x = a - sigma;\n  if (c < 1e-4) {\n    real series = 1 - c * x / 2 + square(c) * (square(x) - 1) / 6;\n    if (series <= 0) return negative_infinity();\n    return lognormal_lpdf(t | meanlog, sigma) + log(log1p(A) / A) + log(series);\n  }\n  return -meanlog + square(sigma) / 2 + cogmod_lognormal_ldiff_Phi(x, c) - log(A);\n}\n\n// [log F, log S] of the accumulator's finishing time, each computed directly\n// on the side where it is the small one. See .lognormal_acc_ltails().\nvector cogmod_lognormal_acc_ltails(real t, real meanlog, real sigma, real A) {\n  real a = (meanlog - log(t)) / sigma;\n  // At A = 0 this is the plain LogNormal, whose two tails are Phi(-a) and\n  // Phi(a). Written that way rather than as lognormal_lcdf()/lognormal_lccdf(),\n  // which are erfc alone and so reach log(0) with non-finite partials around\n  // |a| = 38 - see cogmod_log_Phi() above for what that costs.\n  if (A == 0) return [cogmod_log_Phi(-a), cogmod_log_Phi(a)]';\n  real c = log1p(A) / sigma;\n  if (c < 1e-4) {\n    real r = log1p(A) / A;\n    real corr = c * r * (0.5 + c * (2 * sigma - a) / 6);\n    real lPa = cogmod_log_Phi(a);\n    real lQa = cogmod_log_Phi(-a);\n    real lphi = std_normal_lpdf(a);\n    real lS = fmin(lPa + log1p(corr * exp(lphi - lPa)), 0);\n    real dF = 1 - corr * exp(lphi - lQa);\n    real lF = dF > 0 ? fmin(lQa + log(dF), 0) : negative_infinity();\n    return [lF, lS]';\n  }\n  real lA = log(A);\n  real lD1 = cogmod_lognormal_ldiff_Phi(a, c) - lA;\n  real lD2 = -meanlog + square(sigma) / 2 + log(t)\n             + cogmod_lognormal_ldiff_Phi(a - sigma, c) - lA;\n  if (a < 0) {\n    real lP = cogmod_log_Phi(a + c);\n    real br = 1 + exp(lD1 - lP) - exp(lD2 - lP);\n    if (br <= 0) return [0, negative_infinity()]';\n    real lS = fmin(lP + log(br), 0);\n    return [lS < 0 ? log1m_exp(lS) : negative_infinity(), lS]';\n  }\n  real lQ = cogmod_log_Phi(-a);\n  real R = exp(cogmod_log_Phi(-a - c) - lQ);\n  real br = R - (1 - R) / A + exp(lD2 - lQ);\n  if (br <= 0) return [negative_infinity(), 0]';\n  real lF = fmin(lQ + log(br), 0);\n  return [lF, lF < 0 ? log1m_exp(lF) : negative_infinity()]';\n}\n\n// Above A = 0 the two tails share lD1 and lD2, so building the pair and taking\n// one of them is the cheap way round. At A = 0 they share nothing - each is a\n// single cogmod_log_Phi() of the same standardized time - and the pair would\n// put a whole discarded tail on the autodiff tape for every observation.\n// cogmod_lnr() reads the survival alone, once per trial per loser, so that is\n// the hot path of the family.\nreal cogmod_lognormal_acc_logcdf(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return cogmod_log_Phi((log(t) - meanlog) / sigma);\n  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[1];\n}\n\nreal cogmod_lognormal_acc_logsurv(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return cogmod_log_Phi((meanlog - log(t)) / sigma);\n  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[2];\n}\n\n// Log-likelihood for one observation from the shifted LogNormal model.\n// Y: observed reaction time.\n// mu: mean of the decision time on the log scale (meanlog).\n// sigma: SD of the decision time on the log scale (> 0).\n// sigmabias: start-point range, in units of the threshold offset (>= 0); 0 is the plain LogNormal.\n// ndt: non-decision time, same unit as Y (> 0).\n// poutlier: proportion of responses from the outlier process, in [0, 1].\n//\n// The outlier component is a half Normal with scale 0.2 s. It keeps the density\n// strictly positive below `ndt`, where the shifted decision component has none.\n// That is what removes the hard min-RT boundary and lets `ndt` be estimated\n// directly rather than as a fraction of an observed minimum. The scale is a\n// constant in SECONDS. This family expects reaction times in seconds; give it\n// another unit and the component contributes nothing anywhere in the data,\n// which silently reinstates the min-RT boundary it exists to remove.\n//\n// It is written out rather than called as normal_lpdf() because both of its\n// parameters are constant: written that way, Stan recomputes the normalising\n// constant for every observation on every leapfrog step.\nreal cogmod_lognormal_lpdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    // Parameter checks\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n\n    // The leading constant includes the log(2) that folds the symmetric\n    // Normal onto [0, Inf).\n    real lp_out = 1.3836465597893728 - 12.5 * square(Y);\n    real t_adj  = Y - ndt;\n\n    // Faster than the non-decision time: only the outlier component can have\n    // produced this response.\n    if (t_adj <= 0) return log(poutlier) + lp_out;\n\n    return log_mix(poutlier, lp_out, cogmod_lognormal_acc_ldens(t_adj, mu, sigma, sigmabias));\n}\n\n// Log CDF and log survival of the same mixture, for brms's cens() addition\n// term. See ?rcogmod_invgaussian for what censoring a reaction time means and\n// when it is the right model.\nreal cogmod_lognormal_lcdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n    // Outlier CDF: 2 Phi(Y / s) - 1 = erf(Y / (s sqrt(2)))\n    real lF_out = log(erf(Y * 3.5355339059327369));\n    real t_adj  = Y - ndt;\n    if (t_adj <= 0) return log(poutlier) + lF_out;\n    return log_mix(poutlier, lF_out, cogmod_lognormal_acc_logcdf(t_adj, mu, sigma, sigmabias));\n}\n\nreal cogmod_lognormal_lccdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return 0;\n    // Outlier survival: 2 Phi(-Y / s), through the lower tail (see above)\n    real lS_out = 0.69314718055994529 + std_normal_lcdf(-Y * 5);\n    real t_adj  = Y - ndt;\n    // Not yet past the non-decision time: the decision process cannot have\n    // finished, so its survival is exactly 1.\n    if (t_adj <= 0) return log_mix(poutlier, lS_out, 0);\n    return log_mix(poutlier, lS_out, cogmod_lognormal_acc_logsurv(t_adj, mu, sigma, sigmabias));\n}\n"
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
#> [1] "\n// log(Phi(x)), the one piece of arithmetic every normal tail in this package is\n// built from. Neither of Stan's two routes to it is good enough for both jobs\n// it has here, which is to be right in the far tail *and* to hand back a\n// usable derivative there. All three claims below were measured against\n// central differences of the log probability over 20000 responses.\n//\n// The erfc route - std_normal_lcdf() is not it, but lognormal_lcdf(),\n// lognormal_lccdf() and the log(u1) + log1m(u2 / u1) the LogNormal used to\n// write are - has good partials, to about 4e-6 on a summed gradient of order\n// 1e3. But erfc underflows near x = -38, and then the value is log(0) and the\n// partials are inf or 0/0. That is not a harmless -inf in a mixture:\n// log_mix() in the lpdf stays finite when the decision component is -inf, but\n// reverse mode multiplies the (zero) adjoint into the stored partial, and\n// 0 * inf is NaN, so a single response turns the gradient of the whole model\n// to NaN - 'Gradient evaluated at the initial value is not finite' at the\n// start of a fit, divergent transitions afterwards.\n//\n// std_normal_lcdf() has the range: its value is exact against R's\n// pnorm(log.p = TRUE) as far as x = -1e7. Its partials are not - they sat\n// 1.7e-3 from central differences where the erfc route sat 4e-6 on the LNR,\n// and 2e-4 to 7e-2 on the RDM, which took every tail through it - so it is not\n// a drop-in for the tails of a race, where those partials are the gradient of\n// the drifts and the scales.\n//\n// So: erfc in the body of the distribution, and below x = -25 the asymptotic\n// expansion of the tail,\n//\n//   Phi(x) = phi(x) / (-x) * (1 - 1/x^2 + 3/x^4 - 15/x^6 + 105/x^8 - 945/x^10)\n//\n// whose leading term is the exponent itself. Nothing underflows, the result\n// stays finite and differentiable as far as x = -1e150, and the six terms\n// agree with pnorm(log.p = TRUE) to 4e-16 relative from x = -25 down - the\n// last bit of a double - so the two branches meet with no step in the density.\nreal cogmod_log_Phi(real x) {\n  if (x < -25) {\n    real z = inv_square(x);\n    real series = 1 + z * (-1 + z * (3 + z * (-15 + z * (105 - 945 * z))));\n    return -0.5 * square(x) - log(-x) - 0.91893853320467274 + log(series);\n  }\n  if (x > 0) return log1p(-0.5 * erfc(x * 0.7071067811865476));\n  return log(0.5 * erfc(-x * 0.7071067811865476));\n}\n\n// log(Phi(y + c) - Phi(y)) for c > 0, from whichever tail keeps the two terms\n// from cancelling: the upper one when y > 0, where both CDFs sit near 1, the\n// lower one otherwise. Taken as the larger tail plus log(1 - ratio) in log\n// space, the way .lognormal_ldiff_pnorm() does it - the quotient u2 / u1 this\n// replaces divides two minute numbers and lost its own accuracy long before\n// either underflowed (8e-4 against central differences, against 1e-6 here).\nreal cogmod_lognormal_ldiff_Phi(real y, real c) {\n  real hi;\n  real lo;\n  if (y > 0) {\n    hi = cogmod_log_Phi(-y);\n    lo = cogmod_log_Phi(-(y + c));\n  } else {\n    hi = cogmod_log_Phi(y + c);\n    lo = cogmod_log_Phi(y);\n  }\n  return lo < hi ? hi + log1m_exp(lo - hi) : negative_infinity();\n}\n\n// Log density of the accumulator's finishing time with start-point range A.\n// At A = 0 this is the LogNormal itself, at the LogNormal's cost.\nreal cogmod_lognormal_acc_ldens(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return lognormal_lpdf(t | meanlog, sigma);\n  real a = (meanlog - log(t)) / sigma;\n  real c = log1p(A) / sigma;\n  real x = a - sigma;\n  if (c < 1e-4) {\n    real series = 1 - c * x / 2 + square(c) * (square(x) - 1) / 6;\n    if (series <= 0) return negative_infinity();\n    return lognormal_lpdf(t | meanlog, sigma) + log(log1p(A) / A) + log(series);\n  }\n  return -meanlog + square(sigma) / 2 + cogmod_lognormal_ldiff_Phi(x, c) - log(A);\n}\n\n// [log F, log S] of the accumulator's finishing time, each computed directly\n// on the side where it is the small one. See .lognormal_acc_ltails().\nvector cogmod_lognormal_acc_ltails(real t, real meanlog, real sigma, real A) {\n  real a = (meanlog - log(t)) / sigma;\n  // At A = 0 this is the plain LogNormal, whose two tails are Phi(-a) and\n  // Phi(a). Written that way rather than as lognormal_lcdf()/lognormal_lccdf(),\n  // which are erfc alone and so reach log(0) with non-finite partials around\n  // |a| = 38 - see cogmod_log_Phi() above for what that costs.\n  if (A == 0) return [cogmod_log_Phi(-a), cogmod_log_Phi(a)]';\n  real c = log1p(A) / sigma;\n  if (c < 1e-4) {\n    real r = log1p(A) / A;\n    real corr = c * r * (0.5 + c * (2 * sigma - a) / 6);\n    real lPa = cogmod_log_Phi(a);\n    real lQa = cogmod_log_Phi(-a);\n    real lphi = std_normal_lpdf(a);\n    real lS = fmin(lPa + log1p(corr * exp(lphi - lPa)), 0);\n    real dF = 1 - corr * exp(lphi - lQa);\n    real lF = dF > 0 ? fmin(lQa + log(dF), 0) : negative_infinity();\n    return [lF, lS]';\n  }\n  real lA = log(A);\n  real lD1 = cogmod_lognormal_ldiff_Phi(a, c) - lA;\n  real lD2 = -meanlog + square(sigma) / 2 + log(t)\n             + cogmod_lognormal_ldiff_Phi(a - sigma, c) - lA;\n  if (a < 0) {\n    real lP = cogmod_log_Phi(a + c);\n    real br = 1 + exp(lD1 - lP) - exp(lD2 - lP);\n    if (br <= 0) return [0, negative_infinity()]';\n    real lS = fmin(lP + log(br), 0);\n    return [lS < 0 ? log1m_exp(lS) : negative_infinity(), lS]';\n  }\n  real lQ = cogmod_log_Phi(-a);\n  real R = exp(cogmod_log_Phi(-a - c) - lQ);\n  real br = R - (1 - R) / A + exp(lD2 - lQ);\n  if (br <= 0) return [negative_infinity(), 0]';\n  real lF = fmin(lQ + log(br), 0);\n  return [lF, lF < 0 ? log1m_exp(lF) : negative_infinity()]';\n}\n\n// Above A = 0 the two tails share lD1 and lD2, so building the pair and taking\n// one of them is the cheap way round. At A = 0 they share nothing - each is a\n// single cogmod_log_Phi() of the same standardized time - and the pair would\n// put a whole discarded tail on the autodiff tape for every observation.\n// cogmod_lnr() reads the survival alone, once per trial per loser, so that is\n// the hot path of the family.\nreal cogmod_lognormal_acc_logcdf(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return cogmod_log_Phi((log(t) - meanlog) / sigma);\n  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[1];\n}\n\nreal cogmod_lognormal_acc_logsurv(real t, real meanlog, real sigma, real A) {\n  if (A == 0) return cogmod_log_Phi((meanlog - log(t)) / sigma);\n  return cogmod_lognormal_acc_ltails(t, meanlog, sigma, A)[2];\n}\n\n// Log-likelihood for one observation from the shifted LogNormal model.\n// Y: observed reaction time.\n// mu: mean of the decision time on the log scale (meanlog).\n// sigma: SD of the decision time on the log scale (> 0).\n// sigmabias: start-point range, in units of the threshold offset (>= 0); 0 is the plain LogNormal.\n// ndt: non-decision time, same unit as Y (> 0).\n// poutlier: proportion of responses from the outlier process, in [0, 1].\n//\n// The outlier component is a half Normal with scale 0.2 s. It keeps the density\n// strictly positive below `ndt`, where the shifted decision component has none.\n// That is what removes the hard min-RT boundary and lets `ndt` be estimated\n// directly rather than as a fraction of an observed minimum. The scale is a\n// constant in SECONDS. This family expects reaction times in seconds; give it\n// another unit and the component contributes nothing anywhere in the data,\n// which silently reinstates the min-RT boundary it exists to remove.\n//\n// It is written out rather than called as normal_lpdf() because both of its\n// parameters are constant: written that way, Stan recomputes the normalising\n// constant for every observation on every leapfrog step.\nreal cogmod_lognormal_lpdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    // Parameter checks\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n\n    // The leading constant includes the log(2) that folds the symmetric\n    // Normal onto [0, Inf).\n    real lp_out = 1.3836465597893728 - 12.5 * square(Y);\n    real t_adj  = Y - ndt;\n\n    // Faster than the non-decision time: only the outlier component can have\n    // produced this response.\n    if (t_adj <= 0) return log(poutlier) + lp_out;\n\n    return log_mix(poutlier, lp_out, cogmod_lognormal_acc_ldens(t_adj, mu, sigma, sigmabias));\n}\n\n// Log CDF and log survival of the same mixture, for brms's cens() addition\n// term. See ?rcogmod_invgaussian for what censoring a reaction time means and\n// when it is the right model.\nreal cogmod_lognormal_lcdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return negative_infinity();\n    // Outlier CDF: 2 Phi(Y / s) - 1 = erf(Y / (s sqrt(2)))\n    real lF_out = log(erf(Y * 3.5355339059327369));\n    real t_adj  = Y - ndt;\n    if (t_adj <= 0) return log(poutlier) + lF_out;\n    return log_mix(poutlier, lF_out, cogmod_lognormal_acc_logcdf(t_adj, mu, sigma, sigmabias));\n}\n\nreal cogmod_lognormal_lccdf(real Y, real mu, real sigma, real sigmabias, real ndt, real poutlier) {\n    if (sigma <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {\n      return negative_infinity();\n    }\n    if (Y <= 0) return 0;\n    // Outlier survival: 2 Phi(-Y / s), through the lower tail (see above)\n    real lS_out = 0.69314718055994529 + std_normal_lcdf(-Y * 5);\n    real t_adj  = Y - ndt;\n    // Not yet past the non-decision time: the decision process cannot have\n    // finished, so its survival is exactly 1.\n    if (t_adj <= 0) return log_mix(poutlier, lS_out, 0);\n    return log_mix(poutlier, lS_out, cogmod_lognormal_acc_logsurv(t_adj, mu, sigma, sigmabias));\n}\n"
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
