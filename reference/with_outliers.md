# Include or exclude the outlier component in predictions

Switches the `predict_outliers` flag on a model fitted with
[`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md)
or
[`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md),
controlling whether
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
describe the fitted mixture or the decision process alone.

Predictions **exclude** the outlier component by default, because for
almost every downstream use it is a nuisance: it pulls expected values
toward its own mean (0.16 s) and adds a spike of implausibly fast draws
to posterior predictive samples. It is also a deliberately fixed
regularizer rather than a claim about how guesses are distributed, so
simulating from it means simulating from something the model does not
assert.

`with_outliers()` restores the mixture. The main reason to want it is
[`brms::pp_check()`](https://mc-stan.org/bayesplot/reference/pp_check.html):
on untrimmed data the decision-only predictive has no fast spike to
match the one in the data, which reads as misfit. Use
`pp_check(with_outliers(m))` for a like-for-like check.

The flag is stored **on the model** rather than passed as an argument,
because `brms` and the packages built on it (`insight`, `modelbased`,
`marginaleffects`, `emmeans`) do not forward extra arguments down to a
custom family's prediction methods -
[`posterior_epred()`](https://mc-stan.org/rstantools/reference/posterior_epred.html)
reaches the family method with `prep` and nothing else. Carrying it on
the object is what makes it work through all of them. The same flag can
be set up front with `cogmod_lognormal(predict_outliers = TRUE)`.

[`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html) is
unaffected and has no equivalent switch: the likelihood *is* the
mixture, and dropping a component from it would not be a different
summary of the same model but a different model. One consequence worth
knowing is that
[`posterior_predict()`](https://mc-stan.org/rstantools/reference/posterior_predict.html)
and [`log_lik()`](https://mc-stan.org/rstantools/reference/log_lik.html)
do not describe the same distribution by default. This also desyncs
`loo_pit()`,
[`loo_predict()`](https://mc-stan.org/rstantools/reference/loo-prediction.html)
and
[`bayes_R2()`](https://mc-stan.org/rstantools/reference/bayes_R2.html)
from [`loo()`](https://mc-stan.org/loo/reference/loo.html), not just
hand-rolled checks - anything that compares a simulated replicate
against the likelihood should be run on `with_outliers()`.

## Usage

``` r
with_outliers(object)

without_outliers(object)
```

## Arguments

- object:

  A `brmsfit` fitted with
  [`cogmod_lognormal()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lognormal.md),
  [`cogmod_loggamma()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_loggamma.md)
  or any other family built on the outlier mixture - see the *Supported
  families* section of
  [`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
  for the full list, which includes the choice-and-RT families such as
  [`cogmod_lnr()`](https://dominiquemakowski.github.io/cogmod/reference/rcogmod_lnr.md).

## Value

The model, with the flag set. The fit itself is untouched - only how
predictions are summarised changes.

## Examples

``` r
# f <- bf(RT ~ Condition, ndt ~ 1, poutlier ~ 1, family = cogmod_lognormal())
# m <- brms::brm(f, data = df, stanvars = cogmod_stanvars(f))
#
# # the decision process alone - the default, everywhere downstream
# brms::posterior_epred(m)
# modelbased::estimate_means(m, by = "Condition")
# marginaleffects::avg_predictions(m, by = "Condition")
#
# # the fitted mixture, e.g. for a like-for-like predictive check
# brms::pp_check(with_outliers(m))
# without_outliers(m)  # back to the default
```
