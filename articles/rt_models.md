# RT-only Models

``` r

library(cogmod)
library(easystats)
library(ggplot2)
library(brms)
library(cmdstanr)

options(mc.cores = parallel::detectCores() - 2)
```

## The Data

For this chapter, we will be using the data from [Wagenmakers et al.,
(2008)](https://doi.org/10.1016/j.jml.2007.04.006) - Experiment 1 also
reanalyzed by [Heathcote & Love
(2012)](https://doi.org/10.3389/fpsyg.2012.00292), that contains
responses and response times for several participants in two conditions
(where instructions emphasized either **speed** or **accuracy**). This
dataset is bundled directly with `cogmod` as `wagenmakers2008`. Using
the same procedure as the authors, we excluded all trials with
uninterpretable response time, i.e., responses that are too fast (\<200
ms instead of \<180 ms) or too slow (\>2 sec instead of \>3 sec).

``` r

set.seed(123)  # For reproducibility

df <- cogmod::wagenmakers2008
df <- df[df$RT > 0.2 & df$Participant %in% c(1, 2, 3), ]

# Show 10 first rows
head(df, 10)
#>    Participant Condition    RT Error Frequency
#> 1            1     Speed 0.700 FALSE       Low
#> 2            1     Speed 0.392  TRUE  Very Low
#> 3            1     Speed 0.460 FALSE  Very Low
#> 4            1     Speed 0.455 FALSE  Very Low
#> 5            1     Speed 0.505  TRUE       Low
#> 6            1     Speed 0.773 FALSE      High
#> 7            1     Speed 0.390 FALSE      High
#> 8            1     Speed 0.587  TRUE       Low
#> 9            1     Speed 0.603 FALSE       Low
#> 10           1     Speed 0.435 FALSE      High
```

We are going to first take interest in the response times (RT) of
**Correct** answers only (as we can assume that errors are underpinned
by a different *generative process*).

``` r

df <- df[df$Error == 0, ]
```

``` r

ggplot(df, aes(x = RT, fill = Condition)) +
  geom_histogram(bins = 120, alpha = 0.8, position = "identity") +
  scale_fill_manual(values = c("darkgreen", "darkred")) +
  theme_minimal()
```

![](rt_models_files/figure-html/unnamed-chunk-3-1.png)

## Models

### Normal

A basic linear model.

Code

``` r

f <- bf(
  RT ~ Condition
)

m_normal <- brm(f,
  data = df, 
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_normal <- brms::add_criterion(m_normal, "loo") 

saveRDS(m_normal, file = "models/m_normal.rds")
```

### ExGaussian

The Ex-Gaussian distribution is a popular choice for RT data, as it
separates the central portion of the distribution - captured by the
Gaussian parameters `mu` and `sigma` - from the positively skewed tail,
captured by the exponential parameter `tau`. However, `brms`’s native
[`exgaussian()`](https://paulbuerkner.com/brms/reference/brmsfamily.html)
family does **not** use this “classical” parameterization familiar to
experimental psychologists: its `mu` indexes the mean of the *entire*
distribution (i.e., Gaussian + exponential combined) rather than the
location of the Gaussian component alone. This matters because a change
in the Gaussian location and an opposite change in the exponential tail
can cancel out at the level of the overall mean, so effects estimated on
`brms`’s default `mu` can lead to different (and potentially incorrect)
inferences about the underlying process than effects estimated on the
classical `mu`.

For this reason, `cogmod` provides its own
[`rt_exgaussian()`](https://github.com/DominiqueMakowski/cogmod/reference/rt_exgaussian.md)
custom family (internally relying on Stan’s `exp_mod_normal`
distribution), in which `mu` and `sigma` are the mean and SD of the
Gaussian component and `tau` is the mean of the exponential tail -
directly matching the classical parameterization.

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  family = rt_exgaussian()
)

m_exgauss <- brm(f,
  data = df, 
  stanvars = rt_exgaussian_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_exgauss <- brms::add_criterion(m_exgauss, "loo") 

saveRDS(m_exgauss, file = "models/m_exgauss.rds")
```

### Shifted LogNormal

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_lognormal()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_lognormal <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_lognormal_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_lognormal <- brms::add_criterion(m_lognormal, "loo") 

saveRDS(m_lognormal, file = "models/m_lognormal.rds")
```

### Inverse Gaussian (Shifted Wald)

The Shifted Wald model (also known as the Inverse Gaussian distribution)
is actually equivalent to a one-response version of the Drift Diffusion
Model (DDM) with no between-trial variability in drift rate, starting
point, or non-decision time.

Code

``` r

f <- bf(
  RT ~ Condition,
  bs ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_invgaussian()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_wald <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_invgaussian_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_wald <- brms::add_criterion(m_wald, "loo") 

saveRDS(m_wald, file = "models/m_wald.rds")
```

### Weibull

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_weibull()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_weibull <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_weibull_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_weibull <- brms::add_criterion(m_weibull, "loo") 

saveRDS(m_weibull, file = "models/m_weibull.rds")
```

### LogWeibull (Shifted Gumbel)

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_logweibull()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_logweibull <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_logweibull_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_logweibull <- brms::add_criterion(m_logweibull, "loo") 

saveRDS(m_logweibull, file = "models/m_logweibull.rds")
```

### Inverse Weibull (Shifted Fréchet)

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_invweibull()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_invweibull <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_invweibull_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_invweibull <- brms::add_criterion(m_invweibull, "loo") 

saveRDS(m_invweibull, file = "models/m_invweibull.rds")
```

### Gamma

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_gamma()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_gamma <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_gamma_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_gamma <- brms::add_criterion(m_gamma, "loo") 

saveRDS(m_gamma, file = "models/m_gamma.rds")
```

### Inverse Gamma

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_invgamma()
)

priors <- brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau") |>
  brms::validate_prior(f, data = df)

m_invgamma <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_invgamma_stanvars(),
  init = 0,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_invgamma <- brms::add_criterion(m_invgamma, "loo") 

saveRDS(m_invgamma, file = "models/m_invgamma.rds")
```

### Linear Ballistic Accumulator

Code

``` r

f <- bf(
  RT ~ Condition,
  sigma = 1,
  sigmabias ~ Condition,
  bs ~ Condition,
  tau ~ Condition,
  minrt = min(df$RT),
  family = rt_lba()
)

priors <- c(
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "mu"),
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "sigmabias"),
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "bs"),
  brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau")
  )|>
  brms::validate_prior(f, data = df)

m_lba <- brm(
  f,
  prior = priors,
  data = df, 
  stanvars = rt_lba_stanvars(),
  init = 0.5,
  chains = 4, iter = 500, backend = "cmdstanr"
)

m_lba <- brms::add_criterion(m_lba, "loo") 

saveRDS(m_lba, file = "models/m_lba.rds")
```

## Model Comparison

### Model Fit

We can compare these models together using the `loo` package, which
shows that CHOCO provides a significantly better fit than the other
models.

``` r

loo::loo_compare(m_normal, m_exgauss, m_lognormal, m_wald,
                 m_weibull, m_logweibull, m_invweibull,
                 m_gamma, m_invgamma, m_lba) |>
  parameters(include_ENP = TRUE)
#> Loading required namespace: rstan
#> # Fixed Effects
#> 
#> Name |   LOOIC |   ENP |    ELPD | Difference | Difference_SE |      p
#> ----------------------------------------------------------------------
#> 1    | -4845.1 |  7.29 | 2422.56 |       0.00 |          0.00 |       
#> 2    | -4840.8 |  4.72 | 2420.41 |      -2.15 |          5.95 | 0.717 
#> 3    | -4840.4 |  5.03 | 2420.18 |      -2.38 |          6.00 | 0.691 
#> 4    | -4831.7 |  8.31 | 2415.85 |      -6.71 |          6.22 | 0.281 
#> 5    | -4793.3 |  6.41 | 2396.64 |     -25.93 |          9.55 | 0.007 
#> 6    | -4763.0 |  8.80 | 2381.50 |     -41.06 |         11.08 | < .001
#> 7    | -4728.7 | 10.84 | 2364.36 |     -58.20 |         15.07 | < .001
#> 8    | -4393.2 |  7.53 | 2196.61 |    -225.95 |         24.80 | < .001
#> 9    | -3680.4 |  9.33 | 1840.19 |    -582.37 |         38.96 | < .001
#> 10   | -1904.2 |  7.82 |  952.09 |   -1470.47 |         73.27 | < .001
```

### Sampling Duration

Because each model was fit with only 4 chains, a boxplot of the
per-chain sampling times is not very informative. Instead, we summarize
each model’s sampling duration by the *median* time per chain, annotated
with its value.

As expected, the **Gaussian** model is by far the fastest to sample,
since it relies on `brms`’s built-in (and heavily optimized) Normal
likelihood with no custom Stan code or non-decision time shift involved.
At the other end, the **LBA** is by far the slowest, reflecting the
added cost of its multi-accumulator likelihood. The remaining RT-only
models (ExGaussian, LogNormal, Wald, Weibull, LogWeibull, InvWeibull,
Gamma, and InvGamma) are all relatively comparable to one another, as
they share a similar structure (a simple closed-form density combined
with a non-decision time shift).

``` r

duration <- rbind(
  data_modify(attributes(m_normal$fit)$metadata$time$chain, Model="Gaussian"),
  data_modify(attributes(m_exgauss$fit)$metadata$time$chain, Model="ExGaussian"),
  data_modify(attributes(m_lognormal$fit)$metadata$time$chain, Model="LogNormal"),
  data_modify(attributes(m_wald$fit)$metadata$time$chain, Model="Wald"),
  data_modify(attributes(m_weibull$fit)$metadata$time$chain, Model="Weibull"),
  data_modify(attributes(m_logweibull$fit)$metadata$time$chain, Model="LogWeibull"),
  data_modify(attributes(m_invweibull$fit)$metadata$time$chain, Model="InvWeibull"),
  data_modify(attributes(m_gamma$fit)$metadata$time$chain, Model="Gamma"),
  data_modify(attributes(m_invgamma$fit)$metadata$time$chain, Model="InvGamma"),
  data_modify(attributes(m_lba$fit)$metadata$time$chain, Model="LBA")
) |> 
  data_modify(Model = factor(Model, levels = c("Gaussian", "ExGaussian", "LogNormal", "Wald", "Weibull", "LogWeibull", "InvWeibull", "Gamma", "InvGamma", "LBA")))

duration_median <- aggregate(total ~ Model, data = duration, FUN = median)

duration_median |> 
  ggplot(aes(x = Model, y = total, fill = total)) +
  geom_col() +
  geom_text(aes(label = round(total, 1)), vjust = -0.5, size = 3.5) +
  labs(y = "Median Sampling Duration per Chain (s)", fill = "Duration") +
  scale_fill_gradientn(colors = c("turquoise", "green", "yellow", "gold", "orange", "red", "darkred")) +
  scale_y_sqrt() +
  theme_minimal() 
```

![](rt_models_files/figure-html/unnamed-chunk-16-1.png)

### Posterior Predictive Check

`iterations` controls the actual number of iterations used (e.g., for
the point-estimate) and `keep_iterations` the number included.

Code

``` r


pred <- rbind(
  # estimate_prediction(m_normal, keep_iterations = 50, iterations = 50) |>
  #   reshape_iterations() |>
  #   data_modify(Model = "Normal"),
  estimate_prediction(m_exgauss, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "ExGaussian"),
  estimate_prediction(m_lognormal, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "LogNormal"),
  estimate_prediction(m_wald, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "InvGaussian"),
  estimate_prediction(m_weibull, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "Weibull"),
  estimate_prediction(m_logweibull, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "LogWeibull"),
  estimate_prediction(m_invweibull, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "InvWeibull"),
  estimate_prediction(m_gamma, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "Gamma"),
  estimate_prediction(m_invgamma, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "InvGamma"),
  estimate_prediction(m_lba, keep_iterations = 50, iterations = 50) |>
    reshape_iterations() |>
    data_modify(Model = "LBA") |> 
    data_filter(iter_value < 2)
)

p <- pred |> 
  ggplot(aes(x=iter_value)) +
  geom_histogram(data = df, aes(x=RT, y = after_stat(density), fill = Condition), 
                 position = "identity", bins=120, alpha = 0.8) +
  geom_line(aes(color=Model, group=interaction(Condition, iter_group)), stat="density", alpha=0.2) +
  theme_minimal() +
  theme(axis.text.y = element_blank()) +
  facet_wrap(~Model) +
  coord_cartesian(xlim = c(0, 2)) +
  scale_fill_manual(values = c("darkgreen", "darkred")) +
  scale_color_material_d(guide = "none") +
  labs(x = "RT (s)", y = "Distribution")
p
```

![](rt_models_files/figure-html/unnamed-chunk-17-1.png)
