# cogmod


Cognitive Models for Subjective Scales and Decision Making Tasks in R

## Status

**This package is very much totally exploratory - currently made for my
own needs.** It’s not meant to be stable and robust at this stage. Use
at your own risks.

- If you have suggestions for improvement, please get in touch!
- If you are interested in Sequential Sampling Models, see this amazing
  [Julia
  package](https://github.com/itsdfish/SequentialSamplingModels.jl)

## Installation

``` r
if(!requireNamespace("remotes", quietly = TRUE)) install.packages("remotes")

remotes::install_github("DominiqueMakowski/cogmod")
```

## CHOCO Model

The Choice-Confidence (CHOCO) model is useful to model data from
subjective scales, such as Likert-type or analog scales, in which the
left and the right side correspond to different processes or higher
order categorical responses (e.g., “disagree” vs. “agree”, “true”
vs. “false”). They can be used to jointly model choice (left or right)
and confidence (the degree of left or right).

<details class="code-fold">
<summary>Code</summary>

``` r
library(ggplot2)
library(cogmod)

# Simulate data using rchoco() with two parameter sets
df1 <- rchoco(n = 5000, mu = 0.5, muleft = 0.4, phileft = 3, kleft = 0.9)
df2 <- rchoco(n = 5000, mu = 0.7, muleft = 0.6, phileft = 5, kleft = 0.95)

# Combine data into a single data frame
df <- data.frame(
  value = c(df1, df2),
  group = rep(c("mu=0.5, muleft=0.4, phileft=3, kleft=0.9", 
                "mu=0.7, muleft=0.6, phileft=5, kleft=0.95"), each =5000)
)

# Create the histogram
ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.8, position = "identity", bins = 50) +
  labs(title = "CHOCO Distribution", x = "Value", y = "", fill = "Parameters") +
  theme_minimal() + 
  scale_fill_manual(values = c("#9C27B0", "#FF9800"))
```

</details>

![](man/figures/unnamed-chunk-2-1.png)

## LNR Model

The Log-Normal Race (LNR) model is useful for modeling reaction times
and errors in decision-making tasks. The model assumes that each
accumulator draws a value from a LogNormal distribution (shifted by a
non-decision time τ). The winning accumulator (minimum draw) determines
the observed reaction time and choice.

<details class="code-fold">
<summary>Code</summary>

``` r
# Simulate data using rlnr()
lnr_data <- rlnr(n = 5000, mu = 1, mud = 0.5, sigmazero = 1, sigmad = -0.5, tau = 0.2)

# Create histograms for each choice
ggplot(lnr_data, aes(x = rt, fill = factor(choice))) +
  geom_histogram(alpha = 0.8, position = "identity", bins = 50) +
  labs(title = "LNR Distribution", x = "Reaction Time", y = "Frequency", fill = "Choice") +
  theme_minimal() +
  scale_fill_manual(values = c("#4CAF50", "#FF5722"))
```

</details>

![](man/figures/unnamed-chunk-3-1.png)

## Usage with `brms`

### Simulate Data

``` r
options(mc.cores = parallel::detectCores() - 2)

library(brms)
library(cmdstanr)

df <- brms::rwiener(n = 5000, delta = 0.5, alpha = 1, beta = .3, tau = .25) |> 
  datawizard::data_rename(replacement = c("rt", "choice")) |> 
  datawizard::data_filter(rt < 2)

df |> 
  ggplot(aes(x = rt, fill = factor(choice))) +
  geom_histogram(alpha = 0.8, position = "identity", bins = 100) +
  labs(title = "RT Distribution", x = "Reaction Time", y = "Frequency", fill = "Choice") +
  theme_minimal() +
  scale_fill_manual(values = c("#009688", "#E91E63"))
```

![](man/figures/unnamed-chunk-4-1.png)

### Drift Diffusion Model (DDM)

``` r
f <- bf(rt | dec(choice) ~ 1,
        bs ~ 1,
        bias ~ 1,
        ndt ~ 1)

m_ddm <- brm(f, data = df, family = wiener(link = "identity", link_bs = "log", link_bias = "logit", link_ndt = "log"), 
             chains=4, iter = 500, backend="cmdstanr")

saveRDS(m_ddm, file = "man/figures/m_ddm.rds")
```

``` r
m_ddm <- readRDS("man/figures/m_ddm.rds")

parameters::parameters(m_ddm, component = "all")

# insight::get_predicted(m_ddm)
brms::posterior_predict(m_ddm, ndraws=5, newdata = df[1:4, ], negative_rt = TRUE)
```
