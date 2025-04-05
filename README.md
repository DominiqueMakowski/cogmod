# cogmod


Cognitive Models for Subjective Scales and Decision Making Tasks in R

## Status

**This package is very much totally exploratory - currently made for my
own needs.** It’s not meant to be stable and robust at this stage. Use
at your own risks.

**Please get in touch if you have suggestions for improvement, I am new
to this!**

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
