# Simulated Data Where Linear Models Fail

A simulated repeated-measures experiment in which the two conditions
have, **by construction, the same mean reaction time** while differing
radically in every other respect. Condition `A` has a short non-decision
time (50 ms) and a wide, heavily right-skewed distribution; condition
`B` has a long non-decision time (450 ms) and a narrow, nearly symmetric
one. Twenty participants each contribute 25 trials per condition.

## Usage

``` r
badlm
```

## Format

A data frame with 1,000 rows and 3 variables:

- Participant:

  Participant identifier (factor, `S01`-`S20`).

- Condition:

  Experimental condition, `"A"` or `"B"`.

- RT:

  Reaction time, in seconds.

## Source

Simulated; see `data-raw/badlm.R` for the generating code.

## Details

The dataset exists to demonstrate the limits of the summary-statistics
approach: a linear model - including a correctly specified linear
*mixed* model with a random intercept per participant - finds no effect
of `Condition`, because a difference in shift, in spread and in tail
weight is invisible to a comparison of means. It is used in the *RT
models* vignette and in the `cogmod` paper.

Participants differ in their overall speed through an *additive* offset
(SD = 30 ms) applied identically to both conditions, so that each
participant's true condition effect is exactly zero and a random
intercept is the correctly specified model for the between-participant
variation. The latent offsets are not included in the data.

## Examples

``` r
data(badlm)

# The two conditions have the same mean...
tapply(badlm$RT, badlm$Condition, mean)
#>         A         B 
#> 0.6940667 0.6928665 

# ...but nothing else in common
tapply(badlm$RT, badlm$Condition, sd)
#>          A          B 
#> 0.33150456 0.04731473 

# \donttest{
# A mixed model finds nothing
if (requireNamespace("lme4", quietly = TRUE)) {
  summary(lme4::lmer(RT ~ Condition + (1 | Participant), data = badlm))
}
#> Linear mixed model fit by REML ['lmerMod']
#> Formula: RT ~ Condition + (1 | Participant)
#>    Data: badlm
#> 
#> REML criterion at convergence: -32.3
#> 
#> Scaled residuals: 
#>     Min      1Q  Median      3Q     Max 
#> -2.3271 -0.3793 -0.0506  0.2039  6.1828 
#> 
#> Random effects:
#>  Groups      Name        Variance  Std.Dev.
#>  Participant (Intercept) 0.0005126 0.02264 
#>  Residual                0.0555791 0.23575 
#> Number of obs: 1000, groups:  Participant, 20
#> 
#> Fixed effects:
#>             Estimate Std. Error t value
#> (Intercept)  0.69407    0.01170   59.34
#> ConditionB  -0.00120    0.01491   -0.08
#> 
#> Correlation of Fixed Effects:
#>            (Intr)
#> ConditionB -0.637
# }
```
