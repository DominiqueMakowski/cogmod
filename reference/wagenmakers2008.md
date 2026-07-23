# Lexical Decision Data from Wagenmakers et al. (2008)

Response times and accuracy from a lexical decision task, Experiment 1
in [Wagenmakers et al.
(2008)](https://doi.org/10.1016/j.jml.2007.04.006) (also reanalyzed in
[Heathcote & Love, 2012](https://doi.org/10.3389/fpsyg.2012.00292)). Six
participants judged whether letter strings were valid English words
under instructions that emphasized either speed or accuracy.

## Usage

``` r
wagenmakers2008
```

## Format

A data frame with 10,369 rows and 5 variables:

- Participant:

  Participant identifier (integer, 1-6).

- Condition:

  Instruction condition, `"Speed"` or `"Accuracy"`.

- RT:

  Response time, in seconds.

- Error:

  Logical, whether the response was an error.

- Frequency:

  Word frequency category of the presented stimulus.

## Source

<https://github.com/DominiqueMakowski/CognitiveModels/blob/main/data/wagenmakers2008.csv>

## Examples

``` r
data(wagenmakers2008)
head(wagenmakers2008)
#>   Participant Condition    RT Error Frequency
#> 1           1     Speed 0.700 FALSE       Low
#> 2           1     Speed 0.392  TRUE  Very Low
#> 3           1     Speed 0.460 FALSE  Very Low
#> 4           1     Speed 0.455 FALSE  Very Low
#> 5           1     Speed 0.505  TRUE       Low
#> 6           1     Speed 0.773 FALSE      High
```
