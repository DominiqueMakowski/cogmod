# Warm-start a fit from a previous one: metric, step size and starting values

Takes what warmup produced in a previous fit - the adapted inverse
metric, the step size and the posterior means - and turns it into the
`inv_metric`, `step_size` and `init` arguments of a new
[`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html) call,
so that the new run can get by with a much shorter warmup. The previous
fit can be the **same model** (a refit with more draws, another seed, a
slightly different prior) or a **pilot on a subset of the
participants**: the full model then has more parameters, one
standardized random effect per new participant per group-level term, and
the metric is carried over **by parameter name**, with a sensible filler
for what the pilot never saw.

    ws <- cogmod_warmstart(pilot, data = data)   # the pilot's model, on all the data
    m <- brm(formula, data = data, prior = ..., stanvars = ...,
             init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size,
             warmup = 100, iter = 600, backend = "cmdstanr")

The sampler keeps adapting from the supplied values during whatever
warmup remains, so a poor warm start costs speed, not correctness.
Nothing about it changes the posterior being sampled.

## Usage

``` r
cogmod_warmstart(x, formula = NULL, data = NULL, jitter = 0.05, ...)

cogmod_inv_metric(formula = NULL, data = NULL, warmstart, ...)

cogmod_step_size(formula = NULL, data = NULL, warmstart, ...)

# S3 method for class 'cogmod_warmstart'
as.data.frame(x, row.names = NULL, optional = FALSE, ...)

# S3 method for class 'cogmod_warmstart'
print(x, ...)
```

## Arguments

- x:

  The source: a `brmsfit` fitted with `backend = "cmdstanr"`, a
  `cogmod_warmstart` object, the data frame
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) makes
  of one, or the path to a CSV file holding that data frame.

- formula, data:

  The target model, as they will be passed to
  [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html).
  Either left `NULL` (the default) is taken from the source fit: the
  same formula on new data, the same data under a new formula, or with
  both `NULL` the source fit itself. That requires `x` to be a
  `brmsfit`; a table or file source needs both.

- jitter:

  SD of the noise added to the starting values on the unconstrained
  scale, so that chains start at different points. Smaller than
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)'s
  default because the values come from a converged posterior; `0` gives
  identical starts.

- ...:

  Passed to
  [`brms::make_stancode()`](https://paulbuerkner.com/brms/reference/stancode.html),
  [`brms::make_standata()`](https://paulbuerkner.com/brms/reference/standata.html)
  and [`brms::brm()`](https://paulbuerkner.com/brms/reference/brm.html)
  (with `empty = TRUE`) when the target model is built, for arguments
  such as `data2`.

- warmstart:

  The source, as `x` above: a `brmsfit`, a `cogmod_warmstart` object,
  its data frame, or the path to a CSV file of it.

- row.names, optional:

  Ignored; present for compatibility with the
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  generic.

## Value

An object of class `cogmod_warmstart`: a list with

- `inv_metric`:

  Numeric vector, one variance per unconstrained parameter of the target
  model, in Stan's order.

- `step_size`:

  The step size, averaged over the source's chains.

- `init`:

  A function of one argument, for `brms::brm(init = )`, returning a
  named list of starting values for every parameter the target program
  declares.

- `table`:

  A data frame with one row per unconstrained parameter: `parameter`
  (the Stan label), `group`, `coef` and `level` (for a group-level
  parameter, the grouping factor, the coefficient and, for a
  standardized effect, the level it stands for; otherwise `NA`),
  `inv_metric`, `mean` (the source posterior mean, `NA` where none
  applies) and `step_size` (the same value in every row). This is what
  [`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html)
  returns.

- `counts`, `missing`:

  How many entries came from the source, are new group levels, or have
  no counterpart, and the names of the latter; what
  [`print()`](https://rdrr.io/r/base/print.html) reports.

## What is carried over, and how

Stan adapts one variance per **unconstrained** scalar parameter, in the
order of the program's `parameters` block, and `brms` keeps those
variances (one vector per chain, averaged here) and the step size in the
fit's metadata. To move them to another model they are first labelled
with the Stan parameter names - `Intercept`, `sd_1[1]`, `z_1[1,3]`, and
so on - which are read off the generated program, the same way
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
does it. The target model's labels are built the same way, and the two
are joined:

- **Population-level entries** (coefficients, intercepts, group-level
  SDs and correlations, auxiliary parameters) carry over one to one.

- **Standardized group-level effects** (`z_<k>[m, n]`) are indexed by
  participant, so a pilot participant's entry moves to that
  participant's position in the target model, matched by the level's
  *name*. Participants the pilot never saw get the average of the
  pilot's entries for the same term (they are standardized effects, so
  that is a fair guess) and start at zero.

- **Anything without a counterpart** - a predictor or a group-level term
  the pilot did not have - gets Stan's default variance of 1 and the
  generic starting value
  [`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
  would give it. A count is reported by
  [`print()`](https://rdrr.io/r/base/print.html); a large one usually
  means the two formulas differ more than intended.

Starting values are the pilot's posterior means. `brms` drops the raw
`z_` and Cholesky factors from a saved fit, so these are rebuilt from
what it keeps: the group-level effects `r_`, their SDs and their
correlations. Each chain gets the same values plus a little noise on the
unconstrained scale (`jitter`), so that the chains do not start at one
point and `Rhat` keeps some meaning. Note that tightly initialised
chains are less likely to find a second mode than dispersed ones; if
that is a concern, run the cold start once.

## Same model, or a different one

Whatever is not given is taken from the source fit. With `formula` and
`data` both left `NULL`, the target is the fit itself: the result
reproduces its own adaptation, and is the way to refit with fewer warmup
iterations. With only `data`, the target is the same model on that
data - the pilot-to-full-sample case. With only `formula`, it is that
model on the source's data - a variant of the model, say with one more
predictor, whose shared parameters can start where the first fit left
them. A table or file source carries neither and needs both. Any `brms`
model fitted with the `cmdstanr` backend and the default diagonal metric
can be a source; a `dense_e` fit is refused, and so is the `rstan`
backend, which does not store the adaptation.

The `inv_metric` returned is a plain vector, as `cmdstanr` wants it; the
labels are in `ws$table`. The metric must match the target program
exactly, which is why `formula` and `data` are needed rather than just a
count of participants: `brms` decides the layout from both.

## Storing it, and the three helpers

[`as.data.frame()`](https://rdrr.io/r/base/as.data.frame.html) gives a
small table (one row per unconstrained parameter: label, group and
level, variance, posterior mean, step size) that can be written with
[`utils::write.csv()`](https://rdrr.io/r/utils/write.table.html), so
that a pilot fitted on a laptop can warm-start an array job on a cluster
with a file of a few kilobytes and no `brmsfit` in sight. On the other
side, three helpers with the signature of
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
and
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md) -
the model's formula and data first, the source under `warmstart` - each
give one argument of the
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html) call:

    tab <- read.csv("pilot_warmstart.csv")   # or the path, or the brmsfit itself
    m <- brm(formula, data = data, prior = ..., stanvars = ...,
             init = cogmod_inits(formula, data, warmstart = tab),
             inv_metric = cogmod_inv_metric(formula, data, warmstart = tab),
             step_size = cogmod_step_size(formula, data, warmstart = tab),
             warmup = 100, iter = 600, backend = "cmdstanr")

Each maps the table onto the model `formula` and `data` describe, so it
does not matter which model the table was written for: a table from a
pilot on fewer participants is extended, one written for another formula
falls back to the defaults with a note. When the table was made for this
very model, `tab$inv_metric` and `tab$step_size[1]` are the same numbers
(the step size is one number repeated down the column; a whole column
there would be read as one step size per chain).

## What it is worth

On an LNR and a DDM with participant random intercepts, a pilot on 4 of
8 participants warm-started the full fit to about twice the effective
draws per second of a cold start with a 500-iteration warmup, and four
to six times those of a cold start with the same 100-iteration warmup.
The starting values alone bought nothing: what a short warmup lacks is
an adapted step size and metric, not a good position. Details and the
benchmark are in `vignette("performance")`.

## See also

[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md),
which supplies the starting values used where the source has none, and
`vignette("performance")`.

## Examples

``` r
if (FALSE) { # \dontrun{
# A pilot on some participants, then the full sample
f <- brms::bf(RT | dec(choice) ~ Condition + (1 | id), ndt ~ 1 + (1 | id),
              poutlier ~ 1, family = cogmod_lnr())
pilot <- brms::brm(f, data = df[df$id %in% 1:4, ], prior = cogmod_priors(f, df),
                   init = cogmod_inits(f, df), stanvars = cogmod_stanvars(f),
                   backend = "cmdstanr")
ws <- cogmod_warmstart(pilot, data = df)   # same model, all the data
ws   # how much of the metric came from the pilot
m <- brms::brm(f, data = df, prior = cogmod_priors(f, df), stanvars = cogmod_stanvars(f),
               init = ws$init, inv_metric = ws$inv_metric, step_size = ws$step_size,
               warmup = 100, iter = 600, backend = "cmdstanr")

# Keep it for a cluster
write.csv(as.data.frame(ws), "pilot_warmstart.csv", row.names = FALSE)
# ... and there, one helper per argument, all with the same signature
m <- brms::brm(f, data = df, prior = cogmod_priors(f, df), stanvars = cogmod_stanvars(f),
               init = cogmod_inits(f, df, warmstart = "pilot_warmstart.csv"),
               inv_metric = cogmod_inv_metric(f, df, warmstart = "pilot_warmstart.csv"),
               step_size = cogmod_step_size(f, df, warmstart = "pilot_warmstart.csv"),
               warmup = 100, iter = 600, backend = "cmdstanr")
} # }
```
