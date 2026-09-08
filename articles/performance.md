# Improving Sampling Efficiency and Performance

``` r

library(cogmod)
library(brms)
library(cmdstanr)
```

Evidence accumulation models are slow and expensive to sample, and a
fully mixed model with a comprehensive random effects structure can
easily take hours, days, or weeks for real-life datasets. Below are some
suggestions to help.

## The Model

Below is an example of a mixed DDM, with random intercepts and slopes on
most distributional parameter of interest. We also added a random
intercept to `poutlier`, with the assumption that the proportion of very
fast trials (anticipatory / inhibition failure responses vary across
participants).

``` r

f <- bf(
  RT | dec(Error) ~ Condition + (1 + Condition | Participant),
  boundary ~ Condition + (1 + Condition | Participant),
  bias ~ Condition + (1 + Condition | Participant),
  ndt ~ Condition + (1 + Condition | Participant),
  poutlier ~ 1 + (1 | Participant),
  sigmadrift = 0,
  sigmabias = 0,
  sigmandt = 0,
  family = cogmod_ddm()
)

m <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  chains = 4, iter = 2000, backend = "cmdstanr"
)
```

## Parallelization

The first step is to leverage the compounding benefits of two kinds of
parallelization:

- **Chain parallelization.** Run each of the `chains` on its own core
  via `backend = "cmdstanr"` and
  `options(mc.cores = parallel::detectCores())` (or the `cores` argument
  to [`brm()`](https://paulbuerkner.com/brms/reference/brm.html)). This
  is essentially free and should always be on.
- **Within-chain parallelization (multithreading).** `cmdstanr` can
  split a single chain’s likelihood evaluation across multiple threads
  with `threads = threading(n)`, passed to
  [`brm()`](https://paulbuerkner.com/brms/reference/brm.html) alongside
  `backend = "cmdstanr"`. This has real overhead: Stan has to partition
  the data and reduce the per-thread results back together, so it only
  pays off once the per-observation likelihood is expensive enough
  and/or the dataset is large enough that the reduction overhead is
  small relative to the work being split.

``` r

m <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  chains = 4,
  cores = 4,
  threads = threading(4)
)
```

Combining the two means a machine with, say, 16 cores can run 4 chains
with 4 threads each, or fewer chains with more threads each if warmup is
the bottleneck rather than the number of independent chains needed for
convergence diagnostics.

## Compiler Optimizations

Stan has optimizations that are off by default because they lengthen
compilation and, in one case, remove a safety net. None of them touches
the sampler or the posterior; they only make each gradient evaluation
cheaper. They are passed to
[`cmdstanr::cmdstan_model()`](https://mc-stan.org/cmdstanr/reference/cmdstan_model.html)
through [`brm()`](https://paulbuerkner.com/brms/reference/brm.html)’s
`stan_model_args`:

``` r

m <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  stan_model_args = list(
    stanc_options = list("O1"),
    cpp_options = list(STAN_CPP_OPTIMS = TRUE, STAN_NO_RANGE_CHECKS = TRUE)
  )
)
```

- `O1` is a flag of the Stan compiler (`stanc`). It rewrites the program
  so that vectors and matrices of parameters are stored as one block of
  values and one block of gradients instead of one pair per element,
  which is much friendlier to the automatic differentiation that
  dominates a Stan run. It is well tested, but it is the one option that
  changes the generated C++ nontrivially, which is why Stan has not yet
  made it the default.
- `STAN_CPP_OPTIMS` turns on extra C++ compiler optimizations when the
  program is built. It also rebuilds CmdStan’s own object files under
  those flags the first time it is used, so expect that first
  compilation to take several minutes longer than usual.
- `STAN_NO_RANGE_CHECKS` removes the bounds checks on every vector and
  array access. It is safe once a program is known to run without
  indexing errors, which the `cogmod` families are, and saves a little
  on every observation.

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/metric_dense_e.R)),
which fits the DDM, LBA, LNR and RDM to 900 lexical decision trials from
6 participants, with and without participant random intercepts, we found
that these flags bought nothing on our machine (Windows, the RTools
`gcc`). Gradient evaluations per second ranged from 20% below to 10%
above the default build across the 16 models, more often below it than
above, and compilation took up to 25 s longer. The flags act on the
compiled code only, so this should transfer across models but not
necessarily across compilers: they may still be worth a try on Linux or
macOS.

## Mass Matrix Adaptation (`metric = "dense_e"`)

The sampler’s **metric** (its mass matrix) is the next lever, and it too
leaves the posterior untouched: NUTS is exact under any metric, so the
only thing at stake is speed. During warmup, Stan estimates the
posterior’s scale so that its trajectories can take steps of the right
size in every direction. By default (`diag_e`) it learns **one variance
per parameter** and nothing about how parameters covary. That is the
wrong shape for these models: `boundary` and `ndt` often trade off
against each other (a later start with a lower threshold produces nearly
the same RTs as an earlier start with a higher one), as do drift and
boundary, and the two scale parameters of a race. The posterior is a
thin diagonal ridge, and a sampler that only knows the axis-aligned
widths has to creep along it with many small leapfrog steps per
iteration. `dense_e` estimates the **full covariance** instead, which
rotates and rescales that ridge into something close to a sphere, so
each iteration needs far fewer gradient evaluations and the draws are
less autocorrelated.

The argument goes straight through
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html) to
`cmdstanr` (with the `rstan` backend it is
`control = list(metric = "dense_e")`):

``` r

m <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  chains = 4, cores = 4,
  metric = "dense_e"
)
```

It is not free, and it does not always win. Two costs scale with the
number of parameters `d` that Stan samples. Each leapfrog step now
multiplies by a `d x d` matrix rather than a vector, and warmup has to
pin down `d(d+1)/2` covariances rather than `d` variances from the same
number of warmup iterations. For a handful of population-level
parameters, as in the [decision making
models](https://dominiquemakowski.github.io/cogmod/articles/decision_making.md),
both costs are negligible and the LBA there samples about twice as fast.
For a hierarchy with hundreds of participant-level parameters, the
covariance estimate gets noisy, the per-step cost is felt, and the
default is often the safer choice.

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/metric_dense_e.R)),
we found that this is exactly the case. Each model was sampled under
both metrics with the same data, priors, initial values, seed and 500
warmup iterations; the table gives the effective draws per second of
`dense_e` relative to the default, using the smallest bulk ESS across
all parameters:

| Family | Population-level only (8 to 11 parameters) | With participant random intercepts (22 to 32 parameters) |
|:---|:--:|:--:|
| DDM | 2.5 | 0.4 |
| LBA | 2.7 | 0.1 |
| LNR | 3.3 | 0.6 |
| RDM | 0.7 | 0.1 |

With population-level parameters only, `dense_e` cut the trajectories to
between a quarter (LBA, whose strongest posterior correlation was 0.91,
and RDM) and two thirds of their default length. For the DDM, LBA and
LNR that delivered two and a half to three times the effective draws per
second; for the RDM the shorter trajectories mixed less well per draw
and it came out slightly behind, on a single seed. With random
intercepts for 6 participants on two or three parameters, it lost in
every family and produced more divergences, even though the posterior
correlations were just as strong. What decided the outcome was the
number of parameters, not the strength of the correlations. So the rule
of thumb is: `dense_e` for population-level fits, the default as soon as
random effects enter, even a modest number of them.

## Approximating with Faster Algorithms

Everything so far produces the same posterior faster. The next two
sections trade something for speed. Before committing to a long MCMC
run, it is often worth getting a quick, approximate answer first.
[Pathfinder](https://mc-stan.org/docs/reference-manual/pathfinder.html)
(Zhang et al., [2022](https://arxiv.org/abs/2108.03782)) is a good first
choice and a significant improvement over other Variational Inference
(VI) algorithms like *meanfield ADVI* or *full-rank ADVI*. It is
implemented in `cmdstanr` and can be used in `brms` with the
`algorithm = "pathfinder"` argument.

``` r

m_pathfinder <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  algorithm = "pathfinder",
  chains = 16,
  single_path_draws = 4000,
  max_lbfgs_iters = 8000
)
```

Note that `chains` does not mean MCMC chains here, but controls the
number of Pathfinder **paths**, and more paths give the importance
resampling step more diverse material to draw from. Increasing
`single_path_draws` and `max_lbfgs_iters` similarly buys a better
approximation at very little extra cost. Also, consider using the
`threads` argument to parallelize the Pathfinder draws across multiple
cores.

Pathfinder draws are not a substitute for MCMC posteriors, they
approximate the posterior and are not guaranteed to be well calibrated,
especially for variance components and correlations. However, its
results can potentially be used to **tighten the priors**, which in turn
might help with MCMC convergence and sampling efficiency.

The other cheap approximation `brms` offers through `cmdstanr` is the
**Laplace approximation** (`algorithm = "laplace"`). It is the classical
recipe: find the posterior mode by optimization, take the curvature (the
Hessian) there, and treat the posterior as the multivariate Normal with
that mean and covariance. Everything happens on Stan’s unconstrained
scale, so a parameter on a log link comes out log-Normal rather than
Normal, which is usually the right shape. Pathfinder runs the same kind
of optimizer but does not wait for it to converge: along the
optimization path it builds a low-rank Normal approximation at every
iterate, keeps the one with the best evidence lower bound, and, with
several paths started from different points, pools their draws with
Pareto-smoothed importance resampling.

``` r

m_laplace <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  algorithm = "laplace",
  draws = 4000
)
```

Neither is uniformly better. Laplace uses the full Hessian, so for a
handful of parameters whose posterior is close to Normal on the
unconstrained scale (a non-hierarchical model with a few hundred trials)
it is typically the more accurate of the two, and it is the faster one.
It has three failure modes: the mode is not unique, the Hessian is not
positive definite (a flat or ridge-shaped direction, which the LBA’s
scale indeterminacy produces when nothing is fixed), or the posterior is
far from Normal, as variance components with few groups are.
Pathfinder’s covariance is a low-rank estimate and can be cruder for
small problems, but its multiple paths and importance resampling make it
more robust to these failures, and it reports a Pareto `k` diagnostic
that says when the approximation should not be trusted. It is also what
Stan recommends as a source of MCMC initial values. Both are cheap next
to MCMC (seconds for Laplace, a minute or so for Pathfinder on a mixed
model), so the practical answer is to run both: where they agree with
each other, they are probably both right, and where they disagree, MCMC
is the arbiter.

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/warm_start.R)),
on a DDM with participant random intercepts (6 participants, 900
trials), we found that Laplace took 3 s and Pathfinder 77 s against 4.5
min for MCMC, and that both got the population-level means about right
while neither could be trusted on spread. Pathfinder’s posterior SDs
were two to four times too small for the drift, boundary and `ndt`
intercepts and its smallest random-effect SD five times too large, and
its Pareto `k` of 1.5 said as much. Laplace was worse in the other
direction: the SDs of those same three intercepts were four to thirteen
times too large, and the three participant-level SDs, which MCMC put at
0.6, 0.16 and 0.05, came out at 4, 3.4 and 1.8 with uncertainties as
large as the estimates. The curvature of a hierarchical posterior at its
mode says little about its spread. For the condition effects, the
quantities one usually cares about, both approximations matched MCMC to
within a few hundredths.

## Warm Starts: Reusing Initial Values and Metrics

Warmup produces two things: a position inside the bulk of the posterior,
and an adapted metric and step size. Both are thrown away when the run
ends, and both can be supplied to the next run instead of being learned
again. Stan exposes them as the `init`, `inv_metric` and `step_size`
arguments of the sampler, and
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html) forwards any
argument it does not recognize to `cmdstanr`, so no custom Stan is
involved. The sampler keeps adapting from the supplied values during
whatever warmup remains, so the result is a head start rather than a
fixed setting, and nothing about it changes the posterior being sampled:
a poor warm start costs speed, not correctness. The usual checks
(`Rhat`, divergences, effective sample size) remain the judge of the
run.

### From a previous fit of the same model

The simplest source is a previous fit of the same model: a pilot run, a
fit that needs more draws, or the same model under a slightly different
prior or seed. `brms` keeps the adapted quantities of each chain in the
fit’s metadata (`attr(m$fit, "metadata")`), and
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
reads them back, averaged over the chains, together with the posterior
means as starting values. Passing them to a new run lets its warmup
shrink to what the step size needs:

``` r

ws <- cogmod_warmstart(m)   # inv_metric, step_size and init, from the fit itself

m_more <- brm(fit = m,      # reuse the compiled model
  init = ws$init,
  chains = 4, cores = 4,
  iter = 2100, warmup = 100,
  inv_metric = ws$inv_metric,
  step_size = ws$step_size
)
```

The metric must match the previous run’s shape: a `diag_e` fit stores a
vector of variances, a `dense_e` fit a full matrix, and
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
handles only the former.

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/warm_start.R)),
on a population-level LNR, we found that this halves the cost of a
refit:

| Run | Warmup | Wall time | Min. bulk ESS | ESS per second |
|:---|:--:|:--:|:--:|:--:|
| Reference | 500 | 17 s | 857 | 49 |
| Warm restart, stored metric and step size | 100 | 8 s | 788 | 102 |
| Cold start, same short warmup | 100 | 32 s | 391 | 12 |

### From a pilot fit on some of the participants

A pilot fit on the first few participants, while the model is still
being worked out, is often followed by the same model on the full
sample. That pilot has adapted a metric and a step size and has already
located the population-level parameters, so it is the natural warm start
for the full fit. Stan refuses its metric as is, though. Every
participant adds a standardized random effect (`z_1[1, j]`) to each
group-level term, so the full model has more parameters and the pilot’s
metric, one variance per parameter, has the wrong length. The remedy is
to map it across **by parameter name**: the population-level entries
carry over one to one, a pilot participant’s `z` entries move to that
participant’s slot in the full model, and the participants the pilot
never saw take the average of the pilot’s `z` variances for the same
term (they are standardized effects, so the average is a fair guess).
Initial values are built the same way from the pilot’s posterior means,
with the new participants starting at `z = 0`.
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
does all of this, in two steps. After the pilot, extract what its warmup
learned as a small table, labelled by Stan parameter name, and keep it:

``` r

warmstart <- as.data.frame(cogmod_warmstart(m_pilot))
write.csv(warmstart, "pilot_warmstart.csv", row.names = FALSE)
```

Then, for the full model, hand that table to the three helpers that
mirror
[`cogmod_priors()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_priors.md)
and
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md) -
the model’s formula and data first, the table under `warmstart` - one
per argument of
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html):

``` r

warmstart <- read.csv("pilot_warmstart.csv")

m <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  stanvars = cogmod_stanvars(f),
  init = cogmod_inits(f, df, warmstart = warmstart),
  inv_metric = cogmod_inv_metric(f, df, warmstart = warmstart),
  step_size = cogmod_step_size(f, df, warmstart = warmstart),
  backend = "cmdstanr", chains = 4, cores = 4,
  iter = 600, warmup = 100
)
```

Each helper maps the table onto the model it is given, which is where
the pilot’s participants are matched by name and the new ones filled in,
so the same table serves a refit of the pilot itself and the full sample
alike. It is a few kilobytes, which is how a pilot fitted on a laptop
can warm-start an array job on a cluster with no `brmsfit` in sight; the
pilot fit itself can also stand in for the table wherever `warmstart` is
taken. To see what the mapping did - how many entries came from the
pilot, how many are new participants, how many had no counterpart at
all - build the object for the target model and print it:
`cogmod_warmstart(warmstart, f, df)`, or straight from the fit
`cogmod_warmstart(m_pilot, data = df)`, since whatever it is not given
it takes from the pilot. The same works with a changed formula:
`cogmod_warmstart(m, formula = f2)` is a variant of the model on the
same data, whose shared parameters start where the first fit left them.
Anything the source never had, such as a predictor added to the formula,
gets Stan’s default variance and the generic starting value
[`cogmod_inits()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_inits.md)
would give it; many such entries mean the two formulas differ more than
intended.

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/warm_start_subset.R)),
on the mixed LNR (participant random intercepts on `mu` and `ndt`) and
the mixed DDM (on `mu`, `boundary` and `ndt`) of the previous sections,
a pilot on 4 of 8 participants warm-started the full fit to about twice
the effective draws per second of a cold start with the full warmup, and
four to six times those of a cold start with the same short warmup:

| Family | Run | Participants | Warmup | Wall time | Min. bulk ESS | ESS per second |
|:---|:---|:--:|:--:|:--:|:--:|:--:|
| LNR | Pilot fit | 4 | 500 | 78 s | 550 | 7.1 |
| LNR | Reference, cold start | 8 | 500 | 112 s | 589 | 5.2 |
| LNR | Cold start, short warmup | 8 | 100 | 247 s | 522 | 2.1 |
| LNR | Pilot initial values only | 8 | 100 | 265 s | 569 | 2.2 |
| LNR | [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md): initial values, metric and step size | 8 | 100 | 65 s | 603 | 9.3 |
| DDM | Pilot fit | 4 | 500 | 4.4 min | 337 | 1.29 |
| DDM | Reference, cold start | 8 | 500 | 11.2 min | 276 | 0.41 |
| DDM | Cold start, short warmup | 8 | 100 | 29.2 min | 246 | 0.14 |
| DDM | Pilot initial values only | 8 | 100 | 28.4 min | 298 | 0.18 |
| DDM | [`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md): initial values, metric and step size | 8 | 100 | 4.1 min | 209 | 0.85 |

Two things stand out. First, the initial values alone are worth nothing:
the runs that started from the pilot’s posterior means but let the
sampler find its own metric were as slow as the cold ones. A short
warmup does not hurt because the chains start in the wrong place, they
start close enough either way, but because the step size is left
unadapted, and the sampler then pays for it on every iteration of the
sampling phase, here with wall times two and a half times the
reference’s. Second, the pilot’s metric does not have to be right to
help. Its population-level variances were two to three times the full
fit’s (a posterior narrows as participants are added, and the LNR
pilot’s step size was accordingly smaller, 0.05 against 0.09), yet 100
iterations of adaptation from that start were enough. The
population-level estimates of every run agreed to two decimals. Note
that this is the diagonal metric; a `dense_e` pilot would need the same
treatment on a matrix, with zero covariances for the new participants,
which
[`cogmod_warmstart()`](https://dominiquemakowski.github.io/cogmod/reference/cogmod_warmstart.md)
does not do.

### From a Pathfinder approximation

The second source is an approximation from the previous section.
Pathfinder draws sit inside the bulk of the posterior, which makes them
good initial values, and their covariance is a first estimate of the
metric. Getting the covariance right requires the draws on Stan’s
**unconstrained** scale, which `cmdstanr` can do when the model is
compiled with its model methods, so this pipeline runs at the `cmdstanr`
level. `brms` still writes the program and the data, and the result is
wrapped back into a `brmsfit` at the end, so every post-processing
method works as usual.

``` r

library(cmdstanr)

prior <- cogmod_priors(f, df)
stanvars <- cogmod_stanvars(f)

# 1. brms writes the Stan program and the data; cmdstanr compiles it
scode <- stancode(f, data = df, prior = prior, stanvars = stanvars, backend = "cmdstanr")
sdata <- standata(f, data = df, prior = prior, stanvars = stanvars)
mod <- cmdstan_model(write_stan_file(scode), compile_model_methods = TRUE)

# 2. Pathfinder
pf <- mod$pathfinder(data = sdata, init = cogmod_inits(f, df), num_paths = 8, draws = 1000)

# 3. Its covariance on the unconstrained scale is the starting metric
inv_metric <- cov(pf$unconstrain_draws(format = "draws_matrix"))

# 4. MCMC starting from Pathfinder draws, with that metric, and a shorter warmup
fit <- mod$sample(
  data = sdata, init = pf,
  chains = 4, parallel_chains = 4,
  iter_warmup = 300, iter_sampling = 500,
  metric = "dense_e", inv_metric = inv_metric
)

# 5. Back into brms
m <- brm(f, data = df, prior = prior, stanvars = stanvars, backend = "cmdstanr", empty = TRUE)
m$fit <- read_csv_as_stanfit(fit$output_files(), variables = fit$metadata()$stan_variables, model = mod)
m <- rename_pars(m)
```

In our local benchmarking demo
([script](https://github.com/DominiqueMakowski/cogmod/blob/main/benchmarks/warm_start.R)),
on the same mixed DDM as above, we found that the pipeline works but did
not pay for itself at this size. Started from Pathfinder draws and
covariance and left to adapt for 300 iterations, the sampler finished
its own work in 3.7 min instead of the cold start’s 4.4 min with the
same effective sample size, but Pathfinder’s 77 s ate the difference,
for a net loss of about 10%. The arithmetic favours the warm start as
models grow, since Pathfinder scales gently with the number of
parameters while warmup for a poorly conditioned posterior does not, and
since a bad random initialization on a large model can waste far more
than the approximation costs. What one must **not** do is fix the metric
to the approximation (Stan allows it, by setting `init_buffer` to the
whole warmup so that no adaptation window runs). Pathfinder’s covariance
was, as the previous section showed, several times too small on most
parameters, and a sampler that believes it takes tiny steps and hits the
maximum treedepth: that run took 15 min, diverged on a third of its
transitions, and returned an effective sample size of 10.

## High-Performance Clusters (HPCs)

Chain and thread parallelization is limited by the cores on one machine.
The next step is to spread chains across many machines - typically many
short jobs on an HPC scheduler rather than one long job. Because warmup
dominates total sampling time for these models, recruiting many nodes
for a **short** run (fewer post-warmup iterations each) is usually more
efficient than recruiting a few nodes for a long one: each job pays the
fixed compilation and warmup cost once, but the post-warmup work is what
actually needs to scale with the number of draws you want.

A typical job array (SLURM shown here) runs one or two chains per node,
threading each chain across the cores left over after dividing them
among the node’s chains:

``` r

task_id <- as.numeric(Sys.getenv("SLURM_ARRAY_TASK_ID", unset = "1"))
total_cores <- as.numeric(Sys.getenv("SLURM_CPUS_PER_TASK", unset = "2"))

chains_per_node <- 2
threads_per_chain <- total_cores / chains_per_node  # e.g. 16 cores / 2 chains = 8 threads per chain

warmup <- 1000
iter <- warmup + 200  # short post-warmup draw per node; combined across nodes below

m_node <- brm(f,
  data = df,
  prior = cogmod_priors(f, df),
  init = cogmod_inits(f, df),
  stanvars = cogmod_stanvars(f),
  backend = "cmdstanr",
  chains = chains_per_node,
  cores = chains_per_node,
  threads = threading(threads_per_chain),
  warmup = warmup,
  iter = iter,
  seed = task_id,
  file = paste0("m_node_", task_id, ".rds")
)
```

Each array task produces its own `brmsfit` with a handful of chains.
Once all tasks have finished, combine them into a single fit with
[`brms::combine_models()`](https://paulbuerkner.com/brms/reference/combine_models.html),
which concatenates the post-warmup draws across fits (they must share
the same model structure and data):

``` r

fits <- lapply(list.files(pattern = "^m_node_.*\\.rds$"), readRDS)
m <- do.call(brms::combine_models, fits)
```

The result behaves like a single
[`brm()`](https://paulbuerkner.com/brms/reference/brm.html) fit with
`chains_per_node * n_nodes` chains worth of draws, obtained in roughly
the time a single node’s chains would have taken.

## Future Directions

A different line of work sidesteps the repeated sampling of MCMC/VI
entirely: **amortized inference**, where a neural network is trained
(once, offline, potentially at significant upfront cost) to map observed
data directly to an approximate posterior, so that inference on new
datasets afterwards is close to instantaneous rather than requiring a
fresh MCMC run each time. [BayesFlow](https://bayesflow.org/) is the
most actively developed toolkit in this space, and has already been
applied to evidence accumulation and other cognitive models.
