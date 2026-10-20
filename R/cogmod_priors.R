#' Priors that make a cogmod posterior proper
#'
#' @description
#' Fills in weakly informative priors for the parameters `brms` would otherwise
#' leave flat, or give a default meant for another parameterization, for the
#' model you are actually fitting.
#'
#' **For the families with a non-decision time this is not a convenience.**
#' `brms` gives the intercept of a custom-family parameter it does not
#' recognise a flat, improper prior, and on the `ndt` + `poutlier` mixture the
#' likelihood is exactly flat in two directions: `poutlier` toward 1, where
#' every response belongs to the outlier component and the decision parameters
#' drop out of the density, and `ndt` toward 0, where the gradient with respect
#' to `log(ndt)` vanishes. Flat prior plus infinite flat region is an improper
#' posterior, and the fit does not fail loudly: it returns intercepts around
#' `1e14` with `Rhat` near 2 and an effective sample size of about 5. Several
#' families have flat directions of their own (see *Why these parameters*).
#'
#' @details
#' The function starts from [brms::get_prior()] for the model in hand, edits the
#' rows it has a better prior for, and returns the result of
#' [brms::validate_prior()]: a `brmsprior` object, printed as a table that
#' includes the defaults `brms` will use. Because every row comes from the model
#' `brms` will build, no prior can miss its parameter - `0 + Intercept`,
#' interactions, group-level terms and smooths need no special handling - and a
#' malformed specification errors here rather than inside `brm()`.
#'
#' # Setting your own priors
#'
#' Combine the result with [brms::prior()] entries. `replace = TRUE` replaces
#' one of these defaults; without it the entry is added for another slot:
#'
#' ```r
#' priors <- c(
#'   cogmod_priors(f, df),
#'   brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
#'   replace = TRUE
#' )
#' ```
#'
#' # Supported families
#'
#' The family is read off `formula`, so build it with
#' `brms::bf(..., family = cogmod_lognormal())`. Covered: every family with an
#' `ndt` + `poutlier` component - [cogmod_lognormal()], [cogmod_logstudent()],
#' [cogmod_loggamma()], [cogmod_invgaussian()], [cogmod_exwald()],
#' [cogmod_bisa()], [cogmod_gamma()], [cogmod_invgamma()], [cogmod_weibull()],
#' [cogmod_invweibull()], [cogmod_logweibull()], [cogmod_lba1()] and the
#' choice-and-RT families [cogmod_lnr()], [cogmod_rdm()], [cogmod_lba2()] and
#' [cogmod_ddm()] - plus [cogmod_exgaussian()], [cogmod_geg()] and the
#' rating-scale families [cogmod_choco()], [cogmod_betagate()] and
#' [cogmod_betadiscrete()]. Any other family, or a formula carrying none, gets a
#' message and `brms`'s defaults unchanged, so the call is safe to leave in a
#' script.
#'
#' # What gets set
#'
#' For each parameter in the table under *The priors*:
#'
#' - its **intercept** gets the prior in the table;
#' - its **slopes** get `normal(0, 0.2)` on the link scale, or the wider prior
#'   in the table. On a log or logit link a flat slope prior is not harmless: it
#'   can carry a condition into the flat regions while the intercept stays well
#'   behaved;
#' - its **group-level SDs** get `exponential(1)`, for the same reason, and the
#'   correlations keep `brms`'s LKJ;
#' - its **smooth terms** follow *Smooth terms*, as every smooth does.
#'
#' Parameters not in the table, such as [cogmod_lognormal()]'s `sigma`, keep
#' `brms`'s defaults, and so does the response's own `mu`: its slopes are the
#' effects being estimated, and `brms`'s `student_t()` default is a fair
#' intercept prior on most links. Two exceptions. [cogmod_exgaussian()] and
#' [cogmod_geg()] give `mu` an intercept prior, because there it is a time in
#' seconds on an `identity` link. And in the race families, [cogmod_lnr()],
#' [cogmod_rdm()] and [cogmod_lba2()], `mu` is the rate of one accumulator, so
#' its group-level SD gets the `exponential(1)` the other accumulator's gets
#' and swapping the `dec()` coding does not change what the priors say.
#'
#' A row `brms` left flat is filled. A row `brms` already filled is replaced
#' only where its default was written for a parameter of the same name in one
#' of its own families, and is wrong here: `uniform(0, min_Y)` on an `ndt` left
#' out of `bf()`, the min-RT bound the outlier component exists to remove;
#' `gamma(0.01, 0.01)` on [cogmod_loggamma()]'s `shape`, which must be free to
#' go negative; and `student_t(3, 0, 2.5)` and the like on the `sigma` of
#' [cogmod_exgaussian()] and [cogmod_geg()], the `shape` of [cogmod_geg()] and
#' the `phi` of [cogmod_betagate()] and [cogmod_betadiscrete()] - a Gaussian SD
#' wider than most whole RT distributions, or on [cogmod_betadiscrete()]'s
#' `log` link a precision reaching 2853.
#'
#' # ndt and poutlier
#'
#' `ndt` is in **seconds**, like everything else in these families.
#' `normal(-1.2, 0.5)` on its `log` link centres it on 0.30 s with 95% of the
#' mass between 0.11 and 0.80 s: wide enough for slow populations and demanding
#' responses, and still a fence against `ndt -> 0`.
#'
#' `poutlier` is a proportion. `normal(-6, 1)` on its `logit` link centres it
#' on 0.25%, with 95% between 0.035% and 1.7%, which is where the empirical
#' estimates sit: Miller (2024) puts the outlier proportion below 0.5% across
#' four lexical-decision megastudies, and the untrimmed data of the
#' [outliers article](https://dominiquemakowski.github.io/cogmod/articles/outliers.html)
#' agree. If the data were trimmed before fitting, tighten the prior
#' (`normal(-7, 0.5)`) rather than fixing `poutlier = 0`, which brings the hard
#' min-RT boundary back; see the *Trimmed data* section of [cogmod_lognormal()].
#'
#' # The priors
#'
#' On the link scale for a parameter in `bf()` - for its intercept, or under
#' `0 + Intercept` its coefficient named `Intercept` - and on the natural scale
#' for one left out (see *Parameters left out of the formula*).
#'
#' | Parameter | Families | In `bf()` | Left out | Slopes |
#' | --- | --- | --- | --- | --- |
#' | `ndt` | all with `ndt` | `normal(-1.2, 0.5)` | `lognormal(-1.2, 0.5)` | `normal(0, 0.2)` |
#' | `poutlier` | all with `ndt` | `normal(-6, 1)` | `exponential(280)` | `normal(0, 0.2)` |
#' | `sigmabias` | [cogmod_lognormal()], [cogmod_lnr()] | `normal(0, 1)` | `lognormal(-0.35, 0.75)` | `normal(0, 0.5)` |
#' | `dof` | [cogmod_logstudent()] | `normal(1.8, 0.7)` | `lognormal(1.8, 0.7)` | `normal(0, 0.3)` |
#' | `shape` | [cogmod_loggamma()] | `normal(0, 0.5)` | `normal(0, 0.5)` | `normal(0, 0.2)` |
#' | `boundary` | [cogmod_invgaussian()], [cogmod_exwald()], [cogmod_bisa()] | `normal(0.5, 1)` | `lognormal(0, 0.6)` | `normal(0, 0.5)` |
#' | `sigmadrift` | [cogmod_invgaussian()] | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
#' | `sigmandt` | [cogmod_invgaussian()], [cogmod_ddm()] | `normal(-3, 1)` | `lognormal(-3, 1)` | `normal(0, 0.5)` |
#' | `tau` | [cogmod_exwald()], [cogmod_exgaussian()], [cogmod_geg()] | `normal(-1.5, 0.7)` | `lognormal(-1.5, 0.7)` | `normal(0, 0.5)` |
#' | `sigmabias`, `boundary` | [cogmod_lba1()], [cogmod_rdm()], [cogmod_lba2()] | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
#' | `sigmazero`, `sigmaone` | [cogmod_lnr()], [cogmod_lba2()] | `normal(0, 1)` | `lognormal(-0.7, 0.75)` | `normal(0, 0.5)` |
#' | `nuone` | [cogmod_lnr()] | `normal(0.7, 1.5)` | `normal(0.7, 1.5)` | `normal(0, 0.5)` |
#' | `driftone` | [cogmod_rdm()] | `normal(3, 2)` | `lognormal(1, 0.75)` | `normal(0, 2)` |
#' | `driftone` | [cogmod_lba2()] | `normal(1, 2)` | `normal(1, 2)` | `normal(0, 1.5)` |
#' | `boundary` | [cogmod_ddm()] | `normal(1, 1)` | `lognormal(0.27, 0.5)` | `normal(0, 0.5)` |
#' | `bias` | [cogmod_ddm()] | `logistic(0, 1)` | `brms`'s `beta(1, 1)` | `normal(0, 0.5)` |
#' | `sigmadrift` | [cogmod_ddm()] | `normal(0, 1)` | `lognormal(-1, 0.75)` | `normal(0, 0.5)` |
#' | `sigmabias` | [cogmod_ddm()] | `normal(-2, 1)` | `beta(1, 5)` | `normal(0, 0.5)` |
#' | `mu` | [cogmod_exgaussian()], [cogmod_geg()] | `normal(0.4, 0.25)` | - | `brms`'s |
#' | `sigma` | [cogmod_exgaussian()], [cogmod_geg()] | `normal(-2.3, 0.7)` | `lognormal(-2.3, 0.7)` | `normal(0, 0.5)` |
#' | `shape` | [cogmod_geg()] | `normal(0, 0.5)` | `lognormal(0, 0.5)` | `normal(0, 0.5)` |
#' | `confright`, `confleft` | [cogmod_choco()] | `normal(0, 1)` | `beta(2, 2)` | `normal(0, 0.5)` |
#' | `precright`, `precleft` | [cogmod_choco()] | `normal(2, 1.5)` | `lognormal(0.7, 0.7)` | `normal(0, 0.5)` |
#' | `phi` | [cogmod_betagate()] | `normal(2, 1.5)` | `lognormal(0.7, 0.7)` | `normal(0, 0.5)` |
#' | `phi` | [cogmod_betadiscrete()] | `normal(0.7, 0.8)` | `lognormal(0.7, 0.8)` | `normal(0, 0.5)` |
#' | `pex` | [cogmod_choco()], [cogmod_betagate()] | `normal(-2, 1)` | `beta(2, 12)` | `normal(0, 0.5)` |
#' | `bex` | [cogmod_choco()], [cogmod_betagate()] | `normal(0, 1)` | `beta(2, 2)` | `normal(0, 0.5)` |
#' | `pmid` | [cogmod_choco()] | `normal(-2.5, 1)` | `exponential(9)` | `normal(0, 0.5)` |
#' | `pzero` | [cogmod_betadiscrete()] | `normal(-2.5, 1)` | `exponential(9)` | `normal(0, 0.5)` |
#'
#' # Parameters left out of the formula
#'
#' Writing `ndt ~ 1` and leaving `ndt` out are not the same model to `brms`. A
#' parameter in `bf()`, even as `~ 1`, gets a linear predictor and lives on its
#' **link** scale, reported as `ndt_Intercept` under *Regression
#' Coefficients*. One left out is declared as a plain auxiliary parameter, like
#' `sigma` in `gaussian()`, and lives on its **natural** scale, under *Further
#' Distributional Parameters*, so it needs a prior on that scale. The two
#' columns of the table describe the same belief wherever they can: exactly on
#' a `log` or `identity` link, where `lognormal` is `normal` on the log scale,
#' and matched on the median elsewhere.
#'
#' The exception is a probability whose omission is itself a statement. Leaving
#' `poutlier` out says the data were trimmed or are not expected to contain
#' outliers, and leaving `pmid` or `pzero` out says the event is not expected.
#' Those get an exponential, with its **mode at zero**, which no logit-scale
#' prior can have, and the same median: `exponential(280)` against
#' `plogis(-6) = 0.0025`, and `exponential(9)` against `plogis(-2.5) = 0.076`.
#' It is still a prior, so a genuine spike of fast responses still pulls the
#' rate up.
#'
#' # Why these parameters
#'
#' Each family's own page has the full account. In short, a parameter is in the
#' table for one of five reasons:
#'
#' - **Its link puts a limit at infinity, and the likelihood is flat long
#'   before it** - the `ndt -> 0` failure again. That is the start-point ranges
#'   (`sigmabias`) and the LBA and RDM thresholds that share a ridge with them
#'   (`boundary`), the between-trial variabilities (`sigmadrift`, `sigmandt`),
#'   and [cogmod_logstudent()]'s `dof`, which reaches the LogNormal only as
#'   `dof -> Inf`. Where the design cannot identify a range or a variability,
#'   fixing it in `bf()` (`sigmabias = 0`) is better than estimating it behind
#'   a prior.
#' - **An accumulator that stops winning.** Push a race accumulator's rate down
#'   far enough and the likelihood stops depending on it, or on its SD: `nuone`
#'   and the SDs in [cogmod_lnr()], `driftone` in [cogmod_rdm()] and
#'   [cogmod_lba2()]. Under `dec(Error)` that is the error accumulator, and
#'   wherever errors are rare its posterior sits in the lower tail of the prior:
#'   the fence working, not a prior at odds with the data. `mu` has the mirror
#'   plateau; if the rarely chosen option is on `mu`, copy the prior onto it.
#' - **An arbitrary scale.** [cogmod_lba2()]'s likelihood is exactly constant
#'   along a ray that rescales every evidence parameter. The priors on
#'   `sigmazero` and `sigmaone` make the posterior proper; only fixing one of
#'   them in `bf()` (`sigmazero = 1`) identifies the scale.
#' - **A ridge with another parameter**: `tau` with `ndt` in
#'   [cogmod_exwald()], `shape` with `mu` in [cogmod_geg()], `dof` with
#'   `poutlier` at the fast end of [cogmod_logstudent()]; and
#'   [cogmod_loggamma()]'s `shape`, whose prior also keeps the sampler clear of
#'   `sigma * shape >= 1`, where the density is unbounded at the shift.
#' - **Identified, but `brms` leaves it flat or gives it a default for another
#'   parameterization**: the Wald and DDM thresholds, the DDM's `bias`, and the
#'   ex-Gaussian and rating-scale parameters (below).
#'
#' No prior is set on the shape of [cogmod_weibull()] or [cogmod_gamma()],
#' although a shape below 2 makes their `ndt` gradient unbounded: the
#' likelihood itself prefers that region, by around 100 log units on the data
#' in the [RT models article](https://dominiquemakowski.github.io/cogmod/articles/rt_models.html),
#' so a prior weak enough to be a default could only bias it.
#' `?rcogmod_weibull` sets out what to do instead.
#'
#' # The ex-Gaussian and the GEG
#'
#' [cogmod_exgaussian()] has no `ndt`, but `mu`, `sigma` and `tau` are all
#' lengths of time in seconds, which `brms` cannot know. On the `softplus`
#' link, `sigma` (the Gaussian SD) gets 25 to 330 ms with a median of 96 ms,
#' and `tau` (the mean of the exponential tail) 55 to 630 ms with a median of
#' 201 ms: the range they occupy across the usual simple- and choice-RT tasks.
#' `mu` is on an `identity` link, where `brms`'s default, a `student_t()` with
#' a scale of 2.5 s, spreads the Gaussian's centre over several seconds either
#' way, and under `0 + Intercept` it gets no prior at all. `normal(0.4, 0.25)`
#' puts 95% of it between -0.09 and 0.89 s and deliberately allows negative
#' values: for fast, heavily tailed data the Gaussian component belongs near
#' or below zero, with `tau` carrying the mass. It is the centre of the
#' Gaussian alone, so the mean of the distribution is `mu + tau`.
#' [cogmod_geg()] takes the same rows, and a `shape` centred on the
#' ex-Gaussian it nests (`shape = 1`, with 95% between 0.37 and 2.7): `shape`
#' and `mu` sit on a ridge, and the prior holds it.
#'
#' # The rating-scale families
#'
#' [cogmod_choco()], [cogmod_betagate()] and [cogmod_betadiscrete()] have no
#' infinite flat region of the RT kind, but every parameter arrives either flat
#' or with a `brms` default meant for another parameterization.
#'
#' - **Point masses.** `pmid` is the probability of exactly the midpoint, and
#'   `pzero` of an extra category outside the scale. When the event is absent
#'   from the data the likelihood rises all the way to zero, and under a flat
#'   prior the logit runs off to minus infinity: `poutlier`'s failure, and
#'   `poutlier`'s treatment. `normal(-2.5, 1)` centres it on 8%, with 95% below
#'   37%. To switch one off, fix it in `bf()` (`pmid = 0`); a response exactly
#'   on the midpoint is then impossible, and the fit will not start.
#' - **Extremes.** `pex`, the total probability of a 0 or a 1, is *not* given a
#'   mode at zero: extreme responding is a response style, not a contaminant.
#'   `normal(-2, 1)` centres it on 12%, with 95% between 2% and 49%. `bex`,
#'   which end the extremes favour, is centred on symmetric, which also keeps
#'   it away from 0 and 1, where one gate closes and every response at that end
#'   has zero probability.
#' - **Beta shapes.** `confright` and `confleft`, the Beta means of each side,
#'   are centred on the middle of their half. The precisions are the awkward
#'   ones: below about 1 the Beta is U-shaped and unbounded at both ends, and a
#'   large one concentrates it sharply. `normal(2, 1.5)` on the `softplus` link
#'   spans 0.33 to 4.95, leaving 17% of its mass below 1, since a polarising
#'   item really can produce a U. [cogmod_betadiscrete()]'s `phi` is on a `log`
#'   link: `normal(0.7, 0.8)`, 0.42 to 9.7 with a median of 2.
#'
#' `mu` is left to `brms`, whose `student_t(3, 0, 2.5)` is the standard choice
#' on a `logit` link.
#'
#' # Smooth terms
#'
#' A smooth's SD, `sds`, is not in link units. `brms` writes the penalised part
#' of an `s()` or `t2()` term as `sds * Zs %*% zs` with `zs ~ N(0, 1)`, and how
#' far one unit of `sds` moves the linear predictor depends on the basis `mgcv`
#' builds: about 0.45 link units for `s(x)`, 6 for `s(x, bs = "cr")` and 0.1
#' for `t2(x, z, k = c(5, 5), bs = c("cr", "cr"))`, almost regardless of the
#' data. So each smooth term gets `exponential(rate)`, with the rate set to
#' that term's own scale (read off [brms::make_standata()]): `exponential(1)`
#' on the wiggle measured in link units, median 0.69 and 95% below 3, whatever
#' the basis. A term whose scale cannot be computed gets `exponential(1)`. This
#' applies to every smooth, the response's included, and it matters however
#' many rows there are: a smooth's SD is informed by the handful of
#' coefficients behind each penalty, as a group-level SD is by the number of
#' groups.
#'
#' The unpenalised part of a smooth - the linear trend of an `s(x)`, the
#' marginal trends of a `t2()` - is reported under class `b`, where it would
#' take its parameter's slope prior; but `mgcv` sets the scale of those
#' columns, not the user. On the parameters in the table, each such column
#' gets `normal(0, s / rms)` instead, which is the slope prior `normal(0, s)`
#' stated on the trend at a typical row, `rms` being the column's root mean
#' square.
#'
#' # Centring the priors on a previous fit
#'
#' `warmstart` takes what [cogmod_warmstart()] takes - a `brmsfit`, the table
#' [as.data.frame()] makes of one, or the path to a CSV of it - and replaces
#' each prior with `normal(median, prior_scale * sd)`, from that fit's
#' posterior for the same parameter. It applies to the population-level
#' intercepts and coefficients, the group-level SDs and the parameters left out
#' of the formula; correlations keep their LKJ. A parameter the source never
#' saw keeps the prior it would get anyway, and [cogmod_warmstart()]'s
#' `print()` says how many there are.
#'
#' **This changes the posterior.** Unlike [cogmod_inv_metric()] and
#' [cogmod_step_size()], which only change how the sampler moves, a prior is
#' part of the model. **If the source was fitted to data the new model also
#' contains, it counts those data twice**: a pilot on 5 of 10 participants used
#' to centre the priors for all 10 uses those 5 twice, and the intervals come
#' out too narrow by an amount nothing in the output shows. Use it when the
#' source is an independent data set - last year's sample, another lab's,
#' another session - or for a deliberately sequential analysis.
#'
#' `prior_scale = 1` makes the source's posterior the prior. The default of 3
#' has nine times its variance: enough to keep the sampler in the right region
#' and out of the flat directions, weak enough for data that disagree to win.
#'
#' ```r
#' # last year's sample as the prior for this year's
#' priors <- cogmod_priors(f, df, warmstart = fit_2025)
#' priors <- cogmod_priors(f, df, warmstart = fit_2025, prior_scale = 6)  # weaker
#' ```
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = cogmod_lognormal())`.
#' @param data The data, as passed to `brms::brm()`.
#' @param ... Passed to `brms::get_prior()` and `brms::validate_prior()`, for
#'   arguments such as `data2` or `knots`.
#' @param warmstart Optional. A previous fit to centre the priors on: a
#'   `brmsfit`, a `cogmod_warmstart` object, its data frame, or the path to a
#'   CSV of it, as [cogmod_warmstart()] takes. `NULL` (the default) leaves the
#'   priors where this function put them. **Changes the posterior** - see
#'   *Centring the priors on a previous fit*.
#' @param prior_scale The prior SD, as a multiple of the source's posterior
#'   SD. Only used with `warmstart`. Larger is weaker; 1 would use the
#'   source's posterior as the prior.
#'
#' @return A `brmsprior` object, to pass to `brms::brm(prior = )`.
#'
#' @references
#' - Miller, J. (2024). Estimating the proportions and latencies of reaction
#'     time outliers: A pooling method and case study of lexical decision
#'     tasks. *Behavior Research Methods*, *56*(7), 7280-7306.
#'     \doi{10.3758/s13428-024-02419-y}
#'
#' @seealso [cogmod_inits()], [cogmod_stanvars()], [cogmod_warmstart()]
#'
#' @examples
#' d <- data.frame(RT = rcogmod_lognormal(50, ndt = 0.3, poutlier = 0.02))
#' f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
#'   family = cogmod_lognormal()
#' )
#' cogmod_priors(f, d)
#'
#' # Replace a default, or append a prior for another parameter.
#' priors <- c(
#'   cogmod_priors(f, d),
#'   brms::prior(normal(-2, 0.1), class = "Intercept", dpar = "ndt"),
#'   replace = TRUE
#' )
#'
#' @export
cogmod_priors <- function(formula, data, ..., warmstart = NULL, prior_scale = 3) {
  family <- .cogmod_family(formula)
  fam <- .family_name(family)

  # First, whether the response suits the family at all (milliseconds, say);
  # the header of cogmod_checkdata.R says why the check lives here.
  .cogmod_checkdata(formula, data, family)

  if (isTRUE(fam %in% .OUTLIER_FAMILIES)) {
    out <- .priors_shifted(formula, data, family, ...)
  } else if (isTRUE(fam %in% names(.PRIORS_PLAIN))) {
    spec <- .PRIORS_PLAIN[[fam]]
    out <- .priors_dpars(formula, data, family, spec$prior, spec$override, ...)
  } else {
    message(
      "cogmod_priors() has nothing to add for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'; returning the brms defaults unchanged. ",
      "Currently supported: ", paste(.priors_families(), collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_lognormal())."
    )
    out <- brms::empty_prior()
  }

  # Last, so that it can override any of the above: a previous fit's posterior
  # is a stronger statement about where a parameter is than a default location
  # chosen for the family in general.
  if (!is.null(warmstart)) {
    out <- .priors_warmstart(out, formula, data, family, warmstart, prior_scale, ...)
  }

  args <- list(out, formula, data, ...)
  if (!is.null(family)) args$family <- family
  do.call(brms::validate_prior, args)
}


# Re-centres the priors on a previous fit's posterior: normal(median, scale *
# sd) for every prior row whose parameter the source has a posterior for. The
# join is on the address a get_prior() row carries - class, dpar, coef, group -
# which cogmod_warmstart() attaches to each unconstrained parameter, so the
# prior is about the very parameter the source sampled and no transformation
# is involved (`Intercept` is the centred intercept in both). Rows of `base`
# for anything the source does not know survive.
#' @keywords internal
.priors_warmstart <- function(base, formula, data, family, warmstart, prior_scale, ...) {
  if (!is.numeric(prior_scale) || length(prior_scale) != 1L ||
      !is.finite(prior_scale) || prior_scale <= 0) {
    stop("`prior_scale` must be a single positive number: the multiple of the ",
         "source's posterior SD to use as the prior SD.", call. = FALSE)
  }
  post <- cogmod_warmstart(warmstart, formula = formula, data = data,
                           jitter = 0, ...)$prior
  if (!nrow(post)) {
    message("The warm-start source carries no posterior median and SD, so there ",
            "is nothing to centre the priors on; leaving them as they are. A ",
            "table written before cogmod_priors(warmstart = ) existed has only ",
            "the metric and the means - rewrite it from the fit.")
    return(base)
  }

  args <- list(formula, data = data, ...)
  if (!is.null(family)) args$family <- family
  p <- do.call(brms::get_prior, args)
  # Kept whole: the only place the full coefficient list still exists, which
  # is what says whether a blanket row has been left covering nothing.
  all_rows <- p

  key <- function(d) paste(d$class, d$dpar, d$coef, d$group, sep = "|")
  idx <- match(key(p), key(post))
  hit <- !is.na(idx) & !nzchar(p$resp)
  if (!any(hit)) {
    message("None of the model's prior slots matched a parameter in the warm-start ",
            "source, so the priors are unchanged. Usually this means the source ",
            "was written for a different formula.")
    return(base)
  }

  new <- p[hit, , drop = FALSE]
  new$prior <- sprintf("normal(%s, %s)",
                       .num_str(post$median[idx[hit]]),
                       .num_str(prior_scale * post$sd[idx[hit]]))

  out <- rbind(base[!key(base) %in% key(new), , drop = FALSE], new)
  # Setting every coefficient a blanket row covers leaves the blanket covering
  # nothing, which brms warns about; the same check .priors_dpars() ends with.
  gone <- vapply(seq_len(nrow(out)), .blanket_now_unused, logical(1),
                 p = out, all_rows = all_rows)
  out[!gone, , drop = FALSE]
}


# A number as Stan will read it back. "g" keeps eight significant digits and
# falls to exponent notation only when that is shorter, which Stan accepts.
#' @keywords internal
.num_str <- function(x) trimws(formatC(x, digits = 8, format = "g"))


# Every family cogmod_priors() edits, for the message above.
#' @keywords internal
.priors_families <- function() unique(c(.OUTLIER_FAMILIES, names(.PRIORS_PLAIN)))


# Families with no ndt + poutlier component whose dpars brms still leaves flat,
# or gives a default meant for another parameterization. Same shape as the
# registries' `prior` slot: per dpar, `link` is the intercept prior when the
# dpar is in bf() (link scale), `nat` the prior when it is left out (brms then
# declares it on the natural scale), `slope` the prior on its slopes.
# `override` names the dpars whose modelled intercept is replaced even though
# brms filled it; a left-out dpar's row is always replaced, so needs no
# listing. ?cogmod_priors sets out the reasoning for users; what follows is
# what a maintainer needs on top.
#' @keywords internal
.PRIORS_PLAIN <- list(
  # `tau` arrives flat. `sigma` arrives with student_t(3, 0, 2.5), because
  # brms knows the name, centred on softplus(0) = 0.69 s, hence `override`.
  # `mu` is the response's own intercept, which only a family naming it here
  # touches (see .priors_dpars()); on this identity link brms's default has a
  # scale of 2.5 s, and under `0 + Intercept` there is none. The `nat`
  # lognormals carry the link priors' numbers because softplus(x) ~ exp(x)
  # below x = -2: median 0.100 against 0.096 for `sigma`.
  cogmod_exgaussian = list(
    prior = list(
      mu = c(link = "normal(0.4, 0.25)"),
      sigma = c(link = "normal(-2.3, 0.7)", nat = "lognormal(-2.3, 0.7)",
                slope = "normal(0, 0.5)"),
      tau = c(link = "normal(-1.5, 0.7)", nat = "lognormal(-1.5, 0.7)",
              slope = "normal(0, 0.5)")
    ),
    override = "sigma"
  ),

  # The ex-Gaussian's rows, plus `shape`, which sits on a ridge with `mu`: they
  # correlate about -0.98 at the maximum likelihood on real RTs, because
  # raising the CDF to a power shifts the distribution as well as bending it.
  # On the log link 0 is shape = 1, the ex-Gaussian, and normal(0, 0.5) keeps
  # 95% of the mass in 0.37-2.7, wide enough for the shapes the family exists
  # to reach. brms fills both forms of `shape`, hence `override`.
  cogmod_geg = list(
    prior = list(
      mu = c(link = "normal(0.4, 0.25)"),
      sigma = c(link = "normal(-2.3, 0.7)", nat = "lognormal(-2.3, 0.7)",
                slope = "normal(0, 0.5)"),
      tau = c(link = "normal(-1.5, 0.7)", nat = "lognormal(-1.5, 0.7)",
              slope = "normal(0, 0.5)"),
      shape = c(link = "normal(0, 0.5)", nat = "lognormal(0, 0.5)",
                slope = "normal(0, 0.5)")
    ),
    override = c("sigma", "shape")
  ),

  # The rating-scale families. The point masses `pmid` and `pzero` get
  # poutlier's treatment, a mode at zero when left out: exponential(9) matches
  # plogis(-2.5) = 0.076 on the median and puts 1.2e-4 of its mass above 1,
  # which the declared bound truncates. `pex` deliberately does not, and
  # normal(0, 1) on `bex` also fences bex = 0 or 1, where a gate closes and a
  # response at that end has zero probability. The precisions' `nat` form
  # matches normal(2, 1.5) on the median only (2.01 against 2.13): near x = 2
  # softplus(x) is about x rather than exp(x), so no lognormal is the same
  # distribution. The slopes get normal(0, 0.5) rather than the 0.2 default,
  # an odds ratio of 0.38-2.7 per unit, because these are the effects the
  # models are fitted for. cogmod_choco() needs no `override`: brms knows none
  # of its names, so its precisions arrive flat where `phi` arrives filled.
  cogmod_choco = list(
    prior = list(
      confright = c(link = "normal(0, 1)", nat = "beta(2, 2)",
                    slope = "normal(0, 0.5)"),
      confleft = c(link = "normal(0, 1)", nat = "beta(2, 2)",
                   slope = "normal(0, 0.5)"),
      precright = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
                    slope = "normal(0, 0.5)"),
      precleft = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
                   slope = "normal(0, 0.5)"),
      pex = c(link = "normal(-2, 1)", nat = "beta(2, 12)",
              slope = "normal(0, 0.5)"),
      bex = c(link = "normal(0, 1)", nat = "beta(2, 2)",
              slope = "normal(0, 0.5)"),
      pmid = c(link = "normal(-2.5, 1)", nat = "exponential(9)",
               slope = "normal(0, 0.5)")
    )
  ),

  # `phi` is overridden: brms knows the name from its beta family and supplies
  # student_t(3, 0, 2.5), loose on this softplus link and close to improper on
  # cogmod_betadiscrete()'s log link (95% up to phi = 2853).
  cogmod_betagate = list(
    prior = list(
      phi = c(link = "normal(2, 1.5)", nat = "lognormal(0.7, 0.7)",
              slope = "normal(0, 0.5)"),
      pex = c(link = "normal(-2, 1)", nat = "beta(2, 12)",
              slope = "normal(0, 0.5)"),
      bex = c(link = "normal(0, 1)", nat = "beta(2, 2)",
              slope = "normal(0, 0.5)")
    ),
    override = "phi"
  ),

  cogmod_betadiscrete = list(
    prior = list(
      # Log link, so both forms are one distribution: 0.42 to 9.7, median 2,
      # 19% below 1, where the U and J shapes this family exists for live.
      # Past 10 a seven-point scale has 95% of its mass in three categories.
      phi = c(link = "normal(0.7, 0.8)", nat = "lognormal(0.7, 0.8)",
              slope = "normal(0, 0.5)"),
      pzero = c(link = "normal(-2.5, 1)", nat = "exponential(9)",
                slope = "normal(0, 0.5)")
    ),
    override = "phi"
  )
)


# The two dpars the ndt + poutlier parameterization introduces, and
# cogmod_loggamma()'s `shape`, in the shape of a registry's `prior` slot;
# .priors_shifted() merges them into each family's own.
#
# ndt is in SECONDS, exp(-1.2) = 0.30 s: nothing in these families is unit
# equivariant (see .POUTLIER_SCALE in core_shifted.R). The SD of 0.5 fences
# ndt -> 0 (0.01 s is 6.8 SDs out) without saying where in 0.1-0.8 s the data
# put it; 0.2 was tight enough to cause divergences on data with a 0.6 s
# non-decision time (NEWS, 0.3.2).
#
# poutlier is centred on 0.25%: Miller (2024) puts the outlier rate below
# 0.5%, and hierarchical fits of the untrimmed Illusion Game data sit at -6 to
# -8 on the logit scale (NEWS, 0.3.5). Its `nat` form has its mode at zero -
# leaving the dpar out says the data were trimmed or clean - and the same
# median, exponential(280) against plogis(-6) = 0.0025.
#
# shape has an identity link, so both forms carry one prior, centred on the
# LogNormal and several SDs from sigma * shape = 1, where the density blows up
# at the shift. It is overridden: brms knows the name and supplies
# gamma(0.01, 0.01), positive only, which would truncate the inverse-Weibull
# half of the family at 0.
#' @keywords internal
.SHIFTED_BASE_PRIORS <- list(
  ndt = c(link = "normal(-1.2, 0.5)", nat = "lognormal(-1.2, 0.5)"),
  poutlier = c(link = "normal(-6, 1)", nat = "exponential(280)"),
  shape = c(link = "normal(0, 0.5)", nat = "normal(0, 0.5)")
)


# Priors for the shifted families built on `ndt` + `poutlier`: the two the
# parameterization introduces, plus whatever decision dpars the family named for
# itself in the registry's `prior` slot.
#' @keywords internal
.priors_shifted <- function(formula, data, family, ...) {
  spec <- .mixture_spec(.family_name(family))
  own <- spec$prior
  if (is.null(own)) own <- list()
  # The base entries win, but no family names one of them, so this only decides
  # a collision that cannot currently happen.
  own[names(.SHIFTED_BASE_PRIORS)] <- .SHIFTED_BASE_PRIORS
  .priors_dpars(formula, data, family, own, override = "shape", ...,
                race = isTRUE(spec$race))
}


# Reads get_prior() for the model in hand and fills only the rows brms left
# flat - plus the ones named in `override`, where the default brms supplies is
# worse than nothing - so the result cannot contain a prior matching no
# parameter.
#
# `own` is a named list, one entry per dpar to edit, each a character vector
# with any of `link`, `nat` and `slope`. See .PRIORS_PLAIN above.
#' @keywords internal
.priors_dpars <- function(formula, data, family, own, override = character(0),
                          ..., slope_default = "normal(0, 0.2)", race = FALSE) {
  p <- brms::get_prior(formula, data = data, family = family, ...)
  # Kept whole: once `p` has been filtered down to the rows being set, there is
  # no way left to tell "this blanket row covers a coefficient we left alone"
  # from "this blanket row covers nothing at all".
  all_rows <- p
  dpars <- names(own)
  if (is.null(override)) override <- character(0)

  # Smooths: the SDs follow a rule of their own, on every dpar (.priors_sds()),
  # and the unpenalised columns of a smooth on a dpar the family names get the
  # dpar's slope prior scaled to their basis (.bs_scales()). One
  # make_standata() per distinct term feeds both.
  scales <- .smooth_scales(all_rows, formula, data, ...)
  sds <- .priors_sds(all_rows, scales)
  bs_rms <- .bs_scales(all_rows, scales, dpars)
  bs_row <- p$class == "b" & paste(p$dpar, p$coef, sep = "|") %in% names(bs_rms)

  # A dpar reaches get_prior() in one of two forms. Modelled, even as `~ 1`,
  # it has a linear predictor and rows of class "Intercept" / "b" / "sd" with
  # dpar = "<name>", on the LINK scale. Omitted, brms declares it as a plain
  # auxiliary parameter - the mechanism behind gaussian()'s `sigma` - with
  # class = "<name>" and dpar = "", on the NATURAL scale, so it needs the
  # `nat` prior. Auxiliary rows are replaced whether or not brms filled them,
  # which is what removes the uniform(0, min_Y) on an omitted `ndt` and the
  # defaults brms gives an omitted `shape`, `sigma` or `phi`.
  aux <- p$class %in% dpars & !nzchar(p$dpar)

  link <- p$dpar %in% dpars & !nzchar(p$prior)
  # In a race `mu` is accumulator 0's rate and the dec() coding decides which
  # accumulator that is, so its group-level SD gets what the other
  # accumulator's gets rather than brms's student_t(3, 0, 2.5).
  if (isTRUE(race)) {
    link <- link | (!nzchar(p$dpar) & p$class == "sd" & !nzchar(p$prior))
  }
  # ...unless a proper blanket row already covers them, as brms gives smooths.
  link <- link & !vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  # The modelled counterpart of the auxiliary override: these rows arrive
  # non-empty, so they have to be replaced rather than filled.
  link <- link | (p$dpar %in% override & p$class == "Intercept")

  # `mu` is the response's own linear predictor, reported with an EMPTY dpar
  # (class "Intercept", or "b" with coef "Intercept" under `0 + Intercept`),
  # so neither branch above sees it. A family that names a `mu` prior - the
  # ex-Gaussian and the GEG, whose `mu` is in seconds on an identity link - is
  # asking for its intercept, and only that: the response's slopes are the
  # effects being estimated, and shrinking them is not this function's call.
  resp <- ("mu" %in% dpars) & !nzchar(p$dpar) &
    (p$class == "Intercept" | (p$class == "b" & p$coef == "Intercept"))

  fill <- (aux | link | resp) & p$class != "sds"
  if (!any(fill)) {
    return(if (is.null(sds)) brms::empty_prior() else sds)
  }
  # The dpar to look the prior up under. It is `p$dpar` for everything except
  # the response rows above, whose dpar is empty by construction.
  target <- p$dpar
  target[resp] <- "mu"
  target <- target[fill]
  aux <- aux[fill]
  bs_row <- bs_row[fill]
  p <- p[fill, , drop = FALSE]

  p$prior <- vapply(
    seq_len(nrow(p)),
    function(i) {
      cls <- p$class[i]
      # Natural scale: no link is applied to an auxiliary parameter.
      if (aux[i]) return(.own_prior(own, cls, "nat"))
      intercept <- cls == "Intercept" || (cls == "b" && p$coef[i] == "Intercept")
      if (intercept) {
        .own_prior(own, target[i], "link")
      } else if (cls == "b") {
        # a family may widen the blanket slope prior for a dpar of its own
        slope <- .own_prior(own, target[i], "slope")
        if (!nzchar(slope)) slope <- slope_default
        if (bs_row[i]) {
          # the unpenalised column of a smooth: the same slope prior, stated
          # on the trend in link units rather than on the coefficient
          rms <- bs_rms[[paste(p$dpar[i], p$coef[i], sep = "|")]]
          scaled <- .scale_slope_prior(slope, rms)
          if (nzchar(scaled)) scaled else slope
        } else {
          slope
        }
      } else if (cls == "sd") {
        "exponential(1)"
      } else {
        "" # cor, and anything else brms already gives a proper default
      }
    },
    character(1)
  )
  keep <- nzchar(p$prior)
  bs_row <- bs_row[keep]
  p <- p[keep, , drop = FALSE]

  # get_prior() reports a blanket row (empty coef) alongside one row per
  # coefficient. Setting both makes brms warn that the blanket one is unused, so
  # keep the blanket row and drop what it subsumes. Two exceptions: `b` on a
  # coefficient named "Intercept", which under `0 + Intercept` IS the intercept
  # and carries a different location from the slopes, and the unpenalised
  # column of a smooth, whose slope prior has been rescaled to its basis.
  drop <- vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  drop <- drop & !(p$class == "b" & p$coef == "Intercept") & !bs_row
  kept <- p[!drop, , drop = FALSE]

  # The exception can empty the blanket row out. `ndt ~ 0 + Intercept` has one
  # coefficient, that coefficient is the intercept, and it has just been given a
  # location of its own - leaving the blanket row covering nothing, which brms
  # warns about in turn. So the pair is checked in both directions.
  gone <- vapply(seq_len(nrow(kept)), .blanket_now_unused, logical(1),
                 p = kept, all_rows = all_rows)
  kept <- kept[!gone, , drop = FALSE]
  if (is.null(sds)) kept else rbind(kept, sds)
}


# The wiggliness scale of every smooth (`sds`), on every dpar, as
# exponential(<that term's scale>): exponential(1) on the wiggle in link
# units, whatever basis mgcv built (?cogmod_priors, *Smooth terms*;
# .zs_scales()). The response's own smooths are included: unlike its slopes,
# `sds` is a hyperparameter whose unit the basis sets, and brms's
# student_t(3, 0, 2.5) is as blind to that as a fixed exponential would be.
#
# brms takes one sds prior per term, for all its penalties and `by` levels,
# and warns when its own non-empty blanket row covers nothing. So each dpar's
# blanket row takes the prior of one of its terms, and only the terms whose
# prior differs get a row of their own. Which term is arbitrary, so it is the
# first label in byte order (method = "radix"), not the first row
# get_prior() returns: that order follows the session's collation, which puts
# `s(x)` and `s(x, bs = "cr")` in opposite orders on Windows and on Linux, and
# failed CI once. A term whose scale cannot be computed gets exponential(1).
#' @keywords internal
.priors_sds <- function(p, scales) {
  s <- p[p$class == "sds" & !nzchar(p$resp), , drop = FALSE]
  terms <- unique(s$coef[nzchar(s$coef)])
  if (!length(terms)) return(NULL)
  rate <- vapply(terms, function(tm) {
    v <- scales[[tm]]$sds
    if (is.null(v)) NA_real_ else v
  }, numeric(1))
  prior <- ifelse(is.finite(rate),
                  sprintf("exponential(%s)", .num_str(signif(rate, 3))),
                  "exponential(1)")
  names(prior) <- terms
  rows <- lapply(unique(s$dpar), function(dp) {
    r <- s[s$dpar == dp, , drop = FALSE]
    term <- nzchar(r$coef)
    r$prior[term] <- unname(prior[r$coef[term]])
    blanket <- !term
    if (!any(blanket) || !any(term)) return(r[term, , drop = FALSE])
    first <- r$prior[term][order(r$coef[term], method = "radix")][1]
    r$prior[blanket] <- first
    r[blanket | (term & r$prior != first), , drop = FALSE]
  })
  do.call(rbind, rows)
}


# The scales of every distinct smooth term in a get_prior() table, keyed by
# the label get_prior() gives it, from .smooth_term_scales(). NULL when the
# model has no smooth; a term that cannot be built is NULL in the list. Built
# once per call, since each term costs a make_standata() of its own.
#' @keywords internal
.smooth_scales <- function(p, formula, data, ...) {
  s <- p[p$class == "sds" & !nzchar(p$resp), , drop = FALSE]
  terms <- unique(s$coef[nzchar(s$coef)])
  if (!length(terms)) return(NULL)
  stats::setNames(
    lapply(terms, .smooth_term_scales, formula = formula, data = data, ...),
    terms
  )
}


# The unpenalised columns of the smooths on the dpars a family names, as
# "<dpar>|<coef>" -> root mean square of the column. brms reports them under
# class `b`, so they would take the dpar's slope prior, but mgcv sets their
# scale - an RMS of 0.16 for s(x), 2 for s(x, bs = "cr") - and on an s(x)
# normal(0, 0.2) on the coefficient was normal(0, 0.03) on the trend.
# .priors_dpars() gives each normal(0, s / rms) instead. The response's stay
# with brms, as its slopes do. A column name two smooths share (two smooths of
# one covariate on one dpar) is left to the blanket, since get_prior() reports
# the pair as one row; so is a column whose term could not be built.
#' @keywords internal
.bs_scales <- function(p, scales, dpars) {
  if (is.null(scales)) return(numeric(0))
  s <- p[p$class == "sds" & !nzchar(p$resp) & nzchar(p$coef) &
           p$dpar %in% dpars, , drop = FALSE]
  out <- numeric(0)
  for (dp in unique(s$dpar)) {
    cols <- unlist(lapply(unique(s$coef[s$dpar == dp]), function(tm) {
      scales[[tm]]$bs
    }))
    if (!length(cols)) next
    nm <- names(cols)
    ok <- !is.na(cols) & !nm %in% nm[duplicated(nm)]
    cols <- cols[ok]
    if (!length(cols)) next
    names(cols) <- paste(dp, names(cols), sep = "|")
    out <- c(out, cols)
  }
  out
}


# `normal(0, s)` restated for a column of root mean square `rms`: normal(0,
# s / rms), so that the prior on the coefficient is `normal(0, s)` on the
# coefficient times the column at a typical row. "" if the slope prior is not
# a zero-centred normal, which none of the registries' are, or `rms` is not a
# positive number; the caller then keeps the prior as it stands.
#' @keywords internal
.scale_slope_prior <- function(prior, rms) {
  if (!is.numeric(rms) || length(rms) != 1L || !is.finite(rms) || rms <= 0) {
    return("")
  }
  m <- regmatches(prior, regexec("^normal\\(\\s*0\\s*,\\s*([0-9.eE+-]+)\\s*\\)$",
                                 prior))[[1]]
  if (length(m) != 2L) return("")
  s <- suppressWarnings(as.numeric(m[2]))
  if (!is.finite(s) || s <= 0) return("")
  sprintf("normal(0, %s)", .num_str(signif(s / rms, 3)))
}


# One entry from a family's own `prior` slot in the registry, or "" if the
# family named no prior for that dpar, or named one but not on that scale. The
# slot holds named character vectors, and `[[` on one of those errors rather
# than returning NULL when the name is absent, so both lookups are guarded.
#' @keywords internal
.own_prior <- function(own, dpar, what) {
  if (is.null(own) || !dpar %in% names(own)) return("")
  v <- own[[dpar]]
  if (!what %in% names(v)) return("")
  unname(v[[what]])
}


# Is row i of a get_prior() table subsumed by a blanket row (one with an empty
# `coef`) that already carries a prior? brms warns about the blanket row going
# unused if both are set, so only one of the pair should be filled.
#' @keywords internal
.covered_by_blanket <- function(i, p) {
  if (!nzchar(p$coef[i])) {
    return(FALSE)
  }
  any(
    !nzchar(p$coef) & nzchar(p$prior) & p$class == p$class[i] &
      p$dpar == p$dpar[i] & p$group == p$group[i] & p$resp == p$resp[i]
  )
}


# The mirror image: is blanket row i of the *final* table left with nothing to
# cover, because every coefficient it applies to kept an individual prior?
# `all_rows` is the unfiltered get_prior() table, which is the only place the
# full coefficient list still exists.
#' @keywords internal
.blanket_now_unused <- function(i, p, all_rows) {
  if (nzchar(p$coef[i])) {
    return(FALSE)
  }
  same <- function(d) {
    d$class == p$class[i] & d$dpar == p$dpar[i] & d$group == p$group[i] &
      d$resp == p$resp[i]
  }
  coefs <- all_rows$coef[nzchar(all_rows$coef) & same(all_rows)]
  # No coefficients to cover means this is a blanket over something else
  # entirely (a group-level SD, say), which brms is happy to use as it stands.
  length(coefs) > 0 && all(coefs %in% p$coef[nzchar(p$coef) & same(p)])
}
