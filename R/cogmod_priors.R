#' Priors that make a cogmod posterior proper
#'
#' @description
#' Fills in weakly informative priors for the parameters `brms` would otherwise
#' leave flat, for the model you are actually fitting.
#'
#' **For [cogmod_lognormal()] this is not a convenience.** `brms` assigns a flat,
#' improper prior to the intercept of any custom-family parameter it does not
#' recognise, which there means both `ndt` and `poutlier`. The likelihood has
#' two directions in which it is exactly flat: `poutlier` toward 1, where every
#' response is attributed to the outlier component and `mu`, `sigma` and `ndt`
#' drop out of the density altogether; and `ndt` toward 0, where the model
#' reduces to an unshifted LogNormal and the gradient with respect to
#' `log(ndt)` vanishes. Flat prior plus infinite flat region is an improper
#' posterior. The fit does not fail loudly - it returns intercepts around `1e14`
#' with `Rhat` near 2 and an effective sample size of about 5.
#'
#' The second direction has nothing to do with the mixture; it is inherent to
#' putting a positive shift on a log link, which is why a prior on `poutlier`
#' alone is not enough.
#'
#' @details
#' The function starts from `brms::get_prior()` for the model in hand, edits the
#' rows it knows how to improve, and returns the result of
#' `brms::validate_prior()`. Three things follow.
#'
#' Every row comes from the model `brms` is going to build, so a prior matching
#' no parameter is impossible by construction: `0 + Intercept` formulas,
#' interactions, group-level terms and smooths are all handled because none of
#' them are guessed at.
#'
#' The return value is a `brmsprior` object: the same class returned by
#' [brms::get_prior()] and [brms::prior()]. It prints as a table, including the
#' defaults that `brms` will use, but can be combined with `c()` like any other
#' `brms` prior specification.
#'
#' Passing it through `brms::validate_prior()` means a malformed specification
#' errors here, with the offending row in view, rather than deep inside `brm()`.
#'
#' # Setting your own priors
#'
#' Combine it with [brms::prior()] entries. Set `replace = TRUE` to replace one
#' of these defaults, or omit it to add a prior for a different slot:
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
#' `brms::bf(..., family = cogmod_lognormal())`. Every family built on the direct
#' `ndt` + `poutlier` parameterization is edited: [cogmod_lognormal()],
#' [cogmod_loggamma()], [cogmod_invgaussian()], [cogmod_gamma()], [cogmod_invgamma()],
#' [cogmod_weibull()], [cogmod_invweibull()], [cogmod_logweibull()] and, for the
#' choice-and-RT models, [cogmod_lnr()]. Any other family, or
#' a formula carrying none, is passed through: you get a message and `brms`'s own
#' defaults, unchanged, so the call is always safe to leave in a script.
#'
#' What gets set, on the link scale (`log` for `ndt`, `logit` for `poutlier`,
#' `identity` for `shape`):
#'
#' | class | `ndt` | `poutlier` | `shape` |
#' | --- | --- | --- | --- |
#' | `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-1.2 + log(minrt / 0.3), 0.2)` | `normal(-5, 1)` | `normal(0, 0.5)` |
#' | `b` (slopes) | `normal(0, 0.2)` | `normal(0, 0.2)` | `normal(0, 0.2)` |
#' | `sd`, `sds` | `exponential(1)` | `exponential(1)` | `exponential(1)` |
#'
#' The `ndt` location moves with `minrt`, so the priors stay equivariant to the
#' unit of measurement in the same way the likelihood does; at the default it is
#' `normal(-1.2, 0.2)`, roughly 170 to 300 ms. `poutlier` is a proportion and
#' does not move: `normal(-5, 1)` is centred at about 0.7% and puts roughly 95%
#' of its mass between 0.1% and 5%.
#'
#' `shape` exists only for [cogmod_loggamma()]. `normal(0, 0.5)` is centred on the
#' LogNormal shape and keeps the sampler clear of `sigma * shape >= 1`, where the
#' decision density becomes unbounded at the shift and the likelihood with it.
#'
#' A family may add rows of its own where its likelihood has a flat direction
#' that `brms` would leave improper. Two do:
#'
#' - [cogmod_lba1()]: `sigmabias` and `boundary`. As the start-point range
#'   approaches zero the LBA converges to the recinormal and the likelihood
#'   stops depending on it, and a `softplus` link reaches zero only at minus
#'   infinity. Without those rows `sigmabias` runs off - `softplus(-10.4)`,
#'   `Rhat` 1.69.
#' - [cogmod_lnr()]: `nuone`, `sigmazero` and `sigmaone`. Push an accumulator's
#'   rate far enough down and it never finishes first, so the density depends on
#'   it only through a survival term that has already saturated at 1; past about
#'   `nuone = -6` the log-likelihood is exactly constant, and that accumulator's
#'   `sigma` is unidentified along with it. `mu` - which is `nuzero` - has the
#'   mirror-image plateau, but it is the response's own intercept and `brms`
#'   already gives it a proper `student_t` default, so it is left alone. Model a
#'   rarely-chosen option and it is worth mirroring the `nuone` prior onto it.
#'
#' Both are the same failure as `ndt` and `poutlier`: an infinite flat region
#' under a flat prior. See [cogmod_lba1()] and [cogmod_lnr()].
#'
#' Note that no prior is set on the shape of [cogmod_weibull()] or
#' [cogmod_gamma()], although a shape below 2 makes their `ndt` gradient
#' unbounded. That region is reached because the *likelihood* prefers it, by
#' around 100 log units on the data in `vignette("rt_models")`, so a prior weak
#' enough to be a sensible default cannot move the posterior out of it - only
#' bias it. `?rcogmod_weibull` sets out what to do instead.
#'
#' # Parameters left out of the formula
#'
#' Writing `ndt ~ 1` and omitting `ndt` entirely are not the same thing to
#' `brms`, and the difference matters here. A dpar that appears in `bf()` - even
#' as `~ 1` - gets a linear predictor, so it is estimated on the **link** scale
#' and reported under *Regression Coefficients* as `ndt_Intercept`. A dpar left
#' out is declared as a plain auxiliary parameter instead, the same mechanism as
#' `sigma` for `gaussian()`, estimated on the **natural** scale with no link and
#' reported under *Further Distributional Parameters*.
#'
#' Both forms are filled, with priors on whichever scale the parameter actually
#' lives on:
#'
#' | dpar | in `bf()` (link scale) | omitted (natural scale) |
#' | --- | --- | --- |
#' | `ndt` | `normal(-1.2 + log(minrt / 0.3), 0.2)` | `lognormal(-1.2 + log(minrt / 0.3), 0.2)` |
#' | `poutlier` | `normal(-5, 1)` | `exponential(100)` |
#' | `shape` | `normal(0, 0.5)` | `normal(0, 0.5)` |
#' | `sigmabias`, `boundary` ([cogmod_lba1()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#' | `nuone` ([cogmod_lnr()]) | `normal(0.7, 1.5)` | `normal(0.7, 1.5)` |
#' | `sigmazero`, `sigmaone` ([cogmod_lnr()]) | `normal(0, 1)` | `lognormal(-0.7, 0.75)` |
#'
#' The `ndt` pair describes the same belief twice: `lognormal` is just `normal`
#' on the log scale, written for the untransformed parameter.
#'
#' If the data were trimmed before fitting, tighten this rather than removing
#' it: `normal(-7, 0.5)` asserts essentially no contamination while keeping the
#' density positive below `ndt`. Fixing `poutlier = 0` outright reinstates the
#' hard min-RT boundary that the component exists to remove. See the *Trimmed
#' data* section of [cogmod_lognormal()].
#'
#' `poutlier` is deliberately *not* the same belief twice. Leaving it out of the
#' formula is itself information - you either trimmed the data already or do not
#' expect outliers - so the omitted form puts its **mode at zero**, which a
#' logit-scale prior cannot do at any location. The centre is unchanged:
#' `exponential(100)` has median 0.0069 against `plogis(-5) = 0.0067`. It is
#' still a prior rather than a constraint, so a genuine spike of fast responses
#' will still pull the rate up; to switch the parameter off entirely, trim the
#' data and it will simply sit near zero.
#'
#' The omitted-`ndt` row is the one case where a **non-empty** `brms` default is
#' overridden rather than filled. `brms` recognises the name `ndt` from its own
#' shifted families and supplies `uniform(0, min_Y)` - precisely the min-RT
#' bound that [cogmod_lognormal()]'s parameterization exists to remove, reimposed
#' silently and with a warning about an upper bound on an unbounded parameter.
#' An omitted `poutlier` is left flat over `[0, 1]` by `brms`, which is proper
#' but puts half its mass above 0.5, and an omitted `shape` is flat over the whole
#' real line, which is not proper at all.
#'
#' Slope and group-level priors are deliberately narrow. On a log or a logit
#' link a flat slope prior is not as harmless as it looks, and a group-level SD
#' with no prior can wander far enough for individual groups to reach the flat
#' regions above even when the population intercept is well behaved.
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = cogmod_lognormal())`.
#' @param data The data, as passed to `brms::brm()`.
#' @param ... Passed to `brms::get_prior()` and `brms::validate_prior()`, for
#'   arguments such as `data2` or `knots`.
#'
#' @return A `brmsprior` object, to pass to `brms::brm(prior = )`.
#'
#' @seealso [cogmod_inits()], [cogmod_stanvars()]
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
cogmod_priors <- function(formula, data, ...) {
  family <- .cogmod_family(formula)
  fam <- .family_name(family)

  if (!isTRUE(fam %in% .OUTLIER_FAMILIES)) {
    message(
      "cogmod_priors() has nothing to add for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'; returning the brms defaults unchanged. ",
      "Currently supported: ", paste(.OUTLIER_FAMILIES, collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_lognormal())."
    )
    args <- list(brms::empty_prior(), formula, data, ...)
    if (!is.null(family)) args$family <- family
    return(do.call(brms::validate_prior, args))
  }

  out <- .priors_shifted(formula, data, family, ...)
  brms::validate_prior(out, formula, data, family = family, ...)
}


# Priors for the shifted families built on `ndt` + `poutlier`. Reads get_prior()
# and fills only the rows brms left flat, so the result cannot contain a prior
# matching no parameter.
#' @keywords internal
.priors_shifted <- function(formula, data, family, ...) {
  p <- brms::get_prior(formula, data = data, family = family, ...)
  # Kept whole: once `p` has been filtered down to the rows being set, there is
  # no way left to tell "this blanket row covers a coefficient we left alone"
  # from "this blanket row covers nothing at all".
  all_rows <- p

  # `shape` only exists for cogmod_loggamma(); brms leaves it flat like the others, and
  # a flat prior there lets the sampler wander into sigma * shape >= 1, where the
  # decision density is unbounded at the shift.
  #
  # A family may name further decision dpars of its own in the registry's
  # `prior` slot, for the same reason: a direction the likelihood is flat or
  # near-flat in, which brms would otherwise leave improper. cogmod_lba1() is
  # the case in hand - see the note on its registry entry.
  own <- .mixture_spec(.family_name(family))$prior
  dpars <- unique(c("ndt", "poutlier", "shape", names(own)))

  # A dpar reaches get_prior() in one of two forms, depending on whether it
  # appears in the formula at all.
  #
  #  - *Modelled*, even as `~ 1`: brms builds a linear predictor for it and
  #    reports class "Intercept" / "b" / "sd" with dpar = "<name>". The
  #    parameter lives on the LINK scale (log for ndt, logit for poutlier).
  #
  #  - *Omitted*: brms declares it as a plain auxiliary parameter - the same
  #    mechanism as `sigma` for gaussian() when you do not write `sigma ~ .` -
  #    and reports class = "<name>" with dpar = "". It is declared with the
  #    family's own bounds and NO link applied, so it lives on the NATURAL
  #    scale and needs a different prior, not the same one.
  #
  # Matching on `dpar` alone caught only the first kind, which left an omitted
  # `poutlier` on brms's flat prior over [0, 1] and an omitted `shape` genuinely
  # improper.
  aux <- p$class %in% dpars & !nzchar(p$dpar)

  link <- p$dpar %in% dpars & !nzchar(p$prior)
  # ...unless a proper blanket row already covers them, as brms gives smooths.
  link <- link & !vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)

  # brms recognises the name `shape` from its own gamma/weibull families and
  # supplies gamma(0.01, 0.01) - a prior with support only on the positives, for
  # a parameter that has to be free to go negative (that is the whole inverse
  # Weibull half of the family). Stan would reject every negative proposal,
  # silently truncating `shape` at 0. That row arrives non-empty, so it has to be
  # overridden rather than filled. Same story as the `uniform(0, min_Y)` that
  # brms puts on an omitted `ndt`, handled through `aux` above.
  link <- link | (p$dpar == "shape" & p$class == "Intercept")

  fill <- aux | link
  if (!any(fill)) {
    return(brms::empty_prior())
  }
  aux <- aux[fill]
  p <- p[fill, , drop = FALSE]

  # ndt is a location in time and moves with the timescale; poutlier is a
  # proportion and does not.
  shift <- log(.validate_minrt(.as_minrt(family)) / eval(formals(.validate_minrt)$minrt))
  loc_ndt <- formatC(-1.2 + shift, format = "g", digits = 4, width = 1)

  p$prior <- vapply(
    seq_len(nrow(p)),
    function(i) {
      cls <- p$class[i]
      if (aux[i]) {
        # Natural scale: no link is applied to an auxiliary parameter.
        return(switch(
          cls,
          # The same distribution as normal(loc, 0.2) on log(ndt), written for
          # the untransformed parameter. This deliberately *overrides* a
          # non-empty brms default: brms recognises the name `ndt` from its own
          # shifted families and supplies uniform(0, min_Y), which is exactly
          # the min-RT bound this parameterization exists to remove.
          ndt = sprintf("lognormal(%s, 0.2)", loc_ndt),
          # Leaving `poutlier` out of the formula is itself information: the
          # user either trimmed the data already or does not expect outliers.
          # So the mode sits at zero rather than at a positive rate, unlike the
          # link-scale prior, which cannot put mass at 0 at all. The median
          # matches it regardless - exponential(100) has median 0.0069 against
          # plogis(-5) = 0.0067 - so the centre is unchanged, only the shape.
          # Truncation at 1 is negligible: exp(-100) of the mass is above it.
          # brms's own default here is flat over [0, 1] - proper, but it puts
          # half its mass above 0.5.
          poutlier = "exponential(100)",
          # identity link, so the link-scale prior carries over unchanged
          shape = "normal(0, 0.5)",
          # anything the family named for itself
          .own_prior(own, cls, "nat")
        ))
      }
      intercept <- cls == "Intercept" || (cls == "b" && p$coef[i] == "Intercept")
      if (intercept) {
        switch(
          p$dpar[i],
          ndt = sprintf("normal(%s, 0.2)", loc_ndt),
          # Centred on the LogNormal, and tight enough that sigma * shape = 1 -
          # where the density blows up at the shift - is several SDs away for
          # any realistic sigma.
          shape = "normal(0, 0.5)",
          poutlier = "normal(-5, 1)",
          # anything the family named for itself
          .own_prior(own, p$dpar[i], "link")
        )
      } else if (cls == "b") {
        # a family may widen the blanket slope prior for a dpar of its own
        slope <- .own_prior(own, p$dpar[i], "slope")
        if (nzchar(slope)) slope else "normal(0, 0.2)"
      } else if (cls %in% c("sd", "sds")) {
        "exponential(1)"
      } else {
        "" # cor, and anything else brms already gives a proper default
      }
    },
    character(1)
  )
  p <- p[nzchar(p$prior), , drop = FALSE]

  # get_prior() reports a blanket row (empty coef) alongside one row per
  # coefficient. Setting both makes brms warn that the blanket one is unused, so
  # keep the blanket row and drop what it subsumes. `b` on a coefficient named
  # "Intercept" is the exception: under `0 + Intercept` that IS the intercept and
  # carries a different location from the slopes.
  drop <- vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  drop <- drop & !(p$class == "b" & p$coef == "Intercept")
  kept <- p[!drop, , drop = FALSE]

  # The exception can empty the blanket row out. `ndt ~ 0 + Intercept` has one
  # coefficient, that coefficient is the intercept, and it has just been given a
  # location of its own - leaving the blanket row covering nothing, which brms
  # warns about in turn. So the pair is checked in both directions.
  gone <- vapply(seq_len(nrow(kept)), .blanket_now_unused, logical(1),
                 p = kept, all_rows = all_rows)
  kept[!gone, , drop = FALSE]
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
