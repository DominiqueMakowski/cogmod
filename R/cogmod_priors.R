#' Priors that make a cogmod posterior proper
#'
#' @description
#' Fills in weakly informative priors for the parameters `brms` would otherwise
#' leave flat, for the model you are actually fitting.
#'
#' **For [rt_lognormal()] this is not a convenience.** `brms` assigns a flat,
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
#' `brms::bf(..., family = rt_lognormal())`. Only [rt_lognormal()] is edited.
#' Any other family, or a formula carrying none, is passed through: you get a
#' message and `brms`'s own defaults, unchanged, so the call is always safe to
#' leave in a script.
#'
#' What gets set for `rt_lognormal()`, on the link scale (`log` for `ndt`,
#' `logit` for `poutlier`):
#'
#' | class | `ndt` | `poutlier` |
#' | --- | --- | --- |
#' | `Intercept`, or `b` on a coefficient named `Intercept` | `normal(-1.5 + log(minrt / 0.3), 0.15)` | `normal(-5, 1)` |
#' | `b` (slopes) | `normal(0, 1)` | `normal(0, 1)` |
#' | `sd`, `sds` | `exponential(1)` | `exponential(1)` |
#'
#' The `ndt` location moves with `minrt`, so the priors stay equivariant to the
#' unit of measurement in the same way the likelihood does; at the default it is
#' `normal(-1.5, 0.15)`, roughly 170 to 300 ms. `poutlier` is a proportion and
#' does not move: `normal(-5, 1)` is centred at about 0.7% and puts roughly 95%
#' of its mass between 0.1% and 5%.
#'
#' Slope and group-level priors are deliberately narrow. On a log or a logit
#' link a flat slope prior is not as harmless as it looks, and a group-level SD
#' with no prior can wander far enough for individual groups to reach the flat
#' regions above even when the population intercept is well behaved.
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = rt_lognormal())`.
#' @param data The data, as passed to `brms::brm()`.
#' @param ... Passed to `brms::get_prior()` and `brms::validate_prior()`, for
#'   arguments such as `data2` or `knots`.
#'
#' @return A `brmsprior` object, to pass to `brms::brm(prior = )`.
#'
#' @examples
#' d <- data.frame(RT = rrt_lognormal(50, ndt = 0.3, poutlier = 0.02))
#' f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1,
#'   family = rt_lognormal()
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
  family <- formula[["family"]]
  # brms families carry $name; stats families only $family.
  fam <- if (is.character(family)) family else family[["name"]]
  if (is.null(fam) && !is.null(family)) fam <- family[["family"]]

  if (!identical(fam, "rt_lognormal")) {
    message(
      "cogmod_priors() has nothing to add for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'; returning the brms defaults unchanged. ",
      "Currently supported: rt_lognormal. ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = rt_lognormal())."
    )
    args <- list(brms::empty_prior(), formula, data, ...)
    if (!is.null(family)) args$family <- family
    return(do.call(brms::validate_prior, args))
  }

  out <- .priors_rt_lognormal(formula, data, family, ...)
  brms::validate_prior(out, formula, data, family = family, ...)
}


# Priors for the shifted LogNormal. Reads get_prior() and fills only the rows
# brms left flat, so the result cannot contain a prior matching no parameter.
#' @keywords internal
.priors_rt_lognormal <- function(formula, data, family, ...) {
  p <- brms::get_prior(formula, data = data, family = family, ...)

  fill <- p$dpar %in% c("ndt", "poutlier") & !nzchar(p$prior)
  # ...unless a proper blanket row already covers them, as brms gives smooths.
  fill <- fill & !vapply(seq_len(nrow(p)), .covered_by_blanket, logical(1), p = p)
  if (!any(fill)) {
    return(brms::empty_prior())
  }
  p <- p[fill, , drop = FALSE]

  # ndt is a location in time and moves with the timescale; poutlier is a
  # proportion and does not.
  shift <- log(.validate_minrt(.as_minrt(family)) / eval(formals(.validate_minrt)$minrt))
  loc_ndt <- formatC(-1.5 + shift, format = "g", digits = 4, width = 1)

  p$prior <- vapply(
    seq_len(nrow(p)),
    function(i) {
      cls <- p$class[i]
      intercept <- cls == "Intercept" || (cls == "b" && p$coef[i] == "Intercept")
      if (intercept) {
        if (p$dpar[i] == "ndt") {
          sprintf("normal(%s, 0.15)", loc_ndt)
        } else {
          "normal(-5, 1)"
        }
      } else if (cls == "b") {
        "normal(0, 1)"
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
  p[!drop, , drop = FALSE]
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
