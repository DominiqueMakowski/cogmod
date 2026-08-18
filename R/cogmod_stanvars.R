#' The Stan code a cogmod family needs, read off the model
#'
#' @description
#' Returns the `stanvars` argument for `brms::brm()`, for whichever `cogmod`
#' family the model uses. It is a front end to the per-family
#' `<family>_stanvars()` functions - `cogmod_lognormal_stanvars()`,
#' `cogmod_choco_stanvars()`, `cogmod_ddm_stanvars()` and the rest - which remain available
#' and unchanged.
#'
#' @details
#' `brms` needs the custom likelihood injected into the generated Stan program,
#' and every `cogmod` family ships one. Calling this instead of the family's own
#' function means the family is named once, in `bf()`, rather than twice:
#'
#' ```r
#' f <- brms::bf(RT ~ Condition, ndt ~ Condition,
#'               family = cogmod_lognormal(minrt = 0.25))
#'
#' brms::brm(f, data = df,
#'           prior    = cogmod_priors(f, df),
#'           init     = cogmod_inits(f, df),
#'           stanvars = cogmod_stanvars(f))
#' ```
#'
#' # Why this is not only tidier
#'
#' For the shifted families, `minrt` is *baked into the generated Stan code as a
#' literal*, because a Stan function cannot see the data block. It is therefore
#' possible to fit `cogmod_lognormal(minrt = 0.25)` against Stan code compiled with
#' the default `0.3` - the model runs, and quietly uses a different outlier
#' component from the one the family describes, which nothing downstream will
#' flag.
#'
#' `cogmod_stanvars()` takes `minrt` off the family, so the Stan code and the
#' family cannot disagree. Passing `minrt` explicitly still overrides it, for
#' the rare case where that is wanted.
#'
#' # What it accepts
#'
#' A `brms::bf()` formula carrying the family, the family object itself, or a
#' fitted `brmsfit` (useful for recompiling or for `update()`).
#'
#' @param formula A `brms::bf()` formula carrying the family, a `cogmod` family
#'   object, or a fitted `brmsfit`.
#' @param ... Passed to the family's own `<family>_stanvars()` function, for
#'   arguments such as `minrt`.
#'
#' @return A `stanvars` object, to pass to `brms::brm(stanvars = )`.
#'
#' @seealso [cogmod_priors()], [cogmod_inits()]
#'
#' @examples
#' f <- brms::bf(RT ~ 1, ndt ~ 1, family = cogmod_lognormal(minrt = 0.25))
#' cogmod_stanvars(f)
#'
#' # Equivalent to naming the family a second time:
#' cogmod_lognormal_stanvars(minrt = 0.25)
#'
#' @export
cogmod_stanvars <- function(formula, ...) {
  family <- .cogmod_family(formula)
  fam <- .family_name(family)
  if (!isTRUE(fam %in% .cogmod_families())) {
    stop(
      "cogmod_stanvars() has no Stan code for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'. Supported: ", paste(.cogmod_families(), collapse = ", "), ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_lognormal()).",
      call. = FALSE
    )
  }

  fn <- get(paste0(fam, "_stanvars"), envir = asNamespace("cogmod"),
            mode = "function")
  args <- list(...)
  # The shifted families bake `minrt` into the Stan code as a literal, so it
  # has to come from the family the model is actually using rather than from
  # the default. An explicit argument still wins.
  if ("minrt" %in% names(formals(fn)) && is.null(args$minrt) &&
      !is.null(family$minrt)) {
    args$minrt <- family$minrt
  }
  do.call(fn, args)
}


# The families that ship Stan code. Built from the shifted-family registry so
# the two cannot drift apart.
#' @keywords internal
.cogmod_families <- function() {
  # unique() because the choice families that have been migrated to the outlier
  # mixture (cogmod_lnr, so far) are in .OUTLIER_FAMILIES as well.
  unique(c(
    "cogmod_betadiscrete", "cogmod_betagate", "cogmod_choco", "cogmod_ddm", "cogmod_lba2", "cogmod_lnr", "cogmod_rdm",
    "cogmod_exgaussian", .OUTLIER_FAMILIES
  ))
}
