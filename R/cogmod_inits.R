#' Starting values that keep the sampler out of the flat regions
#'
#' @description
#' Builds an `init` argument for `brms::brm()`, for the custom families whose
#' default starting point is a bad one.
#'
#' **This is not a tuning knob to reach for when sampling looks bad.** For the
#' families with a non-decision time the default start is actively harmful,
#' and the failure is silent: the chain does not error, it simply never moves.
#'
#' @details
#' Two regions of the parameter space have no gradient to escape on, and both
#' are easy to start inside.
#'
#' - **`ndt` too large.** `brms` initialises on the *unconstrained* scale, so
#'   `init = 0` puts `ndt` at `exp(0) = 1` second, above nearly every
#'   sub-second response. Every response is then attributed to the outlier
#'   component, the decision parameters drop out of the density, and their
#'   gradient is exactly zero.
#' - **Shape below 1.** For [cogmod_gamma()] and [cogmod_weibull()] the density
#'   is unbounded at `ndt` whenever the shape is below 1, and `init = 0` starts
#'   a `softplus` shape at `log(2) = 0.69`, inside that region.
#'
#' No single scalar `init` avoids both - `ndt = exp(c)` wants `c` around
#' `-1.6`, `shape = softplus(c)` wants it above `1.9` - so the parameters have
#' to be set separately, which is what this function does.
#'
#' # How it works
#'
#' The parameter declarations are read off [brms::make_stancode()] for the
#' model you are fitting, and their dimensions off [brms::make_standata()], so
#' `0 + Intercept`, interactions, group-level terms and smooths are matched
#' however `brms` names and sizes them. **Every** declared parameter gets a
#' value, which also spares you CmdStan's message about a partial list:
#'
#' - **The family's values** go to the intercept of each distributional
#'   parameter, on its link scale - through the links on the family, so
#'   `cogmod_gamma(link_mu = "log")` is honoured - and, on the natural scale,
#'   to any parameter left out of the formula, which `brms` declares as a plain
#'   auxiliary parameter. Under `0 + Intercept` the value goes to the
#'   coefficient named `Intercept`.
#' - **`ndt`** starts at half the first percentile of the observed response
#'   times: below essentially every response, but on the data's own scale. The
#'   two errors are not symmetric: a start too small only means the shift has
#'   to grow, which the gradient does, while one too large removes the
#'   gradient.
#' - **The error accumulator** of [cogmod_rdm()], [cogmod_lba2()] and
#'   [cogmod_lnr()] starts slower than the correct one (`driftone = 1` against
#'   `mu = 3`, `nuone = 0` against `mu = 0.7`): few data sets are a 50/50
#'   race, and for the Wald and the LBA a start too fast costs far more than
#'   one too slow. Each family's page has the details.
#' - **Everything else** gets a generic value at least as good as Stan's own
#'   `U(-2, 2)`: zero for slopes and standardized effects, just above the bound
#'   for a positive parameter, the midpoint for a doubly bounded one, the
#'   identity for a Cholesky factor. A smooth's SD (`sds`) starts where its
#'   smooth is nearly flat, at 0.1 link units of wiggle, a different value on
#'   every basis; see *Smooth terms* in [cogmod_priors()].
#'
#' Each chain gets a different draw. The noise is added on the unconstrained
#' scale, so a value can never leave its bounds. For a regression coefficient
#' it is divided by the largest row norm of the coefficient block's design
#' columns, wherever that exceeds 1, so that at no row does the SD of what
#' reaches the linear predictor exceed `jitter`: a slope on age in years would
#' otherwise start with an SD of six link units at the oldest participant.
#'
#' # Supported families
#'
#' The family is read off `formula`, so build it with
#' `brms::bf(..., family = cogmod_gamma())`. The families are those of
#' [cogmod_priors()], listed in its *Supported families* section: every family
#' with an `ndt` + `poutlier` component; [cogmod_exgaussian()] and
#' [cogmod_geg()], whose parameters are lengths of time that the default start
#' puts at 0 s or 0.69 s; and the rating-scale families [cogmod_choco()],
#' [cogmod_betagate()] and [cogmod_betadiscrete()], where it describes a scale
#' nobody uses - half the responses on the midpoint or at an endpoint, or a
#' U-shaped Beta. Any other family is an error.
#'
#' @param formula The model formula, as passed to `brms::brm()`. Must carry the
#'   family, i.e. be built with `brms::bf(..., family = cogmod_gamma())`. May
#'   be left `NULL` only when `warmstart` is a `brmsfit`, whose formula is then
#'   used.
#' @param data The data, as passed to `brms::brm()`. May be left `NULL` only
#'   when `warmstart` is a `brmsfit`, whose data are then used.
#' @param jitter SD of the noise added on the unconstrained scale, so that
#'   chains start at different points. One number sets the population-level
#'   blocks (intercepts and slopes; for a slope, divided by how far its design
#'   columns reach) and a fifth of it the group-level and smooth blocks (`z_*`,
#'   `zs_*`, `sd_*`, `sds_*`), where a unit of noise is multiplied through a
#'   scale and a design column, once per group or basis function. Two numbers
#'   set the two tiers, population first. `0` gives identical starts. `NULL`
#'   (the default) means 0.25, or 0.05 with a `warmstart`, whose values come
#'   from a converged posterior.
#' @param warmstart A previous fit to start from instead of the family's
#'   generic values: a `brmsfit`, a [cogmod_warmstart()] object, the data frame
#'   `as.data.frame()` makes of one, or the path to a CSV file of it. The
#'   starting values are its posterior means, mapped onto this model by
#'   parameter name (a pilot on fewer participants included; see
#'   [cogmod_warmstart()]); whatever it has no value for keeps the value this
#'   function would give it anyway.
#' @param ... Passed to `brms::make_stancode()` and `brms::make_standata()`, for
#'   arguments such as `data2`.
#'
#' @return A function of one argument, suitable for `brms::brm(init = )`. Each
#'   call returns a named list of starting values, one for every parameter the
#'   generated Stan program declares.
#'
#' @seealso [cogmod_priors()], [cogmod_stanvars()], [cogmod_warmstart()]
#'
#' @examples
#' d <- data.frame(RT = rcogmod_gamma(50, ndt = 0.3, poutlier = 0.02))
#' f <- brms::bf(RT ~ 1, sigma ~ 1, ndt ~ 1, poutlier ~ 1, family = cogmod_gamma())
#' inits <- cogmod_inits(f, d)
#' inits(1)
#'
#' # The rating-scale families are covered too. `pmid` starts at 0.05 rather
#' # than the logit origin's 0.5, which would put half of every response
#' # exactly on the midpoint of the scale.
#' r <- data.frame(y = rcogmod_choco(50, pmid = 0.05))
#' g <- brms::bf(y ~ 1, pmid ~ 1, family = cogmod_choco())
#' cogmod_inits(g, r, jitter = 0)(1)
#'
#' \donttest{
#' # Fitting needs cmdstanr, which lives outside CRAN - see the package website.
#' if (requireNamespace("cmdstanr", quietly = TRUE) &&
#'     !is.null(cmdstanr::cmdstan_version(error_on_NA = FALSE))) {
#'   m <- brms::brm(f,
#'     data = d, prior = cogmod_priors(f, d),
#'     stanvars = cogmod_stanvars(f), init = cogmod_inits(f, d),
#'     backend = "cmdstanr", chains = 1, iter = 500, refresh = 0
#'   )
#' }
#' }
#'
#' @export
cogmod_inits <- function(
  formula = NULL,
  data = NULL,
  jitter = NULL,
  warmstart = NULL,
  ...
) {
  # A warm start is cogmod_warmstart()'s job: it builds this function's plan
  # for the model, then writes the previous fit's posterior means over it.
  # Only there can `formula` or `data` be left out, taken from the fit.
  if (!is.null(warmstart)) {
    ws <- cogmod_warmstart(
      warmstart,
      formula = formula,
      data = data,
      jitter = if (is.null(jitter)) 0.05 else jitter,
      ...
    )
    return(ws$init)
  }
  if (is.null(formula) || is.null(data)) {
    stop(
      "`formula` and `data` are required, unless `warmstart` is a brmsfit ",
      "to take them from.",
      call. = FALSE
    )
  }
  jitter <- .init_jitter(jitter)

  family <- .cogmod_family(formula)
  fam <- .family_name(family)
  targets <- if (is.null(fam)) NULL else .init_targets(family)
  if (is.null(targets)) {
    stop(
      "cogmod_inits() has nothing to offer for family '",
      if (is.null(fam)) "<none found on the formula>" else fam,
      "'. Supported: ",
      paste(.init_families(), collapse = ", "),
      ". ",
      "The family is read off the formula, so build it with ",
      "bf(..., family = cogmod_gamma()).",
      call. = FALSE
    )
  }
  links <- .family_links(family)

  # Whatever brms decided to call the parameters, and however large it decided
  # to make them, that is what we work from. Warnings are suppressed because
  # only the structure is wanted: the code is generated with brms' own default
  # priors, which is exactly the situation cogmod_priors() exists to fix, and
  # complaining about it here would be noise.
  code <- suppressWarnings(
    brms::make_stancode(formula, data = data, family = family, ...)
  )
  sdata <- suppressWarnings(
    brms::make_standata(formula, data = data, family = family, ...)
  )
  if (!is.null(targets$ndt)) {
    targets$ndt <- .ndt_start(sdata$Y, targets$ndt)
  }
  plan <- .init_plan(.stan_param_decls(code), as.list(sdata), targets, links)
  .init_fun(plan, jitter)
}


# Where `ndt` starts: half the first percentile of the observed responses, in
# seconds. Below the data, so that essentially every response keeps its
# gradient, but on the data's own scale: a fixed start half a second below a
# slow data set's non-decision time costs thousands of log-density units,
# which a cold chain's first trajectory turns into momentum (NEWS, 0.3.2). The
# jitter is on the log scale, so a jittered start stays below the first
# percentile. A per-group refinement was weighed, and data-aware starts for
# every other parameter measured, and both turned down; see "Tried and decided
# against" in AGENTS.md.
#' @keywords internal
.ndt_start <- function(y, fallback = 0.1) {
  y <- y[is.finite(y) & y > 0]
  if (!length(y)) {
    return(fallback)
  }
  0.5 * stats::quantile(y, 0.01, names = FALSE)
}


# The `init` function brm() gets: one draw of the plan per call, each entry
# jittered on the scale Stan samples on and shaped as declared. Shared with
# cogmod_warmstart(), whose plan carries a previous fit's means.
#' @keywords internal
.init_fun <- function(plan, jitter) {
  jitter <- .init_jitter(jitter)
  function(chain_id = 1) {
    out <- lapply(plan, function(e) {
      v <- e$value
      j <- jitter[[e$tier]]
      # per element, for a coefficient block whose design reaches past 1
      if (!is.null(e$scale)) {
        j <- j * e$scale
      }
      if (any(j > 0)) {
        v <- switch(
          e$kind,
          bounds = .jitter_bounded(v, e$lower, e$upper, j),
          sorted = sort(v + stats::rnorm(length(v), 0, j)),
          v # "fixed": a structured value that jitter would invalidate
        )
      }
      if (length(e$dim) > 1) {
        dim(v) <- e$dim
      }
      v
    })
    stats::setNames(out, vapply(plan, `[[`, character(1), "name"))
  }
}


# Natural-scale starting values ------------------------------------------

# Which families cogmod_inits() has targets for.
#' @keywords internal
.init_families <- function() c(.OUTLIER_FAMILIES, names(.INIT_PLAIN))


# Families with no ndt + poutlier component whose dpars the unconstrained
# origin puts somewhere no data set lives. One natural-scale value per dpar,
# as in the registries' `init` slot: .init_value() sends it through the
# family's link for a modelled intercept, and uses it as it stands for a dpar
# left out of bf().
#' @keywords internal
.INIT_PLAIN <- list(
  # In seconds: the Gaussian centre, its SD and the exponential mean. The
  # origin puts `sigma` and `tau` at softplus(0) = 0.69 s, `sigma` alone wider
  # than most whole RT distributions, and `mu` (identity) at 0 s.
  cogmod_exgaussian = list(mu = 0.4, sigma = 0.1, tau = 0.2),

  # The ex-Gaussian it nests: shape = 1 is that model exactly, and the one
  # point on the shape axis known not to be up the shape/mu ridge.
  cogmod_geg = list(mu = 0.4, sigma = 0.1, tau = 0.2, shape = 1),

  # The rating-scale families are a weaker case than the RT ones: no start
  # here is on a flat region, and a chain begun at the origin does move. It
  # first has to walk back a start that describes a rating scale nobody uses,
  # which costs warmup rather than the fit. Each value was chosen on the
  # log-likelihood shortfall per observation of a start against simulated data
  # sets' own parameters, over several response styles.
  #
  # cogmod_choco(). The precisions (softplus) start at 2: the Beta shapes are
  # conf * prec * 2, so the origin's 0.69 is a U, unbounded at both ends of
  # each half - the one shape rating data almost never has - and 2 a gentle
  # hump. The case is the U rather than the average: over five 4000-trial
  # styles the worst shortfall is flat from prec 1 to 2 (1.00, 1.12) and
  # climbs after (2.25 at 6, hence not ?rcogmod_choco's 4). `pmid` starts at
  # 0.05, where sliders with a visible midpoint sit: the origin's 0.5, half of
  # every response on one value, is the costliest part of the default start
  # (983 of its 1871 units on 2000 trials with 5% midpoints). `pex` matters
  # least, since the observed share of 0s and 1s pins it almost at once, and
  # starts at ?rcogmod_choco's 0.1. The rest start at 0.5, where the origin
  # already puts them, and are named so that a non-default link is followed.
  cogmod_choco = list(
    mu = 0.5,
    confright = 0.5,
    precright = 2,
    confleft = 0.5,
    precleft = 2,
    pex = 0.1,
    bex = 0.5,
    pmid = 0.05
  ),

  # cogmod_betagate(). The origin puts `phi` (softplus) at 0.69, a U-shaped
  # Beta (shapes mu * phi * 2), and `pex` at half of all responses on an
  # endpoint. `phi` starts at 1, the Uniform at mu = 0.5, as in
  # cogmod_betadiscrete(). `pex` starts at 0.2, not choco's 0.1: over eight
  # styles the worst shortfall is 1.23 here against 1.65 at 0.1, where two
  # thirds 1s sets it, and 1.36 at the origin. Fitted cold the start made no
  # difference to sampling; it is here so that every family cogmod_priors()
  # covers is covered (benchmarks/betagate_start/).
  cogmod_betagate = list(mu = 0.5, phi = 1, pex = 0.2, bex = 0.5),

  # cogmod_betadiscrete() bins a Beta into k categories, so no start is
  # unbounded. `phi` (log) starts at the origin's 1, the discrete Uniform with
  # mu = 0.5, and is named because brms's real default, init = "random", draws
  # it anywhere in 0.14-7.4, and a concentrated start starves the categories
  # the data occupy: over six 4000-trial styles the worst shortfall is 1.08 at
  # phi = 1, 1.68 at 2 and 2.87 at 4. `pzero`, the hurdle on the extra "0"
  # category, starts at 0.05 rather than the origin's 0.5, half of every
  # response off the scale (worst case 0.38 against 1.01 at 0.5 and 0.53 at
  # 0.02). A scale with no zero category wants `pzero = 0` in bf().
  cogmod_betadiscrete = list(mu = 0.5, phi = 1, pzero = 0.05)
)


# The natural-scale value each dpar starts at, or NULL for a family with none.
# The ndt + poutlier families take theirs from their registry entry
# (.mixture_spec()) plus the two the parameterization introduces; the rest
# come from .INIT_PLAIN.
#' @keywords internal
.init_targets <- function(family) {
  fam <- .family_name(family)
  if (is.null(fam)) {
    return(NULL)
  }
  if (fam %in% .OUTLIER_FAMILIES) {
    if (is.character(family)) {
      stop(
        "cogmod_inits() needs the family object rather than its name, ",
        "because the links ride on it. Build the formula with ",
        "bf(..., family = ",
        fam,
        "()).",
        call. = FALSE
      )
    }
    # In seconds. The `ndt` entry is a fallback: cogmod_inits() replaces it
    # with .ndt_start() of the response once it has the data, and only a
    # response with no finite positive value leaves it in place.
    return(c(.mixture_spec(fam)$init, list(ndt = 0.1, poutlier = 0.02)))
  }
  if (!fam %in% names(.INIT_PLAIN)) {
    return(NULL)
  }
  .INIT_PLAIN[[fam]]
}


# Turning declarations into values ----------------------------------------

# One entry per declared parameter: the value on its declared scale, its shape,
# its bounds, and how it may be jittered.
#
# A declaration whose dimensions cannot be resolved is dropped rather than
# guessed at. That costs only the CmdStan message about a partial init list,
# whereas a wrong dimension is a hard error at startup.
#' @keywords internal
.init_plan <- function(decls, sdata, targets, links) {
  plan <- lapply(decls, function(d) {
    dims <- vapply(d$dims, .eval_stan_expr, numeric(1), sdata = sdata)
    if (anyNA(dims)) {
      return(NULL)
    }
    dims <- as.integer(dims)
    n <- if (length(dims)) prod(dims) else 1L
    bounds <- .stan_bounds(d$bounds, sdata)
    lower <- bounds[["lower"]]
    upper <- bounds[["upper"]]

    entry <- list(
      name = d$name,
      dim = dims,
      lower = lower,
      upper = upper,
      kind = "bounds",
      tier = .init_tier(d$name)
    )

    # Structured types: a valid value is not just a number in a range, so the
    # bounds-based default and the jitter are both skipped.
    structured <- switch(
      d$type,
      cholesky_factor_corr = ,
      cholesky_factor_cov = ,
      corr_matrix = ,
      cov_matrix = as.vector(diag(dims[1])),
      simplex = rep(1 / dims[1], dims[1]),
      unit_vector = c(1, rep(0, dims[1] - 1)),
      NULL
    )
    if (!is.null(structured)) {
      entry$value <- structured
      entry$kind <- "fixed"
      if (
        d$type %in%
          c(
            "cholesky_factor_corr",
            "cholesky_factor_cov",
            "corr_matrix",
            "cov_matrix"
          )
      ) {
        entry$dim <- c(dims[1], dims[1])
      }
      return(entry)
    }
    # An ordered vector has to stay ordered, so it is jittered and re-sorted.
    if (d$type %in% c("ordered", "positive_ordered")) {
      base <- if (identical(d$type, "ordered")) 0 else 0.5
      entry$value <- base + seq_len(n) * 0.5
      entry$kind <- "sorted"
      return(entry)
    }

    entry$value <- .init_value(d, n, lower, upper, targets, links, sdata)
    entry$scale <- .jitter_scale(d, n, sdata)
    entry
  })
  plan[!vapply(plan, is.null, logical(1))]
}


# How far the jitter on a block of regression coefficients may move the
# linear predictor. The noise on a coefficient is multiplied by its design
# column, so at the row where the columns reach furthest the block moves the
# predictor with SD jitter * ||X[i, ]||. The block's jitter is divided by that
# largest row norm wherever it exceeds 1, which keeps the SD at or under the
# tier's value at every row and leaves a dummy column or a standardized
# covariate alone. Dividing by the column SD instead would not do: the columns
# that bit - a t2()'s unpenalised `bs_*`, age in years - have SDs under 1, and
# it is their extreme rows that put responses below `ndt` (NEWS, 0.3.4).
#
# `b` is matched to the matrix brms builds from `X`: centred and without its
# Intercept column when the formula has an intercept (the parameter is then
# one shorter than `X`), as it stands otherwise, where an `Intercept` column
# under `0 + Intercept` keeps the full jitter. `bs` is matched to `Xs`, which
# mgcv has already centred. Returns one factor per element, or NULL for
# anything else, which keeps the tier's jitter as it is.
#' @keywords internal
.jitter_scale <- function(d, n, sdata) {
  if (!identical(d$type, "vector") || length(d$dims) != 1) {
    return(NULL)
  }
  nm <- d$name
  design <- if (identical(nm, "b")) {
    "X"
  } else if (startsWith(nm, "b_")) {
    sub("^b_", "X_", nm)
  } else if (identical(nm, "bs")) {
    "Xs"
  } else if (startsWith(nm, "bs_")) {
    sub("^bs_", "Xs_", nm)
  } else {
    return(NULL)
  }
  X <- sdata[[design]]
  if (is.null(X)) {
    return(NULL)
  }
  X <- as.matrix(X)
  icpt <- if (is.null(colnames(X))) {
    rep(FALSE, ncol(X))
  } else {
    colnames(X) %in% "Intercept"
  }
  out <- rep(1, n)
  if (any(icpt) && n == ncol(X) - sum(icpt)) {
    slope <- rep(TRUE, n)
    Z <- X[, !icpt, drop = FALSE]
    Z <- sweep(Z, 2, colMeans(Z))
  } else if (n == ncol(X)) {
    slope <- !icpt
    Z <- X[, slope, drop = FALSE]
  } else {
    return(NULL)
  }
  if (!any(slope) || !nrow(Z)) {
    return(NULL)
  }
  reach <- max(sqrt(rowSums(Z^2)))
  if (is.finite(reach) && reach > 1) {
    out[slope] <- 1 / reach
  }
  out
}


# The value for a plain real / vector / matrix parameter, preferring a
# family-specific target where one applies.
#' @keywords internal
.init_value <- function(d, n, lower, upper, targets, links, sdata) {
  scalar <- identical(d$type, "real") && !length(d$dims)

  # brms calls the main parameter's intercept `Intercept`, and every other
  # dpar's `Intercept_<dpar>`. A dpar omitted from the formula becomes a bare
  # auxiliary parameter named after itself, on the natural scale. Only a scalar
  # `real` is ever a dpar: brms also declares vectors named after its own
  # machinery (`vector[Ks] bs` for a smooth, for one), and matching those by
  # name would hand a spline coefficient a dpar's starting value.
  if (scalar) {
    dpar <- if (identical(d$name, "Intercept")) {
      "mu"
    } else if (startsWith(d$name, "Intercept_")) {
      sub("^Intercept_", "", d$name)
    } else {
      NA_character_
    }
    if (!is.na(dpar) && dpar %in% names(targets)) {
      v <- .apply_link(targets[[dpar]], links[dpar])
      if (!is.na(v)) return(v)
    }
    if (d$name %in% names(targets)) {
      return(targets[[d$name]]) # natural scale, no link applied
    }
  }

  # Under `0 + Intercept` brms declares no separate intercept: the column named
  # `Intercept` inside the (uncentered) design matrix carries it, so the target
  # goes to that element of `b`. The centered parameterization drops the first
  # column instead, which is how the two are told apart - there `b` is one
  # shorter than the design matrix and every element really is a slope.
  if (
    identical(d$type, "vector") &&
      length(d$dims) == 1 &&
      (identical(d$name, "b") || startsWith(d$name, "b_"))
  ) {
    dpar <- if (identical(d$name, "b")) "mu" else sub("^b_", "", d$name)
    out <- rep(0, n)
    cn <- colnames(sdata[[
      if (identical(dpar, "mu")) "X" else paste0("X_", dpar)
    ]])
    if (length(cn) == n && dpar %in% names(targets) && "Intercept" %in% cn) {
      v <- .apply_link(targets[[dpar]], links[dpar])
      if (!is.na(v)) out[match("Intercept", cn)] <- v
    }
    return(out)
  }

  # A smooth's SD starts near flat, at 0.1 link units of wiggle: 0.1 /
  # .zs_scales() in units of `sds`, about 1 on a t2(..., bs = c("cr", "cr"))
  # and 0.017 on s(x, bs = "cr"). The generic lower-bound value, 0.25, is a
  # plausible group-level SD but a lot of wiggle for a spline. With `zs` at
  # zero under the hierarchical jitter, the starting smooth stays flat. Each
  # penalty gets its own scale; 0.05 is the fallback for a block whose scale
  # cannot be read.
  if (startsWith(d$name, "sds_")) {
    sc <- .zs_scales(sdata, sub("^sds_", "", d$name))
    if (length(sc) == n && !anyNA(sc)) {
      return(unname(0.1 / sc))
    }
    return(rep(0.05, n))
  }

  rep(.default_value(lower, upper), n)
}


# Where no target applies: zero for a free parameter (a slope of zero and a
# standardized group-level effect of zero are both sensible), just inside the
# bound for a constrained one, the midpoint when bounded on both sides.
#' @keywords internal
.default_value <- function(lower, upper) {
  if (is.na(lower) && is.na(upper)) {
    return(0)
  }
  if (is.na(upper)) {
    return(lower + 0.25)
  }
  if (is.na(lower)) {
    return(upper - 0.25)
  }
  lower + (upper - lower) / 2
}


# Jitter on the scale Stan samples on, so the result cannot leave the bounds:
# additive when free, multiplicative on the distance to a single bound, and on
# the logit scale between two.
#' @keywords internal
.jitter_bounded <- function(v, lower, upper, jitter) {
  e <- stats::rnorm(length(v), 0, jitter)
  if (is.na(lower) && is.na(upper)) {
    return(v + e)
  }
  if (is.na(upper)) {
    return(lower + pmax(v - lower, 1e-8) * exp(e))
  }
  if (is.na(lower)) {
    return(upper - pmax(upper - v, 1e-8) * exp(e))
  }
  w <- pmin(pmax((v - lower) / (upper - lower), 1e-8), 1 - 1e-8)
  lower + (upper - lower) * stats::plogis(stats::qlogis(w) + e)
}


# Which of the two jitters a declared parameter gets. brms names the
# standardized group-level effects `z_<k>`, the smooth coefficients
# `zs_<k>_<j>` and their scales `sd_<k>` and `sds_<k>_<j>`; Gaussian processes
# follow with `zgp_`, `sdgp_` and `lscale_`. A unit of noise on one of those is
# multiplied by a scale and a design column, and every group and basis
# function draws its own, so the predictor at one row moves by the sum of
# many: at the full 0.25, a tensor-smooth model's starts moved by one to two
# and a half link units (NEWS, 0.3.3). Intercepts and population-level slopes
# keep the full jitter, which is where between-chain dispersion makes Rhat
# mean something.
#' @keywords internal
.init_tier <- function(name) {
  if (grepl("^(z|zs|sd|sds|zgp|sdgp|lscale)_", name)) {
    "hierarchical"
  } else {
    "population"
  }
}


# `jitter` as the user gives it - one SD for everything, or one per tier - as
# the named pair .init_fun() reads. Idempotent, so it can be applied at the
# door and again inside.
#' @keywords internal
.init_jitter <- function(jitter) {
  if (is.null(jitter)) {
    jitter <- 0.25
  }
  if (
    !is.numeric(jitter) ||
      !length(jitter) %in% 1:2 ||
      anyNA(jitter) ||
      any(jitter < 0)
  ) {
    stop(
      "`jitter` must be one non-negative number, or two: population-level ",
      "then hierarchical.",
      call. = FALSE
    )
  }
  if (length(jitter) == 1) {
    jitter <- c(jitter, jitter / 5)
  }
  c(population = jitter[[1]], hierarchical = jitter[[2]])
}


# Reading the Stan program ------------------------------------------------

# Every declaration in the `parameters` block, as name / type / dimension
# expressions / bounds expression. Dimensions and bounds are left as text
# because they are Stan expressions over the data block (`Kc`, `knots_1[1]`),
# and are resolved against make_standata() output later.
#' @keywords internal
.stan_param_decls <- function(code) {
  lines <- strsplit(code, "\n")[[1]]
  start <- grep("^\\s*parameters\\s*\\{", lines)
  if (!length(start)) {
    return(list())
  }
  close <- grep("^\\s*\\}\\s*$", lines)
  close <- close[close > start[1]]
  if (!length(close)) {
    return(list())
  }
  block <- sub("//.*$", "", lines[(start[1] + 1):(close[1] - 1)])

  stmts <- trimws(strsplit(paste(block, collapse = " "), ";", fixed = TRUE)[[
    1
  ]])
  out <- lapply(stmts[nzchar(stmts)], .parse_stan_decl)
  out[!vapply(out, is.null, logical(1))]
}


#' @keywords internal
.STAN_PARAM_TYPES <- c(
  "real",
  "vector",
  "row_vector",
  "matrix",
  "simplex",
  "ordered",
  "positive_ordered",
  "unit_vector",
  "cholesky_factor_corr",
  "cholesky_factor_cov",
  "corr_matrix",
  "cov_matrix"
)


# One declaration. Hand-rolled rather than a single regex because dimensions
# nest (`vector[knots_1[1]] zs_1_1`), which a regex cannot bracket-match.
# Anything unrecognised returns NULL and is skipped.
#' @keywords internal
.parse_stan_decl <- function(s) {
  s <- trimws(s)
  dims <- character(0)

  if (grepl("^array\\s*\\[", s)) {
    g <- .take_bracket(sub("^array\\s*", "", s))
    if (is.null(g)) {
      return(NULL)
    }
    dims <- .split_top_commas(g$inside)
    s <- trimws(g$rest)
  }

  type <- regmatches(s, regexpr("^[A-Za-z_][A-Za-z0-9_]*", s))
  if (!length(type) || !type %in% .STAN_PARAM_TYPES) {
    return(NULL)
  }
  s <- trimws(substring(s, nchar(type) + 1L))

  bounds <- ""
  if (startsWith(s, "<")) {
    i <- regexpr(">", s, fixed = TRUE)
    if (i < 0) {
      return(NULL)
    }
    bounds <- substring(s, 2L, i - 1L)
    s <- trimws(substring(s, i + 1L))
  }

  if (startsWith(s, "[")) {
    g <- .take_bracket(s)
    if (is.null(g)) {
      return(NULL)
    }
    dims <- c(dims, .split_top_commas(g$inside))
    s <- trimws(g$rest)
  }

  if (!grepl("^[A-Za-z_][A-Za-z0-9_]*$", s)) {
    return(NULL)
  }
  list(name = s, type = type, dims = dims, bounds = bounds)
}


# Contents of the bracketed group `s` opens with, plus whatever follows it.
#' @keywords internal
.take_bracket <- function(s) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(chars) || chars[1] != "[") {
    return(NULL)
  }
  depth <- 0L
  for (i in seq_along(chars)) {
    if (chars[i] == "[") {
      depth <- depth + 1L
    }
    if (chars[i] == "]") {
      depth <- depth - 1L
      if (depth == 0L) {
        return(list(
          inside = paste(chars[seq_len(i - 1L)[-1]], collapse = ""),
          rest = paste(chars[-seq_len(i)], collapse = "")
        ))
      }
    }
  }
  NULL
}


# Split on the commas that separate arguments, ignoring those inside brackets
# or parentheses.
#' @keywords internal
.split_top_commas <- function(s) {
  chars <- strsplit(s, "", fixed = TRUE)[[1]]
  if (!length(chars)) {
    return(character(0))
  }
  depth <- 0L
  grp <- 1L
  which_grp <- integer(length(chars))
  for (i in seq_along(chars)) {
    if (chars[i] %in% c("[", "(")) {
      depth <- depth + 1L
    }
    if (chars[i] %in% c("]", ")")) {
      depth <- depth - 1L
    }
    if (chars[i] == "," && depth == 0L) {
      grp <- grp + 1L
      next # which_grp stays 0, so the separator itself is dropped
    }
    which_grp[i] <- grp
  }
  keep <- which_grp > 0L
  parts <- split(chars[keep], factor(which_grp[keep], levels = seq_len(grp)))
  trimws(unname(vapply(parts, paste, character(1), collapse = "")))
}


# `lower` and `upper` out of a `<...>` bounds expression, as numbers.
#' @keywords internal
.stan_bounds <- function(bounds, sdata) {
  out <- c(lower = NA_real_, upper = NA_real_)
  if (!nzchar(bounds)) {
    return(out)
  }
  for (part in .split_top_commas(bounds)) {
    kv <- regmatches(part, regexec("^\\s*(lower|upper)\\s*=\\s*(.+)$", part))[[
      1
    ]]
    if (length(kv) == 3) out[[kv[2]]] <- .eval_stan_expr(kv[3], sdata)
  }
  out
}


# A Stan expression over the data block, as a number. `enclos = baseenv()` so
# that arithmetic resolves but nothing from the caller's workspace does.
#' @keywords internal
.eval_stan_expr <- function(expr, sdata) {
  v <- tryCatch(
    eval(parse(text = expr), envir = sdata, enclos = baseenv()),
    error = function(e) NULL
  )
  if (!is.numeric(v) || length(v) != 1 || is.na(v)) {
    return(NA_real_)
  }
  as.numeric(v)
}


# Natural scale -> link scale. Returns NA for a link with no inverse here,
# which leaves the parameter on its generic default rather than erroring: an
# unhelpful starting value is a smaller problem than a function that refuses to
# run.
#' @keywords internal
.apply_link <- function(x, link) {
  if (is.null(link) || is.na(link)) {
    return(NA_real_)
  }
  out <- switch(
    link,
    identity = x,
    log = log(x),
    logm1 = log(x - 1),
    log1p = log1p(x),
    inverse = 1 / x,
    sqrt = sqrt(x),
    logit = stats::qlogis(x),
    probit = ,
    probit_approx = stats::qnorm(x),
    cloglog = log(-log1p(-x)),
    cauchit = stats::qcauchy(x),
    softplus = log(expm1(x)),
    # squareplus(z) = (z + sqrt(z^2 + 4)) / 2, whose inverse is x - 1 / x
    squareplus = x - 1 / x,
    tan_half = tan(x / 2),
    NA_real_
  )
  if (!is.numeric(out) || length(out) != 1 || !is.finite(out)) NA_real_ else out
}
