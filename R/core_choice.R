# Shared machinery for the shifted choice + RT families
# =====================================================
#
# The counterpart of core_shifted.R for the models that produce a *choice* as
# well as a reaction time. The structure is the same - a decision-time
# distribution shifted by `ndt`, mixed with a half Student-t outlier component
# of weight `poutlier` and scale `minrt` - with one genuine difference, which is
# why these cannot simply be added to `.SHIFTED`.
#
# A choice model has a *defective* density per response option `k`, with
# `sum_k integral f_k(t) dt = 1`. A contaminant process therefore has to produce
# a choice as well as an RT. The construction used here is a guess: the choice
# is uniform over the `K` options and the RT comes from the same half Student-t
# as in the RT-only families.
#
#   f(t, k) = poutlier * (1 / K) * g(t) + (1 - poutlier) * f_k(t - ndt)
#
# The `1 / K` is what keeps this normalized:
#
#   integral sum_k [ p*(1/K)*g(t) + (1-p)*f_k(t-ndt) ] dt
#     = p * sum_k (1/K) * integral g  +  (1-p) * sum_k integral f_k
#     = p * 1 + (1-p) * 1 = 1
#
# Without it the total comes to `1 + p` rather than `1`, which is silent: every
# individual density is still positive and still looks reasonable.
#
# See ?rcogmod_lognormal for the account of why `ndt` is expressed directly
# rather than as a fraction of an observed minimum, what the outlier component
# is for, and why `minrt` is a constant on the family rather than a dpar. The
# outlier component itself (`.dcontam()`, `.validate_minrt()`, `.minrt()`,
# `.predict_outliers()`, `with_outliers()`, `p_outlier()`) is shared with the
# RT-only families and lives in core_shifted.R.


# The registry ------------------------------------------------------------

# Each entry describes the *decision* component only.
#
#  dpars/links/lb/ub : the distributional parameters, before ndt and poutlier
#  K                 : number of response options
#  vars              : the brms `vars` the family needs for the response index
#  stan_check        : Stan expression that is TRUE for invalid parameters
#  stan_dens         : Stan expression for the decision log-density at `t_adj`,
#                      which may branch on `dec`
#  ldens             : R log-density of the decision component, at t > 0, given
#                      the response index `k`
#  rng               : R sampler; returns data.frame(rt, response), unshifted
#  init              : natural-scale starting values, for cogmod_inits()
#  label             : human-readable name, used in the generated Stan comments
#' @keywords internal
.CHOICE <- list(
  cogmod_lnr = list(
    # `mu` is nuzero: brms requires the first dpar of a custom family to be
    # called `mu`. Both nu are minus the log-scale mean, so larger = faster.
    dpars = c("mu", "nuone", "sigmazero", "sigmaone"),
    links = c("identity", "identity", "softplus", "softplus"),
    lb = c(NA, NA, 0, 0), ub = c(NA, NA, NA, NA),
    K = 2L,
    vars = "dec[n]",
    stan_check = "sigmazero <= 0 || sigmaone <= 0",
    stan_dens = paste0(
      "dec == 0\n",
      "      ? lognormal_lpdf(t_adj | -mu, sigmazero)",
      " + lognormal_lccdf(t_adj | -nuone, sigmaone)\n",
      "      : lognormal_lpdf(t_adj | -nuone, sigmaone)",
      " + lognormal_lccdf(t_adj | -mu, sigmazero)"
    ),
    ldens = function(t, k, p) {
      win_ml <- ifelse(k == 0, -p$mu, -p$nuone)
      win_sd <- ifelse(k == 0, p$sigmazero, p$sigmaone)
      los_ml <- ifelse(k == 0, -p$nuone, -p$mu)
      los_sd <- ifelse(k == 0, p$sigmaone, p$sigmazero)
      # lower.tail = FALSE gives the log-survival of the loser directly, which
      # is stable where log(1 - exp(log F)) is not.
      stats::dlnorm(t, meanlog = win_ml, sdlog = win_sd, log = TRUE) +
        stats::plnorm(t, meanlog = los_ml, sdlog = los_sd,
                      lower.tail = FALSE, log.p = TRUE)
    },
    rng = function(n, p) {
      d0 <- stats::rlnorm(n, meanlog = -p$mu, sdlog = p$sigmazero)
      d1 <- stats::rlnorm(n, meanlog = -p$nuone, sdlog = p$sigmaone)
      data.frame(rt = pmin(d0, d1), response = as.numeric(d0 >= d1))
    },
    init = list(mu = 0.7, nuone = 0.7, sigmazero = 0.5, sigmaone = 0.5),
    dpar_doc = c(
      "dec: the observed choice, 0 or 1.",
      paste("mu: nuzero, the processing speed of accumulator 0",
            "(meanlog = -mu, so larger is faster)."),
      "nuone: the same for accumulator 1.",
      "sigmazero, sigmaone: the log-scale SD of each accumulator (> 0)."
    ),
    label = "Log-Normal Race"
  )
)


# The choice families built on this machinery. Folded into .OUTLIER_FAMILIES in
# core_shifted.R, so with_outliers(), p_outlier(), cogmod_priors() and
# cogmod_inits() all reach them.
#' @keywords internal
.CHOICE_FAMILIES <- names(.CHOICE)


#' @keywords internal
.choice_spec <- function(name) {
  spec <- .CHOICE[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the shifted choice families. Supported: ",
      paste0(.CHOICE_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Either registry, for the code that is genuinely family-agnostic: p_outlier()
# and cogmod_inits(). The univariate machinery keeps using .shifted_spec(), so
# it cannot be handed a choice family by accident - the entries are not
# interchangeable, since a choice family's `ldens` needs the response index.
#' @keywords internal
.mixture_spec <- function(name) {
  spec <- .SHIFTED[[name]]
  if (is.null(spec)) spec <- .CHOICE[[name]]
  if (is.null(spec)) {
    stop(
      "'", name, "' is not one of the families built on the outlier mixture. ",
      "Supported: ", paste0(.OUTLIER_FAMILIES, "()", collapse = ", "), ".",
      call. = FALSE
    )
  }
  spec
}


# Family constructor ------------------------------------------------------

# Builds the brms custom family for `name`, appending `ndt` and `poutlier` to
# whatever decision dpars the registry lists. Same as .shifted_family() plus the
# `vars` the response index arrives through.
#' @keywords internal
.choice_family <- function(name, links = NULL, predict_outliers = FALSE,
                           minrt = 0.3) {
  spec <- .choice_spec(name)
  dec_links <- if (is.null(links)) spec$links else links
  if (length(dec_links) != length(spec$dpars)) {
    stop("Expected ", length(spec$dpars), " link functions for ", name, "().",
         call. = FALSE)
  }
  fam <- brms::custom_family(
    name = name,
    dpars = c(spec$dpars, "ndt", "poutlier"),
    links = c(dec_links, "log", "logit"),
    lb = c(spec$lb, 0, 0),
    ub = c(spec$ub, NA, 1),
    type = "real",
    vars = spec$vars
  )
  # Both ride on the family because that is the only thing brms carries down to
  # a custom family's prediction methods, and because a dpar left out of the
  # formula would be estimated rather than defaulted. See ?rcogmod_lognormal.
  fam$predict_outliers <- isTRUE(predict_outliers)
  invisible(.validate_minrt(minrt))
  fam$minrt <- minrt
  fam
}


# Stan code ---------------------------------------------------------------

# Generates the `<name>_lpdf` Stan function: the same mixture skeleton for every
# choice family, with only the decision density swapped in. `minrt` is baked in
# as a literal because Stan functions cannot see the data block, and because a
# dpar would be estimated whenever the user left it out of the formula.
#
# brms appends `vars` after the dpars, so the response index is the *last*
# argument and is an `int`.
#' @keywords internal
.choice_lpdf <- function(name, minrt = 0.3, prelude = "") {
  spec <- .choice_spec(name)
  if (!is.null(spec$prelude)) prelude <- get(spec$prelude)
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(.validate_minrt(minrt), format = "g", digits = 15, width = 1)

  # As in .shifted_lpdf(): the outlier term has no parameter in it, so its
  # normalising constant is folded into a literal rather than left for Stan to
  # recompute on every leapfrog step. The `-log(K)` that makes the guess uniform
  # over the response options goes into the same constant - at K = 2 it cancels
  # the log(2) that folds the Student-t onto [0, Inf) exactly.
  nu <- .POUTLIER_DF
  lc <- log(2) + lgamma((nu + 1) / 2) - lgamma(nu / 2) -
    0.5 * log(nu * pi) - log(.validate_minrt(minrt)) - log(spec$K)
  lp_out <- sprintf(
    "%s - %s * log1p(square(Y / %s) / %s)",
    formatC(lc, format = "g", digits = 17, width = 1),
    formatC((nu + 1) / 2, format = "g", digits = 15, width = 1),
    scale,
    formatC(nu, format = "g", digits = 15, width = 1)
  )
  args <- paste(
    c(sprintf("real %s", c("Y", spec$dpars, "ndt", "poutlier")), "int dec"),
    collapse = ", "
  )
  dpar_doc <- if (is.null(spec$dpar_doc)) "" else {
    paste0(paste0("// ", spec$dpar_doc, collapse = "\n"), "\n")
  }
  note <- if (is.null(spec$note)) "" else {
    paste0("//\n", paste0("// ", spec$note, collapse = "\n"), "\n")
  }
  sprintf(
    "%s
// Log-likelihood for one observation from the shifted %s model.
// Y: observed reaction time.
%s// ndt: non-decision time, same unit as Y (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is a guess: the choice is uniform over the %s response
// options and the RT is a half Student-t with 3 degrees of freedom and scale
// %s (= minrt). It keeps the density strictly positive below `ndt`, where the
// shifted decision component has none. That is what removes the hard min-RT
// boundary and lets `ndt` be estimated directly rather than as a fraction of an
// observed minimum. Because the scale is tied to `minrt`, the likelihood is
// equivariant to the unit Y is measured in.
//
// The 1 / %s is what keeps the total density summing to one over the response
// options; without it it comes to 1 + poutlier. The leading constant below
// carries it, along with the normalising constant of the Student-t and the
// log(2) that folds it onto [0, Inf) - all of them constant, so writing them
// out saves Stan recomputing them for every observation on every leapfrog step.
%sreal %s_lpdf(%s) {
    // Parameter checks
    if (%s || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (dec < 0 || dec > %s) return negative_infinity();
    if (Y <= 0) return negative_infinity();

    real lp_out = %s;
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    real lp_dec = %s;
    return log_mix(poutlier, lp_out, lp_dec);
}
",
    prelude, spec$label, dpar_doc, spec$K, scale, spec$K, note, name, args,
    spec$stan_check, spec$K - 1L, lp_out, spec$stan_dens
  )
}


# Shared expose helper: compiles the generated lpdf and hands it back as an R
# function, for checking the Stan code against the R density.
#' @keywords internal
.choice_expose <- function(name, minrt = 0.3) {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {
",
    .choice_lpdf(name, .as_minrt(minrt)),
    "}"
  )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions[[paste0(name, "_lpdf")]]
}


# Density, RNG ------------------------------------------------------------

# Recycle the decision dpars plus ndt/poutlier and the response to a common
# length. `...` holds the decision dpars, by name.
#' @keywords internal
.prepare_choice <- function(name, x = NULL, n = NULL, response = NULL, ndt,
                            poutlier, ...) {
  spec <- .choice_spec(name)
  dec <- list(...)
  missing <- setdiff(spec$dpars, names(dec))
  if (length(missing)) {
    stop("Missing parameter(s) for ", name, "(): ",
         paste(missing, collapse = ", "), ".", call. = FALSE)
  }
  dec <- dec[spec$dpars]

  # Checked against the registry's own bounds, so the R functions reject exactly
  # what `stan_check` rejects and the two cannot disagree. `lb` is strict
  # throughout the registry, which is why the comparison is `<=`.
  for (j in seq_along(spec$dpars)) {
    d <- spec$dpars[j]
    if (!is.na(spec$lb[j]) && any(dec[[d]] <= spec$lb[j], na.rm = TRUE)) {
      stop("`", d, "` must be greater than ", spec$lb[j], ".", call. = FALSE)
    }
    if (!is.na(spec$ub[j]) && any(dec[[d]] >= spec$ub[j], na.rm = TRUE)) {
      stop("`", d, "` must be less than ", spec$ub[j], ".", call. = FALSE)
    }
  }
  if (any(ndt < 0, na.rm = TRUE)) stop("`ndt` must be non-negative.")
  if (any(poutlier < 0 | poutlier > 1, na.rm = TRUE)) {
    stop("`poutlier` must be in [0, 1].")
  }
  # The parentheses are load-bearing: `%in%` binds tighter than `-`, so
  # `response %in% seq_len(K) - 1L` compares against 1:K and then subtracts
  # from the logical.
  options <- seq_len(spec$K) - 1L
  if (!is.null(response)) {
    if (any(!(response %in% options), na.rm = TRUE)) {
      stop("`response` must be one of ", paste(options, collapse = ", "), ".",
           call. = FALSE)
    }
  }

  lens <- c(vapply(dec, length, integer(1)), length(ndt), length(poutlier),
            length(response))
  if (!is.null(x)) {
    if (is.null(response)) {
      stop("`response` must be provided alongside `x`.", call. = FALSE)
    }
    m <- max(length(x), lens)
    if (m == 0) stop("At least one input vector must have non-zero length.")
  } else if (!is.null(n)) {
    if (length(n) > 1) n <- length(n)
    if (length(n) != 1 || n < 0 || n != floor(n)) {
      stop("n must be a single non-negative integer.")
    }
    m <- max(n, lens)
    if (n == 0) m <- 0
  } else {
    stop("Either 'x' or 'n' must be provided.")
  }

  params <- lapply(dec, rep_len, m)
  params$ndt <- rep_len(ndt, m)
  params$poutlier <- rep_len(poutlier, m)
  if (!is.null(x)) params$x <- rep_len(x, m)
  if (!is.null(response)) params$response <- rep_len(response, m)
  params$ndraws <- m
  params
}


# Log-density of the decision component alone, at t (which may be <= 0, giving
# -Inf) for response `k`. Shape is preserved, so this works on the
# draws x observations matrices p_outlier() builds.
#' @keywords internal
.ldec_choice <- function(name, t, k, p) {
  spec <- .choice_spec(name)
  ok <- is.finite(t) & t > 0
  # pmax() only keeps the density functions from being handed a non-positive
  # number; those entries are overwritten with -Inf immediately afterwards.
  ld <- spec$ldens(pmax(t, 1e-300), k, p)
  ld[!ok] <- -Inf
  ld[is.na(ld)] <- -Inf
  dim(ld) <- dim(t)
  ld
}


# Log-density of the outlier component: the half-t, thinned by the K options the
# guess is spread over.
#' @keywords internal
.lout_choice <- function(x, K, minrt = 0.3) {
  .dcontam(x, minrt = minrt, log = TRUE) - log(K)
}


# Mixture density, shared by every choice family's d*() function.
#' @keywords internal
.dchoice <- function(name, x, response, ndt, poutlier, minrt = 0.3,
                     log = FALSE, ...) {
  spec <- .choice_spec(name)
  params <- tryCatch(
    .prepare_choice(name, x = x, response = response, ndt = ndt,
                    poutlier = poutlier, ...),
    error = function(e) {
      warning(conditionMessage(e), ". Returning 0 density / -Inf log-density.")
      list(ndraws = length(x), error = TRUE)
    }
  )
  if (!is.null(params$error)) {
    return(rep(ifelse(log, -Inf, 0), params$ndraws))
  }

  lp_out <- .lout_choice(params$x, spec$K, minrt = minrt)
  lp_dec <- .ldec_choice(name, params$x - params$ndt, params$response, params)

  ld <- .log_mix(params$poutlier, lp_out, lp_dec)
  ld[is.na(ld)] <- -Inf
  if (log) ld else exp(ld)
}


# Mixture RNG, shared by every choice family's r*() function. Returns a data
# frame with `rt` and `response`, so an outlier draw carries a choice as well as
# a time - drawn uniformly over the K options, matching the density.
#' @keywords internal
.rchoice <- function(name, n, ndt, poutlier, minrt = 0.3, ...) {
  spec <- .choice_spec(name)
  params <- .prepare_choice(name, n = n, ndt = ndt, poutlier = poutlier, ...)
  m <- params$ndraws

  rt <- numeric(m)
  response <- numeric(m)
  if (m == 0) return(data.frame(rt = rt, response = response))

  is_out <- stats::runif(m) < params$poutlier

  if (any(is_out)) {
    n_out <- sum(is_out)
    rt[is_out] <- abs(
      .validate_minrt(minrt) * stats::rt(n_out, df = .POUTLIER_DF)
    )
    response[is_out] <- sample(seq_len(spec$K) - 1L, n_out, replace = TRUE)
  }
  if (any(!is_out)) {
    keep <- !is_out
    p <- lapply(params[spec$dpars], function(v) v[keep])
    sim <- spec$rng(sum(keep), p)
    rt[keep] <- sim$rt + params$ndt[keep]
    response[keep] <- sim$response
  }
  data.frame(rt = rt, response = response)
}


# brms methods ------------------------------------------------------------

# Pull the decision dpars for observation i (or all of them) off a brmsprep.
#' @keywords internal
.dpars_from_prep_choice <- function(name, prep, i = NULL) {
  spec <- .choice_spec(name)
  out <- lapply(spec$dpars, function(d) {
    if (is.null(i)) brms::get_dpar(prep, d) else brms::get_dpar(prep, d, i = i)
  })
  stats::setNames(out, spec$dpars)
}


# The response index brms carried through `vars`, for observation i.
#' @keywords internal
.dec_from_prep <- function(prep, i = NULL) {
  if (!"dec" %in% names(prep$data)) {
    stop("Decision variable 'dec' not found in prep$data. Add it to the ",
         "formula with `RT | dec(response) ~ ...`.", call. = FALSE)
  }
  if (is.null(i)) prep$data$dec else prep$data$dec[i]
}


#' @keywords internal
.log_lik_choice <- function(name, i, prep) {
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  dec <- .dpars_from_prep_choice(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)
  response <- .dec_from_prep(prep, i)

  n_draws <- max(vapply(c(dec, list(ndt, poutlier)), length, integer(1)))
  if (n_draws == 0) return(numeric(0))

  ll <- do.call(.dchoice, c(
    list(name = name, x = rep(y, length.out = n_draws),
         response = rep(response, length.out = n_draws), ndt = ndt,
         poutlier = poutlier, minrt = .minrt(prep), log = TRUE),
    dec
  ))
  ll[is.na(ll)] <- -Inf
  ll
}


# Returns a draws x 2 matrix: the RT and the simulated choice. That is the shape
# brms expects from a custom family whose response has more than one column.
#' @keywords internal
.posterior_predict_choice <- function(name, i, prep, predict_outliers = NULL) {
  dec <- .dpars_from_prep_choice(name, prep, i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }
  n_draws <- max(vapply(c(dec, list(ndt)), length, integer(1)))

  out <- do.call(.rchoice, c(
    list(name = name, n = n_draws, ndt = ndt, poutlier = poutlier,
         minrt = .minrt(prep)),
    dec
  ))
  as.matrix(out)
}
