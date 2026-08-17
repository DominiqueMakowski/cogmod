#' @rdname rrt_lognormal
#' @param link_mu,link_sigma,link_ndt,link_poutlier Link functions for the
#'   parameters.
#' @param predict_outliers Logical; whether `posterior_predict()` and
#'   `posterior_epred()` should include the outlier component. `FALSE` (the
#'   default in `rt_lognormal()`) fixes `poutlier` to zero for prediction, so
#'   predictions describe the decision process alone; the likelihood is always
#'   the full mixture either way. On the prediction methods themselves the
#'   default is `NULL`, which defers to the flag carried on the model - see
#'   [with_outliers()] to change it after fitting. See Details.
#' @export
rt_lognormal <- function(
  link_mu = "identity",
  link_sigma = "softplus",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE
) {
  fam <- brms::custom_family(
    name = "rt_lognormal",
    dpars = c("mu", "sigma", "ndt", "poutlier"),
    links = c(link_mu, link_sigma, link_ndt, link_poutlier),
    lb = c(NA, 0, 0, 0), # sigma > 0, ndt > 0, poutlier >= 0
    ub = c(NA, NA, NA, 1), # poutlier <= 1
    type = "real" # Continuous outcome variable (RT)
  )
  # The flag rides on the family because that is the only thing brms carries
  # down to a custom family's prediction methods. See the Details section of
  # ?rrt_lognormal, and with_outliers() for flipping it after fitting.
  fam$predict_outliers <- isTRUE(predict_outliers)
  fam
}


#' @keywords internal
.rt_lognormal_lpdf <- function() {
  "
// Log-likelihood for one observation from the shifted LogNormal distribution.
// Y: observed reaction time (seconds).
// mu: mean of the decision time on the log scale (meanlog).
// sigma: SD of the decision time on the log scale (> 0).
// ndt: non-decision time in seconds (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is fixed at a half Student-t with scale 0.4 and 3
// degrees of freedom. Its role is to keep the density strictly positive below
// `ndt`, where the shifted decision component has none. That is what removes
// the hard min-RT boundary and lets `ndt` be estimated directly instead of as a
// fraction of an observed minimum. The half-t is flat at the origin, so the
// fastest responses are not starved of density, and its tails are heavy enough
// to cover the whole plausible RT range.
real rt_lognormal_lpdf(real Y, real mu, real sigma, real ndt, real poutlier) {
    // Parameter checks
    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();

    // log(2) folds the symmetric Student-t onto [0, Inf).
    real lp_out = log(2) + student_t_lpdf(Y | 3, 0, 0.4);
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response. Finite, and smooth in ndt, because the LogNormal
    // decision density approaches zero with vanishing derivatives at t_adj = 0.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    return log_mix(poutlier, lp_out, lognormal_lpdf(t_adj | mu, sigma));
    }
"
}


#' @rdname rrt_lognormal
#' @export
rt_lognormal_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
    "functions {
",
    .rt_lognormal_lpdf(),
    "}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_lognormal_lpdf
}


#' @rdname rrt_lognormal
#' @export
rt_lognormal_stanvars <- function() {
  brms::stanvar(scode = .rt_lognormal_lpdf(), block = "functions")
}



# brms methods ------------------------------------------------------------

#' @rdname rrt_lognormal
#' @inheritParams lnr
#' @export
log_lik_rt_lognormal <- function(i, prep) {
  # Extract observation
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)

  n_draws <- max(length(mu), length(sigma), length(ndt), length(poutlier))
  if (n_draws == 0) return(numeric(0))

  # Compute log-likelihood using the vectorized drt_lognormal function
  ll <- drt_lognormal(
    x = rep(y, length.out = n_draws),
    mu = mu,
    sigma = sigma,
    ndt = ndt,
    poutlier = poutlier,
    log = TRUE
  )

  ll[is.na(ll)] <- -Inf
  ll
}


# Resolve whether predictions should include the outlier component. An explicit
# argument wins; otherwise the flag carried on the family is used (set by
# rt_lognormal() or with_outliers()), defaulting to excluding it.
#' @keywords internal
.predict_outliers <- function(predict_outliers, prep) {
  if (!is.null(predict_outliers)) {
    return(isTRUE(predict_outliers))
  }
  isTRUE(prep$family$predict_outliers)
}


#' @rdname rrt_lognormal
#' @export
posterior_predict_rt_lognormal <- function(
  i,
  prep,
  predict_outliers = NULL,
  ...
) {
  # Get parameters for observation i across all draws
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }

  n_draws <- max(length(mu), length(sigma), length(ndt))

  # Simulate using rrt_lognormal (vectorized over its parameters)
  out <- rrt_lognormal(
    n = n_draws,
    mu = mu,
    sigma = sigma,
    ndt = ndt,
    poutlier = poutlier
  )

  # Return as a matrix (draws x 1)
  as.matrix(out)
}


#' @rdname rrt_lognormal
#' @export
posterior_epred_rt_lognormal <- function(prep, predict_outliers = NULL) {
  # Extract draws (matrices: draws x observations)
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  ndt <- brms::get_dpar(prep, "ndt")

  # Expected RT of the decision process alone
  mean_dec <- exp(mu + sigma^2 / 2) + ndt
  if (!.predict_outliers(predict_outliers, prep)) {
    return(mean_dec)
  }

  # E[RT] of the mixture
  poutlier <- brms::get_dpar(prep, "poutlier")
  (1 - poutlier) * mean_dec + poutlier * .mcontam()
}


#' Include or exclude the outlier component in predictions
#'
#' @description
#' Switches the `predict_outliers` flag on a model fitted with [rt_lognormal()],
#' controlling whether `posterior_predict()` and `posterior_epred()` describe the
#' fitted mixture or the decision process alone.
#'
#' Predictions **exclude** the outlier component by default, because for almost
#' every downstream use it is a nuisance: it pulls expected values toward its own
#' mean (0.441 s) and adds a spike of implausibly fast draws to posterior
#' predictive samples. It is also a deliberately fixed regularizer rather than a
#' claim about how guesses are distributed, so simulating from it means
#' simulating from something the model does not assert.
#'
#' `with_outliers()` restores the mixture. The main reason to want it is
#' `brms::pp_check()`: on untrimmed data the decision-only predictive has no fast
#' spike to match the one in the data, which reads as misfit. Use
#' `pp_check(with_outliers(m))` for a like-for-like check.
#'
#' The flag is stored **on the model** rather than passed as an argument, because
#' `brms` and the packages built on it (`insight`, `modelbased`,
#' `marginaleffects`, `emmeans`) do not forward extra arguments down to a custom
#' family's prediction methods - `posterior_epred()` reaches the family method
#' with `prep` and nothing else. Carrying it on the object is what makes it work
#' through all of them. The same flag can be set up front with
#' `rt_lognormal(predict_outliers = TRUE)`.
#'
#' `log_lik()` is unaffected and has no equivalent switch: the likelihood *is*
#' the mixture, and dropping a component from it would not be a different summary
#' of the same model but a different model. One consequence worth knowing is that
#' `posterior_predict()` and `log_lik()` no longer describe the same
#' distribution by default, so a hand-rolled LOO-PIT check should use
#' `with_outliers()`.
#'
#' @param object A `brmsfit` fitted with the [rt_lognormal()] family.
#'
#' @return The model, with the flag set. The fit itself is untouched - only how
#'   predictions are summarised changes.
#'
#' @examples
#' # m <- brms::brm(bf(RT ~ Condition, ndt ~ 1, poutlier ~ 1,
#' #                   family = rt_lognormal()),
#' #                data = df, stanvars = rt_lognormal_stanvars())
#' #
#' # # the decision process alone - the default, everywhere downstream
#' # brms::posterior_epred(m)
#' # modelbased::estimate_means(m, by = "Condition")
#' # marginaleffects::avg_predictions(m, by = "Condition")
#' #
#' # # the fitted mixture, e.g. for a like-for-like predictive check
#' # brms::pp_check(with_outliers(m))
#' # without_outliers(m)  # back to the default
#'
#' @export
with_outliers <- function(object) {
  .set_predict_outliers(object, TRUE)
}


#' @rdname with_outliers
#' @export
without_outliers <- function(object) {
  .set_predict_outliers(object, FALSE)
}


#' @keywords internal
.set_predict_outliers <- function(object, value) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.", call. = FALSE)
  }
  if (!identical(object$family$name, "rt_lognormal")) {
    stop("`object` must have been fitted with the rt_lognormal() family.",
         call. = FALSE)
  }
  # brms rebuilds prep$family from object$formula$family, so the flag has to
  # live there to survive prepare_predictions(). object$family is set too, so
  # that inspecting the fit tells the same story.
  object$formula$family$predict_outliers <- value
  object$family$predict_outliers <- value
  object
}


#' Per-trial outlier probabilities
#'
#' @description
#' Posterior probability that each response was generated by the outlier
#' component rather than by the decision process, for a model fitted with
#' [rt_lognormal()]. This is the mixture *responsibility*
#' `poutlier * g(rt) / (poutlier * g(rt) + (1 - poutlier) * f(rt - ndt))`,
#' averaged over posterior draws.
#'
#' Responses faster than `ndt` come out at 1, responses in the heart of the
#' distribution near 0, and responses in either tail somewhere in between - the
#' model discriminates by evidence rather than by a cutoff. A response in the
#' middle can still be an outlier; a low probability means the data cannot tell,
#' not that the trial is clean.
#'
#' @param object A `brmsfit` fitted with the [rt_lognormal()] family.
#' @param summary Logical; if `TRUE` (default) returns a data frame with one row
#'   per observation. If `FALSE`, returns the full draws x observations matrix.
#'
#' @return A data frame with columns `rt`, `p_outlier`, and `fast` (whether the
#'   response falls below the median RT), or a matrix if `summary = FALSE`.
#'
#' @examples
#' # m <- brms::brm(bf(RT ~ 1, ndt ~ 1, poutlier ~ 1, family = rt_lognormal()),
#' #                data = df, stanvars = rt_lognormal_stanvars())
#' # head(p_outlier(m))
#'
#' @export
p_outlier <- function(object, summary = TRUE) {
  if (!inherits(object, "brmsfit")) {
    stop("`object` must be a brmsfit.")
  }
  prep <- brms::prepare_predictions(object)
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in the fitted model.")
  }
  y <- prep$data$Y
  n <- length(y)

  mu <- brms::get_dpar(prep, "mu")
  n_draws <- if (is.matrix(mu)) nrow(mu) else brms::ndraws(object)
  as_mat <- function(x) {
    if (is.matrix(x)) x else matrix(x, nrow = n_draws, ncol = n)
  }

  mu <- as_mat(mu)
  sigma <- as_mat(brms::get_dpar(prep, "sigma"))
  ndt <- as_mat(brms::get_dpar(prep, "ndt"))
  poutlier <- as_mat(brms::get_dpar(prep, "poutlier"))
  Y <- matrix(y, nrow = n_draws, ncol = n, byrow = TRUE)

  dens_out <- .dcontam(Y)
  x_adj <- Y - ndt
  dens_dec <- ifelse(
    x_adj > 0,
    stats::dlnorm(pmax(x_adj, 1e-300), meanlog = mu, sdlog = sigma),
    0
  )

  num <- poutlier * dens_out
  r <- num / (num + (1 - poutlier) * dens_dec)
  r[!is.finite(r)] <- 1 # both components vanish: attribute to the outlier

  if (!summary) return(r)

  data.frame(
    rt = y,
    p_outlier = colMeans(r),
    fast = y < stats::median(y)
  )
}
