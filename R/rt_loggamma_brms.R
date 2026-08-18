#' @rdname rrt_loggamma
#' @param link_mu,link_sigma,link_shape,link_ndt,link_poutlier Link functions for the
#'   parameters. `shape` is unconstrained and takes an `identity` link, so that the
#'   LogNormal (`shape = 0`) sits in the interior of its range rather than at a
#'   boundary.
#' @param predict_outliers Logical; whether `posterior_predict()` and
#'   `posterior_epred()` should include the outlier component. `FALSE` (the
#'   default in `rt_loggamma()`) fixes `poutlier` to zero for prediction, so
#'   predictions describe the decision process alone; the likelihood is always
#'   the full mixture either way. On the prediction methods themselves the
#'   default is `NULL`, which defers to the flag carried on the model - see
#'   [with_outliers()] to change it after fitting.
#' @export
rt_loggamma <- function(
  link_mu = "identity",
  link_sigma = "softplus",
  link_shape = "identity",
  link_ndt = "log",
  link_poutlier = "logit",
  predict_outliers = FALSE,
  minrt = 0.3
) {
  fam <- brms::custom_family(
    name = "rt_loggamma",
    dpars = c("mu", "sigma", "shape", "ndt", "poutlier"),
    links = c(link_mu, link_sigma, link_shape, link_ndt, link_poutlier),
    lb = c(NA, 0, NA, 0, 0), # sigma > 0, ndt > 0, poutlier >= 0
    ub = c(NA, NA, NA, NA, 1), # poutlier <= 1
    type = "real" # Continuous outcome variable (RT)
  )
  # Both ride on the family for the same reasons as in rt_lognormal(): it is the
  # only thing brms carries down to a custom family's prediction methods, and a
  # dpar left out of the formula would be estimated rather than defaulted.
  fam$predict_outliers <- isTRUE(predict_outliers)
  invisible(.validate_minrt(minrt)) # validate before storing
  fam$minrt <- minrt
  fam
}


#' @keywords internal
.rt_loggamma_lpdf <- function(minrt = 0.3) {
  # The outlier scale is baked in as a literal rather than passed as a dpar:
  # Stan functions cannot see the data block, and a dpar would be estimated
  # whenever the user left it out of the formula.
  # width = 1 so formatC does not pad the literal out with spaces
  scale <- formatC(
    .validate_minrt(minrt),
    format = "g",
    digits = 15,
    width = 1
  )
  sprintf(
    "
// Normalizing constant of the standardized log-gamma with shape k = 1 / shape^2,
// with its leading k already cancelled against the equal and opposite k hidden
// in the kernel below. Doing that cancellation on paper is what keeps the
// density accurate near shape = 0, where both terms reach 1e8 by |shape| = 1e-4.
// Below |shape| = 1e-3 the Stirling series is used, where 1 / (12k) = shape^2 / 12; the
// two branches agree to about 1e-9 where they meet, so the density is smooth
// through the switch and finite at shape = 0 exactly.
real rt_loggamma_lconst(real shape) {
  if (abs(shape) < 1e-3) {
    return -0.5 * log(2 * pi()) - square(shape) / 12 + pow(shape, 6) / 360;
  }
  real k = inv_square(shape);
  return (k - 0.5) * log(k) - lgamma(k) - k;
}

// (shape * w - expm1(shape * w)) / shape^2, the log-gamma kernel with the same k removed.
// At shape = 0 this is exactly -w^2 / 2, the Normal kernel, which is why shape = 0
// gives back the LogNormal rather than approaching it.
real rt_loggamma_lkernel(real shape, real w) {
  real u = shape * w;
  if (abs(u) < 1e-6) {
    return -0.5 * square(w) - shape * pow(w, 3) / 6 - square(shape) * pow(w, 4) / 24;
  }
  return (u - expm1(u)) / square(shape);
}

// Log-likelihood for one observation from the shifted Log-Gamma distribution.
// Y: observed reaction time.
// mu: location of the decision time on the log scale.
// sigma: scale of the decision time on the log scale (> 0).
// shape: shape of the log-gamma on the log-RT scale (unconstrained). shape = 0 is the
//    LogNormal, shape = sigma the Gamma, shape = 1 the Weibull.
// ndt: non-decision time, same unit as Y (> 0).
// poutlier: proportion of responses from the outlier process, in [0, 1].
//
// The outlier component is a half Student-t with 3 degrees of freedom and scale
// %s (= minrt). Its role is to keep the density strictly positive below `ndt`,
// where the shifted decision component has none. That is what removes the hard
// min-RT boundary and lets `ndt` be estimated directly instead of as a fraction
// of an observed minimum. Because the scale is tied to `minrt`, the likelihood
// is equivariant to the unit Y is measured in.
//
// Note that when sigma * shape >= 1 the decision density is unbounded at `ndt`, so
// the likelihood is unbounded too; the prior on shape is what keeps the sampler out
// of that region. See ?rrt_loggamma.
real rt_loggamma_lpdf(real Y, real mu, real sigma, real shape, real ndt, real poutlier) {
    // Parameter checks
    if (sigma <= 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (Y <= 0) return negative_infinity();

    // log(2) folds the symmetric Student-t onto [0, Inf).
    real lp_out = log(2) + student_t_lpdf(Y | 3, 0, %s);
    real t_adj  = Y - ndt;

    // Faster than the non-decision time: only the outlier component can have
    // produced this response.
    if (t_adj <= 0) return log(poutlier) + lp_out;

    real w = (log(t_adj) - mu) / sigma;
    real lp_dec = -log(sigma) - log(t_adj) + rt_loggamma_lconst(shape)
                  + rt_loggamma_lkernel(shape, w);

    return log_mix(poutlier, lp_out, lp_dec);
    }
",
    scale,
    scale
  )
}


#' @rdname rrt_loggamma
#' @export
rt_loggamma_lpdf_expose <- function(minrt = 0.3) {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
    "functions {
",
    .rt_loggamma_lpdf(.as_minrt(minrt)),
    "}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_loggamma_lpdf
}


#' @rdname rrt_loggamma
#' @export
rt_loggamma_stanvars <- function(minrt = 0.3) {
  brms::stanvar(
    scode = .rt_loggamma_lpdf(.as_minrt(minrt)),
    block = "functions"
  )
}



# brms methods ------------------------------------------------------------

#' @rdname rrt_loggamma
#' @inheritParams lnr
#' @export
log_lik_rt_loggamma <- function(i, prep) {
  # Extract observation
  if (!"Y" %in% names(prep$data)) {
    stop("Outcome variable 'Y' not found in prep$data.")
  }
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- brms::get_dpar(prep, "poutlier", i = i)

  n_draws <- max(
    length(mu), length(sigma), length(shape), length(ndt), length(poutlier)
  )
  if (n_draws == 0) return(numeric(0))

  # Compute log-likelihood using the vectorized drt_loggamma function
  ll <- drt_loggamma(
    x = rep(y, length.out = n_draws),
    mu = mu,
    sigma = sigma,
    shape = shape,
    ndt = ndt,
    poutlier = poutlier,
    minrt = .minrt(prep),
    log = TRUE
  )

  ll[is.na(ll)] <- -Inf
  ll
}


#' @rdname rrt_loggamma
#' @export
posterior_predict_rt_loggamma <- function(
  i,
  prep,
  predict_outliers = NULL,
  ...
) {
  # Get parameters for observation i across all draws
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  shape <- brms::get_dpar(prep, "shape", i = i)
  ndt <- brms::get_dpar(prep, "ndt", i = i)
  poutlier <- if (.predict_outliers(predict_outliers, prep)) {
    brms::get_dpar(prep, "poutlier", i = i)
  } else {
    0
  }

  n_draws <- max(length(mu), length(sigma), length(shape), length(ndt))

  # Simulate using rrt_loggamma (vectorized over its parameters)
  out <- rrt_loggamma(
    n = n_draws,
    mu = mu,
    sigma = sigma,
    shape = shape,
    ndt = ndt,
    poutlier = poutlier,
    minrt = .minrt(prep)
  )

  # Return as a matrix (draws x 1)
  as.matrix(out)
}


#' @rdname rrt_loggamma
#' @export
posterior_epred_rt_loggamma <- function(prep, predict_outliers = NULL) {
  # Extract draws (matrices: draws x observations)
  mu <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  shape <- brms::get_dpar(prep, "shape")
  ndt <- brms::get_dpar(prep, "ndt")

  # Expected RT of the decision process alone. Inf where sigma * abs(shape) >= 1
  # with shape < 0, i.e. where the power-law right tail has no mean.
  #
  # .mloggamma() recycles to a flat vector, so the draws x observations shape
  # has to be put back. get_dpar() returns a scalar for a dpar that is constant
  # across draws and observations, so the shape is taken from whichever one is
  # a matrix rather than from `mu`.
  dims <- NULL
  for (d in list(mu, sigma, shape, ndt)) {
    if (is.matrix(d)) {
      dims <- dim(d)
      break
    }
  }
  mean_dec <- .mloggamma(mu = mu, sigma = sigma, shape = shape)
  if (!is.null(dims) && length(mean_dec) == prod(dims)) {
    dim(mean_dec) <- dims
  }
  mean_dec <- mean_dec + ndt

  if (!.predict_outliers(predict_outliers, prep)) {
    return(mean_dec)
  }

  # E[RT] of the mixture
  poutlier <- brms::get_dpar(prep, "poutlier")
  (1 - poutlier) * mean_dec + poutlier * .mcontam(.minrt(prep))
}
