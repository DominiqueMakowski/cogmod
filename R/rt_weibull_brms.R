#' Shifted Weibull Model
#'
#' Provides the necessary functions to use a shifted Weibull distribution
#' as a custom family in `brms`. This version is reparameterized to use
#' `tau` (proportion of non-decision time relative to minimum RT) and `minrt`
#' (minimum possible RT), where the non-decision time `ndt = tau * minrt`.
#' The distribution is parameterized by the shape (`mu`), scale (`sigma`),
#' as well as `tau`, and `minrt` (which is typically fixed).
#'
#' @param link_mu,link_sigma,link_tau,link_minrt Link functions for the parameters.
#'
#' @export
rt_weibull <- function(link_mu = "softplus", link_sigma = "softplus",
                       link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_weibull",
    dpars = c("mu", "sigma", "tau", "minrt"),
    links = c(link_mu, link_sigma, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: mu > 0, sigma > 0, tau >= 0, minrt > 0
    ub = c(NA, NA, 1, NA), # Upper bound: tau <= 1
    type = "real" # Continuous outcome variable (RT)
  )
}

#' @keywords internal
.rt_weibull_lpdf <- function() {
"// Log-likelihood for a single observation from the Shifted Weibull distribution.
// Y: observed reaction time.
// mu: shape parameter of the Weibull distribution (> 0).
// sigma: scale parameter of the Weibull distribution (> 0).
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (> 0).
real rt_weibull_lpdf(real Y, real mu, real sigma, real tau, real minrt) {
    // Parameter checks
    if (mu <= 0 || sigma <= 0 || tau < 0 || tau > 1 || minrt <= 0) return negative_infinity();

    // Compute non-decision time and adjusted time
    real ndt   = tau * minrt;
    real t_adj = Y - ndt;
    if (t_adj <= 0) return negative_infinity(); // Density is 0 if Y <= ndt

    // Use Stan's built-in Weibull lpdf
    return weibull_lpdf(t_adj | mu, sigma);
    }
"
}

#' @rdname rt_weibull
#' @export
rt_weibull_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {",
.rt_weibull_lpdf(),
"}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_weibull_lpdf
}

#' @rdname rt_weibull
#' @export
rt_weibull_stanvars <- function() {
  brms::stanvar(scode = .rt_weibull_lpdf(), block = "functions")
}

# brms methods ------------------------------------------------------------

#' @rdname rt_weibull
#' @inheritParams lnr
#' @export
log_lik_rt_weibull <- function(i, prep) {
  # Extract observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)
  scale <- brms::get_dpar(prep, "sigma", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Determine number of draws
  n_draws <- length(shape)
  if (n_draws == 0) return(numeric(0))

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y, length.out = n_draws)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Calculate log-likelihood using R's dweibull (vectorized)
  ll <- stats::dweibull(x = y_vec - ndt, shape = shape, scale = scale, log = TRUE)

  # Ensure correct handling for y <= ndt
  ll[y_vec <= ndt] <- -Inf

  # Ensure no other NaN/NA values
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll
}

#' @rdname rt_weibull
#' @export
posterior_predict_rt_weibull <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)
  scale <- brms::get_dpar(prep, "sigma", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Number of posterior draws
  n_draws <- length(shape)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Simulate using rweibull (vectorized) and add the shift
  weibull_part <- stats::rweibull(n = n_draws, shape = shape, scale = scale)
  final_out <- weibull_part + ndt

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}

#' @rdname rt_weibull
#' @export
posterior_epred_rt_weibull <- function(prep) {
  # Extract draws for the necessary parameters (matrices: draws x observations)
  shape <- brms::get_dpar(prep, "mu")
  scale <- brms::get_dpar(prep, "sigma")
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedWeibull] = scale * Gamma(1 + 1/shape) + ndt
  epred <- scale * gamma(1 + 1 / shape) + ndt

  epred
}
