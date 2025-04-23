#' Shifted Inverse Gamma Model for brms
#'
#' Provides the necessary functions to use a shifted Inverse Gamma distribution
#' as a custom family in `brms`. This version is parameterized using
#' `tau` (proportion of non-decision time relative to minimum RT) and `minrt`
#' (minimum possible RT), where the non-decision time `ndt = tau * minrt`.
#' The distribution is parameterized by the shape (`mu` in brms, alpha in standard InvGamma)
#' and scale (`sigma` in brms, beta in standard InvGamma),
#' as well as `tau`, and `minrt` (which is typically fixed).
#'
#' @param link_mu,link_sigma,link_tau,link_minrt Link functions for the parameters.
#'
#' @return A `brms::custom_family` object.
#' @export
rt_invgamma <- function(link_mu = "softplus", link_sigma = "softplus",
                        link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_invgamma",
    dpars = c("mu", "sigma", "tau", "minrt"), # mu = shape (alpha), sigma = scale (beta)
    links = c(link_mu, link_sigma, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: shape > 0, scale > 0, tau >= 0, minrt > 0
    ub = c(NA, NA, 1, NA), # Upper bound: tau <= 1
    type = "real" # Continuous outcome variable (RT)
  )
}


#' @keywords internal
.rt_invgamma_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Inverse Gamma distribution.
// Y: observed reaction time.
// mu: shape parameter (alpha > 0).
// sigma: scale parameter (beta > 0).
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (> 0).
real rt_invgamma_lpdf(real Y, real mu, real sigma, real tau, real minrt) {
    // Parameter checks
    if (mu <= 0 || sigma <= 0 || tau < 0 || tau > 1 || minrt <= 0) return negative_infinity();

    // Compute non-decision time and adjusted time
    real ndt   = tau * minrt;
    real t_adj = Y - ndt;
    if (t_adj <= 0) return negative_infinity(); // Density is 0 if Y <= ndt

    // Use Stan's built-in inv_gamma lpdf (shape, scale)
    return inv_gamma_lpdf(t_adj | mu, sigma);
}
"
}

#' Expose Stan function for Shifted Inverse Gamma LPDF
#'
#' Makes the Stan implementation of the Shifted Inverse Gamma log probability density
#' function available in R for testing or direct use. Requires `cmdstanr`.
#'
#' @return The exposed Stan function `rt_invgamma_lpdf`.
#' @export
rt_invgamma_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {
", .rt_invgamma_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_invgamma_lpdf
}


#' Stan variables for Shifted Inverse Gamma
#'
#' Provides the Stan code block necessary for the Shifted Inverse Gamma distribution
#' to be used within `brms`.
#'
#' @return A `brms::stanvar` object containing the Stan function code.
#' @export
rt_invgamma_stanvars <- function() {
  brms::stanvar(scode = .rt_invgamma_lpdf(), block = "functions")
}



# brms methods ------------------------------------------------------------

#' Log-likelihood for Shifted Inverse Gamma
#'
#' Calculates the log-likelihood for the Shifted Inverse Gamma distribution for
#' use with `brms`. Uses base R functions.
#'
#' @param i Observation index.
#' @param prep A `brms::brmsprep` object.
#' @return A vector of log-likelihood values for observation `i` across posterior draws.
#' @export
log_lik_rt_invgamma <- function(i, prep) {
  # Extract observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)    # mu = shape (alpha)
  scale <- brms::get_dpar(prep, "sigma", i = i) # sigma = scale (beta)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Determine number of draws
  n_draws <- length(shape)
  if (n_draws == 0) return(numeric(0))

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y, length.out = n_draws)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Calculate adjusted time
  t_adj <- y_vec - ndt

  # Calculate log-likelihood using base R functions (vectorized)
  # log(f(x | alpha, beta)) = alpha * log(beta) - lgamma(alpha) - (alpha + 1) * log(x) - beta / x
  # Use pmax to avoid issues with log(0) or log(<0) if t_adj is non-positive
  # (although these cases should be handled by the subsequent check)
  log_t_adj <- log(pmax(t_adj, .Machine$double.eps))
  ll <- shape * log(scale) - lgamma(shape) - (shape + 1) * log_t_adj - scale / t_adj

  # Ensure correct handling for y <= ndt or invalid parameters
  ll[t_adj <= 0 | shape <= 0 | scale <= 0 | tau < 0 | tau > 1 | minrt <= 0] <- -Inf

  # Ensure no other NaN/NA values (e.g., from invalid parameter combinations)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll
}


#' Posterior predictions for Shifted Inverse Gamma
#'
#' Generates posterior predictions for the Shifted Inverse Gamma distribution for
#' use with `brms`. Uses base R functions.
#'
#' @param i Observation index.
#' @param prep A `brms::brmsprep` object.
#' @param ... Additional arguments (unused).
#' @return A matrix of posterior predictions (draws x 1).
#' @export
posterior_predict_rt_invgamma <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)    # mu = shape (alpha)
  scale <- brms::get_dpar(prep, "sigma", i = i) # sigma = scale (beta)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Number of posterior draws
  n_draws <- length(shape)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Simulate using 1 / rgamma (vectorized) and add the shift
  # InvGamma(shape, scale) is equivalent to 1 / Gamma(shape, rate = scale)
  # Handle cases where parameters might be invalid for rgamma (e.g., shape/scale <= 0)
  gamma_shape <- pmax(shape, .Machine$double.eps)
  gamma_rate <- pmax(scale, .Machine$double.eps)
  invgamma_part <- 1 / stats::rgamma(n = n_draws, shape = gamma_shape, rate = gamma_rate)
  final_out <- invgamma_part + ndt

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}


#' Posterior expected value for Shifted Inverse Gamma
#'
#' Calculates the posterior expected value (mean) for the Shifted Inverse Gamma
#' distribution for use with `brms`.
#'
#' @param prep A `brms::brmsprep` object.
#' @return A matrix of posterior expected values (draws x observations).
#' @export
posterior_epred_rt_invgamma <- function(prep) {
  # Extract draws for the necessary parameters (matrices: draws x observations)
  shape <- brms::get_dpar(prep, "mu")    # mu = shape (alpha)
  scale <- brms::get_dpar(prep, "sigma") # sigma = scale (beta)
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedInvGamma] = E[InvGamma] + ndt = scale / (shape - 1) + ndt (for shape > 1)
  # Set to Inf if shape <= 1, as the mean is undefined or infinite.
  epred <- ifelse(shape > 1, scale / (shape - 1), Inf) + ndt

  # Handle potential invalid parameters leading to NA/NaN
  epred[shape <= 0 | scale <= 0 | tau < 0 | tau > 1 | minrt <= 0] <- NA_real_

  epred
}
