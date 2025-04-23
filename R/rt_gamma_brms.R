#' Shifted Gamma Model for brms
#'
#' Provides the necessary functions to use a shifted Gamma distribution
#' as a custom family in `brms`. This version is parameterized using
#' `tau` (proportion of non-decision time relative to minimum RT) and `minrt`
#' (minimum possible RT), where the non-decision time `ndt = tau * minrt`.
#' The distribution is parameterized by the shape (`mu` in brms, alpha in standard Gamma)
#' and scale (`sigma` in brms, beta in standard Gamma),
#' as well as `tau`, and `minrt` (which is typically fixed).
#'
#' @param link_mu,link_sigma,link_tau,link_minrt Link functions for the parameters.
#'
#' @return A `brms::custom_family` object.
#' @export
rt_gamma <- function(link_mu = "softplus", link_sigma = "softplus",
                     link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_gamma",
    dpars = c("mu", "sigma", "tau", "minrt"), # mu = shape, sigma = scale
    links = c(link_mu, link_sigma, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: shape > 0, scale > 0, tau >= 0, minrt > 0
    ub = c(NA, NA, 1, NA), # Upper bound: tau <= 1
    type = "real" # Continuous outcome variable (RT)
  )
}


#' @keywords internal
.rt_gamma_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Gamma distribution.
// Y: observed reaction time.
// mu: shape parameter (alpha > 0).
// sigma: scale parameter (beta > 0).
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (> 0).
real rt_gamma_lpdf(real Y, real mu, real sigma, real tau, real minrt) {
    // Parameter checks
    if (mu <= 0 || sigma <= 0 || tau < 0 || tau > 1 || minrt <= 0) return negative_infinity();

    // Compute non-decision time and adjusted time
    real ndt   = tau * minrt;
    real t_adj = Y - ndt;
    if (t_adj <= 0) return negative_infinity(); // Density is 0 if Y <= ndt

    // Use Stan's built-in gamma lpdf (shape, rate = 1/scale)
    return gamma_lpdf(t_adj | mu, 1 / sigma);
}
"
}

#' @rdname rt_gamma
#' @export
rt_gamma_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {
", .rt_gamma_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_gamma_lpdf
}


#' @rdname rt_gamma
#' @export
rt_gamma_stanvars <- function() {
  brms::stanvar(scode = .rt_gamma_lpdf(), block = "functions")
}



# brms methods ------------------------------------------------------------

#' @rdname rt_gamma
#' @export
log_lik_rt_gamma <- function(i, prep) {
  # Extract observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)    # mu = shape
  scale <- brms::get_dpar(prep, "sigma", i = i) # sigma = scale
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

  # Calculate log-likelihood using R's dgamma (vectorized)
  # Note: dgamma uses shape and scale OR shape and rate. We use shape and scale.
  ll <- stats::dgamma(x = t_adj, shape = shape, scale = scale, log = TRUE)

  # Ensure correct handling for y <= ndt or invalid parameters
  ll[t_adj <= 0 | shape <= 0 | scale <= 0 | tau < 0 | tau > 1 | minrt <= 0] <- -Inf

  # Ensure no other NaN/NA values (e.g., from invalid parameter combinations)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll
}


#' @rdname rt_gamma
#' @inheritParams lnr
#' @export
posterior_predict_rt_gamma <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  shape <- brms::get_dpar(prep, "mu", i = i)    # mu = shape
  scale <- brms::get_dpar(prep, "sigma", i = i) # sigma = scale
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Number of posterior draws
  n_draws <- length(shape)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Simulate using rgamma (vectorized) and add the shift
  # Note: rgamma uses shape and scale OR shape and rate. We use shape and scale.
  # Handle cases where parameters might be invalid for rgamma (e.g., shape/scale <= 0)
  # We generate NA for invalid parameters and rely on subsequent steps if needed,
  # though ideally priors/links prevent this.
  gamma_part <- stats::rgamma(n = n_draws, shape = pmax(shape, .Machine$double.eps), scale = pmax(scale, .Machine$double.eps))
  final_out <- gamma_part + ndt

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}


#' @rdname rt_gamma
#' @export
posterior_epred_rt_gamma <- function(prep) {
  # Extract draws for the necessary parameters (matrices: draws x observations)
  shape <- brms::get_dpar(prep, "mu")    # mu = shape
  scale <- brms::get_dpar(prep, "sigma") # sigma = scale
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedGamma] = E[Gamma] + ndt = shape * scale + ndt
  epred <- shape * scale + ndt

  # Handle potential invalid parameters leading to NA/NaN, though less likely for epred
  epred[shape <= 0 | scale <= 0 | tau < 0 | tau > 1 | minrt <= 0] <- NA_real_

  epred
}
