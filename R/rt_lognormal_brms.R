#' Shifted Lognormal Model
#'
#' Provides the necessary functions to use a shifted lognormal distribution
#' as a custom family in `brms`. This version is reparameterized to use
#' `tau` (proportion of non-decision time relative to minimum RT) and `minrt`
#' (minimum possible RT), where the non-decision time `ndt = tau * minrt`.
#' The distribution is parameterized by the mean (`meanlog`, named `mu` in brms)
#' and standard deviation (`sigma`) of the distribution on the log scale (`sdlog`),
#' as well as `tau`, and `minrt` (which is typically fixed).
#'
#' @param link_mu,link_sigma,link_tau,link_minrt Link functions for the parameters.
#'
#' @export
rt_lognormal <- function(link_mu = "identity", link_sigma = "softplus",
                         link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_lognormal",
    dpars = c("mu", "sigma", "tau", "minrt"),
    links = c(link_mu, link_sigma, link_tau, link_minrt), # Updated links
    lb = c(NA, 0, 0, 0), # Lower bounds: sigma > 0, tau >= 0, minrt > 0
    ub = c(NA, NA, 1, NA), # Upper bound: tau <= 1
    type = "real" # Continuous outcome variable (RT)
  )
}





#' @keywords internal
.rt_lognormal_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Lognormal distribution.
// Y: observed reaction time.
// mu: mean of the distribution on the log scale (meanlog).
// sigma: standard deviation of the distribution on the log scale (> 0).
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (> 0).
real rt_lognormal_lpdf(real Y, real mu, real sigma, real tau, real minrt) {
    // Parameter checks
    if (sigma <= 0 || tau < 0 || tau > 1 || minrt <= 0) return negative_infinity();

    // Compute non-decision time and adjusted time
    real ndt   = tau * minrt;
    real t_adj = Y - ndt;
    if (t_adj <= 0) return negative_infinity(); // Density is 0 if Y <= ndt

    // Use Stan's built-in lognormal lpdf
    return lognormal_lpdf(t_adj | mu, sigma);
    }
"
}

#' @rdname rt_lognormal
#' @export
rt_lognormal_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {
", .rt_lognormal_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_lognormal_lpdf # Corrected function name
}


#' @rdname rt_lognormal
#' @export
rt_lognormal_stanvars <- function() {
  brms::stanvar(scode = .rt_lognormal_lpdf(), block = "functions")
}



# brms methods ------------------------------------------------------------

#' @rdname rt_lognormal
#' @inheritParams lnr
#' @export
log_lik_rt_lognormal <- function(i, prep) { # Renamed function
  # Extract observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  # Get parameters for observation i across all draws
  meanlog <- brms::get_dpar(prep, "mu", i = i)
  sdlog   <- brms::get_dpar(prep, "sigma", i = i)
  tau     <- brms::get_dpar(prep, "tau", i = i) # Get tau
  minrt   <- brms::get_dpar(prep, "minrt", i = i) # Get minrt

  # Determine number of draws
  n_draws <- length(meanlog)
  if (n_draws == 0) return(numeric(0))

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y, length.out = n_draws)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt # Calculate ndt from tau and minrt

  # Calculate log-likelihood using R's dlnorm (vectorized)
  ll <- stats::dlnorm(x = y_vec - ndt, meanlog = meanlog, sdlog = sdlog, log = TRUE)

  # Ensure correct handling for y <= ndt
  ll[y_vec <= ndt] <- -Inf

  # Ensure no other NaN/NA values
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll
}


#' @rdname rt_lognormal
#' @export
posterior_predict_rt_lognormal <- function(i, prep, ...) { # Renamed function
  # Get parameters for observation i across all draws
  meanlog <- brms::get_dpar(prep, "mu", i = i)
  sdlog   <- brms::get_dpar(prep, "sigma", i = i)
  tau     <- brms::get_dpar(prep, "tau", i = i) # Get tau
  minrt   <- brms::get_dpar(prep, "minrt", i = i) # Get minrt

  # Number of posterior draws
  n_draws <- length(meanlog)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt # Calculate ndt

  # Simulate using rlnorm (vectorized) and add the shift
  lognormal_part <- stats::rlnorm(n = n_draws, meanlog = meanlog, sdlog = sdlog)
  final_out <- lognormal_part + ndt # Add calculated ndt

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}


#' @rdname rt_lognormal
#' @export
posterior_epred_rt_lognormal <- function(prep) { # Renamed function
  # Extract draws for the necessary parameters (matrices: draws x observations)
  meanlog <- brms::get_dpar(prep, "mu")
  sdlog   <- brms::get_dpar(prep, "sigma")
  tau     <- brms::get_dpar(prep, "tau") # Get tau
  minrt   <- brms::get_dpar(prep, "minrt") # Get minrt

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt # Calculate ndt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedLognormal] = E[Lognormal] + ndt
  epred <- exp(meanlog + sdlog^2 / 2) + ndt # Add calculated ndt

  epred
}
