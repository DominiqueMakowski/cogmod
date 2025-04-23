#' @keywords internal
.rt_invgaussian_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Wald distribution.
// Calculation is done directly without a separate inv_gaussian helper.
// Y: observed reaction time.
// mu: drift rate (equivalent to drift in R functions). Must be positive.
// bs: decision threshold. Must be positive.
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (used to scale tau).
real rt_invgaussian_lpdf(real Y, real mu, real bs, real tau, real minrt) {
  // Parameter checks
  if (mu <= 0 || bs <= 0 || tau < 0 || tau > 1 || minrt < 0) return negative_infinity();

  // Compute non-decision time and adjusted time
  real ndt   = tau * minrt;
  real t_adj = Y - ndt;
  if (t_adj <= 0) return negative_infinity();

  // Inverse-Gaussian parameters
  real ig_mu     = bs / mu;
  real ig_lambda = square(bs);

  // Log-density
  real log_term = 0.5 * (log(ig_lambda) - (log(2) + log(pi())) - 3 * log(t_adj));
  real diff     = t_adj - ig_mu;
  real exponent = ig_lambda * square(diff) / (2 * square(ig_mu) * t_adj);

  return log_term - exponent;
}
"
}


#' @rdname rrt_invgaussian
#' @export
rt_invgaussian_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {
", .rt_invgaussian_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_invgaussian_lpdf
}


#' @rdname rrt_invgaussian
#' @export
rt_invgaussian_stanvars <- function() {
  brms::stanvar(scode = .rt_invgaussian_lpdf(), block = "functions")
}


#' @rdname rrt_invgaussian
#' @param link_mu,link_bs,link_tau,link_minrt Link functions for the parameters.
#' @export
rt_invgaussian <- function(link_mu = "softplus", link_bs = "softplus", link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_invgaussian",
    dpars = c("mu", "bs", "tau", "minrt"),
    links = c(link_mu, link_bs, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: mu>0, bs>0, tau>=0, minrt>=0
    ub = c(NA, NA, 1, NA), # Upper bound: tau<=1
    type = "real" # Continuous outcome variable (RT)
  )
}


# brms methods ------------------------------------------------------------

#' @rdname rrt_invgaussian
#' @inheritParams rbetagate
#' @export
log_lik_rt_invgaussian <- function(i, prep, ...) {
  # Extract observation
  # Assuming the response variable is named 'Y' in the data
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_) # Return NA if outcome is NA

  # Get parameters for observation i across all draws
  mu    <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Get minrt parameter

  # Determine number of draws
  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0)) # Handle case with no draws

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y, length.out = n_draws)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Calculate log-likelihood using the vectorized drt_invgaussian function
  # Note: drt_invgaussian uses 'drift', 'bs', 'ndt'
  ll <- drt_invgaussian(x = y_vec, drift = mu, bs = bs, ndt = ndt, log = TRUE, ...)

  # Ensure no NaN/NA values (drt_invgaussian should return -Inf for zero density)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll # Return the vector of log-likelihoods for all draws
}


#' @rdname rrt_invgaussian
#' @inheritParams rbetagate
#' @export
posterior_predict_rt_invgaussian <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  mu    <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Number of posterior draws
  n_draws <- length(mu)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Simulate using rrt_invgaussian (vectorized)
  # Note: rrt_invgaussian uses 'drift', 'bs', 'ndt'
  final_out <- rrt_invgaussian(n = n_draws, drift = mu, bs = bs, ndt = ndt)

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}


#' @rdname rrt_invgaussian
#' @inheritParams rbetagate
#' @export
posterior_epred_rt_invgaussian <- function(prep) {
  # Extract draws for the necessary parameters (matrices: draws x observations)
  mu    <- brms::get_dpar(prep, "mu")
  bs <- brms::get_dpar(prep, "bs")
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedWald] = E[InverseGaussian] + ndt
  # E[InverseGaussian(mean=bs/mu, shape=bs^2)] = bs / mu
  epred <- (bs / mu) + ndt

  epred # Return the matrix of posterior expectations (draws x observations)
}