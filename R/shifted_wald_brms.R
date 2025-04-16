#' @keywords internal
.shifted_wald_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Wald distribution.
// Calculation is done directly without a separate inv_gaussian helper.
// Y: observed reaction time.
// mu: drift rate (equivalent to nu in R functions). Must be positive.
// alpha: decision threshold. Must be positive.
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (used to scale tau).
real shifted_wald_lpdf(real Y, real mu, real alpha, real tau, real minrt) {
  real eps = 1e-10;
  
  if (mu <= 0.0 || alpha <= 0.0 || tau < 0.0 || tau > 1.0 || minrt < 0.0 ||
      is_nan(mu) || is_inf(mu) ||
      is_nan(alpha) || is_inf(alpha) ||
      is_nan(tau) || is_inf(tau) ||
      is_nan(minrt) || is_inf(minrt))
    return negative_infinity();

  real ndt = tau * minrt;
  real t_adj = Y - ndt;
  if (t_adj <= eps) return negative_infinity();

  real ig_mean = alpha / mu;
  real ig_lambda = square(alpha);
  // Precompute constant to save time on every evaluation.
  real log_two_pi = log(2 * pi());
  return 0.5 * (log(ig_lambda) - log_two_pi - 3 * log(t_adj)) - (ig_lambda * square(t_adj - ig_mean)) / (2 * square(ig_mean) * t_adj);
}
"
}


#' @rdname rshifted_wald
#' @examples
#' \dontrun{
#' # You can expose the lpdf function as follows:
#' if (requireNamespace("cmdstanr", quietly = TRUE)) {
#'   sw_lpdf <- shifted_wald_lpdf_expose()
#'   # Example call (using nu, alpha, ndt for clarity, map to mu, alpha, tau*minrt)
#'   sw_lpdf(Y = 0.5, mu = 3, alpha = 0.5, tau = 0.5, minrt = 0.4) # ndt = 0.2
#' }
#' }
#' @export
shifted_wald_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {
", .shifted_wald_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$shifted_wald_lpdf
}


#' @rdname rshifted_wald
#' @export
shifted_wald_stanvars <- function() {
  brms::stanvar(scode = .shifted_wald_lpdf(), block = "functions")
}


#' @rdname rshifted_wald
#' @param link_mu,link_alpha,link_tau,link_minrt Link functions.
#' @export
#' @examples
#' \dontrun{
#' # Example brms formula using shifted_wald family
#' # bf(rt ~ condition + (1|subject),
#' #    alpha ~ 1,
#' #    tau ~ condition,
#' #    minrt ~ 1,  # Often fixed or estimated per subject
#' #    family = shifted_wald())
#' }
shifted_wald <- function(link_mu = "log", link_alpha = "log",
                         link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "shifted_wald",
    dpars = c("mu", "alpha", "tau", "minrt"),
    links = c(link_mu, link_alpha, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: mu>0, alpha>0, tau>=0, minrt>=0
    ub = c(NA, NA, 1, NA), # Upper bound: tau<=1
    type = "real" # Continuous outcome variable (RT)
  )
}


# brms methods ------------------------------------------------------------

#' @rdname rshifted_wald
#' @inheritParams rbext
#' @export
log_lik_shifted_wald <- function(i, prep) {
  # Extract observation
  # Assuming the response variable is named 'Y' in the data
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_) # Return NA if outcome is NA

  # Get parameters for observation i across all draws
  mu    <- brms::get_dpar(prep, "mu", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Get minrt parameter

  # Determine number of draws
  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0)) # Handle case with no draws

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y, length.out = n_draws)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Calculate log-likelihood using the vectorized dshifted_wald function
  # Note: dshifted_wald uses 'nu', 'alpha', 'ndt'
  ll <- dshifted_wald(x = y_vec, nu = mu, alpha = alpha, ndt = ndt, log = TRUE)

  # Ensure no NaN/NA values (dshifted_wald should return -Inf for zero density)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll # Return the vector of log-likelihoods for all draws
}


#' @rdname rshifted_wald
#' @inheritParams rbext
#' @export
posterior_predict_shifted_wald <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  mu    <- brms::get_dpar(prep, "mu", i = i)
  alpha <- brms::get_dpar(prep, "alpha", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Number of posterior draws
  n_draws <- length(mu)

  # Calculate non-decision time (vectorized)
  ndt <- tau * minrt

  # Simulate using rshifted_wald (vectorized)
  # Note: rshifted_wald uses 'nu', 'alpha', 'ndt'
  final_out <- rshifted_wald(n = n_draws, nu = mu, alpha = alpha, ndt = ndt)

  # Return as a matrix (draws x 1)
  as.matrix(final_out)
}


#' @rdname rshifted_wald
#' @inheritParams rbext
#' @export
posterior_epred_shifted_wald <- function(prep) {
  # Extract draws for the necessary parameters (matrices: draws x observations)
  mu    <- brms::get_dpar(prep, "mu")
  alpha <- brms::get_dpar(prep, "alpha")
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  # Calculate non-decision time (matrix: draws x observations)
  ndt <- tau * minrt

  # Calculate the expectation (mean) for each draw and observation
  # E[ShiftedWald] = E[InverseGaussian] + ndt
  # E[InverseGaussian(mean=alpha/mu, shape=alpha^2)] = alpha / mu
  epred <- (alpha / mu) + ndt

  # Check for potential division by zero or invalid results (should be handled by priors/links)
  epred[!is.finite(epred)] <- NA # Or handle more gracefully if needed

  epred # Return the matrix of posterior expectations (draws x observations)
}