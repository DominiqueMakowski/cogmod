# Stanvars ----------------------------------------------------------------

#' @keywords internal
.lnr_lpdf <- function() {
"
// Log-likelihood for a single observation from the Log-Normal Race model.
// Uses the 'nu' parameterization where meanlog = -nu (higher nu = faster).
// Y: observed reaction time.
// dec: decision indicator (0 or 1).
// mu: nuzero (named mu as per Stan requirement). processing speed parameter for accumulator 0.
// nuone: processing speed parameter for accumulator 1.
// sigmazero: log-space standard deviation for accumulator 0.
// sigmaone: log-space standard deviation for accumulator 1.
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (used to scale tau).
real lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real tau, real minrt, int dec) {
  real eps = 1e-9; // Small epsilon for numerical stability

  // --- 1. Input Validation ---
  // Check parameters are finite and within valid ranges
  if (sigmazero <= 0.0 || sigmaone <= 0.0 || tau < 0.0 || tau > 1.0 || minrt < 0.0) return negative_infinity();
  
  // Check decision indicator is valid
  if (dec != 0 && dec != 1) return negative_infinity();

  // --- 2. Derived Parameters ---
  // Ensure sigmas are positive
  real sig0 = fmax(sigmazero, eps);
  real sig1 = fmax(sigmaone, eps);
  // Calculate non-decision time
  real ndt = tau * minrt;
  // Calculate adjusted time (decision time)
  real t_adj = Y - ndt;
  // Check if adjusted time is valid (must be positive)
  if (t_adj < eps) return negative_infinity();

  // Calculate meanlog parameters from nu parameters
  real meanlog0 = -mu;
  real meanlog1 = -nuone;
  real log_lik;

  // --- 3. Log-likelihood using built-in lognormal functions ---
  // Calculate log-likelihood based on the winning accumulator (dec)
  if (dec == 0) {
    // Accumulator 0 finished first
    log_lik = lognormal_lpdf(t_adj | meanlog0, sig0) + // Log-PDF of winner
              lognormal_lccdf(t_adj | meanlog1, sig1); // Log-Survival of loser
  } else {
    // Accumulator 1 finished first
    log_lik = lognormal_lpdf(t_adj | meanlog1, sig1) + // Log-PDF of winner
              lognormal_lccdf(t_adj | meanlog0, sig0); // Log-Survival of loser
  }

  // Return negative infinity if log-likelihood is NaN or Inf
  return (is_nan(log_lik) || is_inf(log_lik)) ? negative_infinity() : log_lik;
}
"
}

#' @rdname rlnr
#' @examples
#' \dontrun{
#' # You can expose the lpdf function as follows:
#' insight::check_if_installed("cmdstanr")
#' lnr_lpdf <- lnr_lpdf_expose()
#' # Example call with nu parameterization (higher nu = faster)
#' lnr_lpdf(Y = 0.5, mu = 0.5, nuone = 0.2, sigmazero = 1.0, sigmaone = 0.8,
#'          tau = 0.4, minrt = 0.2, dec = 0) # ndt = 0.4 * 0.2 = 0.08
#' }
#'
#' @export
lnr_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
"functions {
", .lnr_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$lnr_lpdf
}




#' @rdname rlnr
#' @export
lnr_stanvars <- function() {
  brms::stanvar(scode = .lnr_lpdf(), block = "functions")
}


#' @rdname rlnr
#' @param link_nuzero,link_nuone Link function for the nu parameters.
#' @param link_sigmazero,link_sigmaone Link function for the sigma parameters.
#' @param link_tau Link function for the tau parameter (non-decision time proportion).
#' @param link_minrt Link function for the minrt parameter (minimum RT scale).
#' @export
lnr <- function(link_nuzero = "identity", link_nuone = "identity",
                link_sigmazero = "softplus", link_sigmaone = "softplus",
                link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "lnr",
    dpars = c("mu", "nuone", "sigmazero", "sigmaone", "tau", "minrt"),
    links = c(link_nuzero, link_nuone, link_sigmazero, link_sigmaone, link_tau, link_minrt),
    lb = c(NA, NA, 0, 0, 0, 0), # Lower bounds for sigma > 0, tau >= 0, minrt >= 0
    ub = c(NA, NA, NA, NA, 1, NA), # Upper bound for tau <= 1
    type = "real", # Response variable type
    vars = "dec[n]" # Required additional data variable for the decision
  )
}



# # brms --------------------------------------------------------------------

#' @rdname rlnr
#' @importFrom brms get_dpar
log_lik_lnr <- function(i, prep) {
  # Extract observation (response variable, usually RT)
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_) # Return NA if response is missing

  # Get parameters for observation i across all posterior draws
  nuzero <- brms::get_dpar(prep, "mu", i = i)
  nuone <- brms::get_dpar(prep, "nuone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Get decision indicator for observation i
  response <- prep$data$dec[i] # Assumes 'dec' column exists in the data
  # Basic check for valid response coding (0 or 1)
  if (!response %in% c(0, 1)) {
    warning("Response ('dec') must be 0 or 1. Found: ", response, " for observation ", i, ". Returning -Inf log-likelihood.")
    # Return -Inf for all draws for this observation
    return(rep(-Inf, length(nuzero)))
  }

  # Compute log-likelihood using the R density function dlnr
  # dlnr is already vectorized over its parameters
  ll <- dlnr(x = y, nuzero = nuzero, nuone = nuone,
             sigmazero = sigmazero, sigmaone = sigmaone,
             ndt = ndt, response = response, log = TRUE)

  ll # Return vector of log-likelihoods (one per posterior draw)
}



#' @rdname rlnr
#' @inheritParams brms::posterior_predict.brmsfit
#' @importFrom brms get_dpar
#' @export
posterior_predict_lnr <- function(i, prep, ...) {
  # Get parameters for observation i across all posterior draws
  nuzero <- brms::get_dpar(prep, "mu", i = i)
  nuone <- brms::get_dpar(prep, "nuone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Generate predictions using the R random number generator rlnr
  # rlnr is already vectorized over its parameters
  n_draws <- length(nuzero) # Number of posterior draws
  sim_data <- rlnr(n = n_draws, nuzero = nuzero, nuone = nuone,
                   sigmazero = sigmazero, sigmaone = sigmaone,
                   ndt = ndt)

  # Return simulated data as a matrix (draws x variables)
  # The variables are 'rt' and 'response' from rlnr output
  as.matrix(sim_data)
}




# Note: It's often better not to define posterior_epred for complex models like race models
# where the expectation is computationally prohibitive to calculate accurately within
# this context. Users needing expected values are typically better off calculating them
# manually from the posterior draws using more optimized integration or simulation methods
# outside the main brms fitting/post-processing functions.
# However, if you must provide one, calculating P(response=1) is slightly less complex
# than E[RT]. Here's how it could be implemented, with a strong caveat about performance.

#' @rdname rlnr
#' @export
posterior_epred_lnr <- function(prep) {
  stop("Computing posterior_epred for complex models like race models is computationally prohibitive",
       " and the user is better off computing 'prediction' rather than 'expectation' and",
       " using the posterior draws to compute any quantities of interest.")
}


# --- Commented out internal helper for epred ---
# #' Calculate Marginal Probability for LNR (Internal Helper)
# #'
# #' Integrates the joint density dlnr to get the marginal probability of a response.
# #' NOTE: This involves numerical integration and can be slow.
# #' @keywords internal
# .posterior_epred_lnr_marginal_p <- function(nuzero, nuone, sigmazero, sigmaone, ndt, response, upper_limit = Inf, subdivisions = 100) {
#   # Vectorize integration over parameters
#   mapply(function(nz, no, sz, so, n) {
#     # Define the integrand for a single set of parameters
#     integrand <- function(rt_vec) {
#       dlnr(rt_vec, nuzero = nz, nuone = no, sigmazero = sz, sigmaone = so, ndt = n, response = response)
#     }
#     # Integrate from ndt to upper_limit
#     lower_bound <- n + .Machine$double.eps # Start integration just above ndt
#     if (lower_bound >= upper_limit) return(0) # Handle case where ndt >= upper_limit
#
#     prob <- tryCatch(
#       stats::integrate(integrand, lower = lower_bound, upper = upper_limit,
#                        subdivisions = subdivisions, stop.on.error = FALSE)$value,
#       error = function(e) {
#         warning("Integration failed in .marginal_p_lnr: ", e$message)
#         NA_real_ # Return NA on integration failure
#       }
#     )
#     # Clamp probability between 0 and 1 to handle potential integration errors
#     max(0, min(1, prob, na.rm = TRUE), na.rm = TRUE) # Added na.rm
#   }, nuzero, nuone, sigmazero, sigmaone, ndt, SIMPLIFY = TRUE)
# }
#
# # --- Commented out epred implementation using the helper ---
# #' @rdname rlnr
# #' @export
# #' @section Posterior Epred:
# #'   Note: Calculating the posterior expected prediction for the LNR model
# #'   typically involves numerical integration for each posterior draw and
# #'   observation, which can be computationally intensive and slow. This
# #'   function calculates the marginal probability of response 1, P(response=1).
# #'   Calculating the expected reaction time E[RT] is generally too complex
# #'   for this context. Use with caution, especially for large datasets or
# #'   many posterior draws. Consider calculating expectations manually post-hoc
# #'   if needed.
# posterior_epred_lnr <- function(prep) {
#   # Extract parameters for ALL observations (matrices: draws x observations)
#   nuzero <- brms::get_dpar(prep, "nuzero")
#   nuone <- brms::get_dpar(prep, "nuone")
#   sigmazero <- brms::get_dpar(prep, "sigmazero")
#   sigmaone <- brms::get_dpar(prep, "sigmaone")
#   tau <- brms::get_dpar(prep, "tau")
#   minrt <- brms::get_dpar(prep, "minrt")
#
#   # Calculate ndt (matrix: draws x observations)
#   ndt <- tau * minrt
#
#   # Determine a reasonable upper integration limit.
#   upper_limit <- Inf # Using Inf, accepting potential slowness/integration issues.
#
#   # Get dimensions
#   n_draws <- nrow(nuzero)
#   n_obs <- ncol(nuzero)
#   epred_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)
#
#   # Loop through observations (columns)
#   for (j in 1:n_obs) {
#       # Pass columns (vectors of parameters for this observation across draws)
#       epred_matrix[, j] <- .posterior_epred_lnr_marginal_p(
#           nuzero = nuzero[, j],
#           nuone = nuone[, j],
#           sigmazero = sigmazero[, j],
#           sigmaone = sigmaone[, j],
#           ndt = ndt[, j],
#           response = 1, # Calculate P(response=1)
#           upper_limit = upper_limit
#       )
#   }
#
#   # Return the matrix of probabilities (draws x observations)
#   epred_matrix
# }