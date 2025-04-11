# Stanvars ----------------------------------------------------------------

#' @keywords internal
.lnr_lpdf <- function() {
"
// Log-likelihood for a single observation from the reparameterized Log-Normal Race model.
// Y: observed reaction time.
// dec: decision indicator (0 or 1).
// mu: baseline accumulator mean (in log-space) for choice 0.
// mudelta: log-deviation for the mean of choice 1.
// sigmazero: baseline accumulator standard deviation (in log-space) for choice 0.
// sigmadelta: log-deviation for the standard deviation of choice 1.
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (used to scale tau).
real lnr_lpdf(real Y, real mu, real mudelta, real sigmazero, real sigmadelta, real tau, real minrt, int dec) {
  real eps = 1e-8; // Use a slightly smaller epsilon

  // --- 1. Input Parameter Validation ---
  // Basic checks for plausible parameter values enforced by model/priors,
  // but good practice for standalone function robustness.
  // Check for non-finite values using is_inf() and is_nan()
  if (sigmazero <= 0.0 || tau < 0.0 || tau > 1.0 || minrt < 0.0 ||
      (is_inf(mu) || is_nan(mu)) ||
      (is_inf(mudelta) || is_nan(mudelta)) ||
      (is_inf(sigmazero) || is_nan(sigmazero)) || // Also check sigmazero itself
      (is_inf(sigmadelta) || is_nan(sigmadelta)) ||
      (is_inf(tau) || is_nan(tau)) ||
      (is_inf(minrt) || is_nan(minrt))) {
    return negative_infinity();
  }
   // Check decision validity
  if (dec != 0 && dec != 1) {
    return negative_infinity();
  }

  // --- 2. Calculate Derived Parameters ---
  // Ensure positive standard deviations
  real sig0 = fmax(sigmazero, eps);
  // Calculate sigma1 carefully
  real exp_sigmadelta = exp(sigmadelta);
  // Check for potential overflow/NaN from exp() using is_inf() and is_nan()
  if (is_inf(exp_sigmadelta) || is_nan(exp_sigmadelta)) return negative_infinity();
  real sig1 = fmax(sig0 * exp_sigmadelta, eps);

  // Calculate non-decision time
  real ndt = tau * minrt;

  // Calculate adjusted RT (time available for the race process)
  real rt_adj = Y - ndt;

  // Check if adjusted RT is valid (must be positive)
  if (rt_adj < eps) {
    return negative_infinity();
  }

  // --- 3. Calculate Log-Likelihood Components ---
  real mu0 = mu;
  real mu1 = mu + mudelta;
  real log_lik;

  // --- 4. Calculate log-likelihood ---
  // Case 1: Response 0 (dec == 0) ---
  if (dec == 0) {
    // Accumulator 0 finished first
    // Component 1: Likelihood of accumulator 0 finishing at rt_adj
    real log_pdf0 = lognormal_lpdf(rt_adj | mu0, sig0);

    // Component 2: Likelihood of accumulator 1 finishing *after* rt_adj (survival function)
    real log_cdf1 = lognormal_lcdf(rt_adj | mu1, sig1);

    // Check for impossible scenario (CDF1 is ~1, meaning acc 1 finished before rt_adj)
    if (log_cdf1 > -eps) { // Use eps for consistency
        return negative_infinity();
    }

    // Calculate log survival probability for accumulator 1 stably
    real log_surv1;
    // If log_cdf1 is very small, log(1-cdf1) is approx log(1) = 0
    // log(1 - exp(-37)) is approx 0 within double precision
    if (log_cdf1 < -37.0) {
        log_surv1 = 0.0;
    } else {
        log_surv1 = log1m_exp(log_cdf1); // Use log1m_exp for stability: log(1 - CDF1)
    }

    log_lik = log_pdf0 + log_surv1;

    // Case 2: Response 1 (dec == 1) ---
  } else {
    // Accumulator 1 finished first
    // Component 1: Likelihood of accumulator 1 finishing at rt_adj
    real log_pdf1 = lognormal_lpdf(rt_adj | mu1, sig1);

    // Component 2: Likelihood of accumulator 0 finishing *after* rt_adj (survival function)
    real log_cdf0 = lognormal_lcdf(rt_adj | mu0, sig0);

    // Check for impossible scenario (CDF0 is ~1)
     if (log_cdf0 > -eps) { // Use eps for consistency
        return negative_infinity();
    }

    // Calculate log survival probability for accumulator 0 stably
    real log_surv0;
    if (log_cdf0 < -37.0) {
        log_surv0 = 0.0;
    } else {
        log_surv0 = log1m_exp(log_cdf0); // log(1 - CDF0)
    }

    log_lik = log_pdf1 + log_surv0;
  }

  // Final check for non-finite results (e.g., from lognormal_lpdf itself) using is_inf() and is_nan()
  if (is_inf(log_lik) || is_nan(log_lik)) {
      return negative_infinity();
  }

  return log_lik;
}
"
}

#' @rdname rlnr
#' @examples
#' # You can expose the lpdf function as follows:
#' # lnr_lpdf <- lnr_lpdf_expose()
#' # lnr_lpdf(Y = 0.5, mu = 0, mudelta = 0, sigmazero = 1, sigmadelta = 0.5,
#' #          tau = 0.1, minrt = 0.2, dec = 1)
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
#' @param link_mu Link function for the `mu` parameter in the custom family.
#'   Determines how `mu` is transformed in the model. Default: "identity".
#' @param link_mudelta Link function for the `mudelta` parameter in the custom family.
#'   Determines how `mudelta` is transformed in the model. Default: "identity".
#' @param link_sigmazero Link function for the `sigmazero` parameter in the custom family.
#'   Ensures `sigmazero` remains positive. Default: "identity".
#' @param link_sigmadelta Link function for the `sigmadelta` parameter in the custom family.
#'   Determines how `sigmadelta` is transformed in the model. Default: "identity".
#' @param link_tau Link function for the `tau` parameter in the custom family.
#' @param link_minrt Link function for the `minrt` parameter in the custom family.
#' @export
lnr <- function(link_mu = "identity", link_mudelta = "identity",
                link_sigmazero = "softplus", link_sigmadelta = "identity",
                link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "lnr",
    dpars = c("mu", "mudelta", "sigmazero", "sigmadelta", "tau", "minrt"),
    links = c(link_mu, link_mudelta, link_sigmazero, link_sigmadelta, link_tau, link_minrt),
    lb = c(NA, NA, 0, NA, 0, 0),
    ub = c(NA, NA, NA, NA, 1, NA),    # Ensure tau stays <= 1
    vars = "dec[n]"
  )
}



# # brms --------------------------------------------------------------------

#' @rdname rlnr
log_lik_lnr <- function(i, prep) {
  # Extract observation
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA)

  # Get parameters
  mu <- brms::get_dpar(prep, "mu", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time
  ndt <- tau * minrt

  # Get decision
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning("Response must be 0 or 1. Got: ", response)
    return(rep(-Inf, length(mu)))
  }

  # Compute log-likelihood
  ll <- dlnr(x = y, mu = mu, mudelta = mudelta,
           sigmazero = sigmazero, sigmadelta = sigmadelta,
           ndt = ndt, response = response, log = TRUE)

  ll
}



#' @rdname rlnr
#' @inheritParams choco
#' @export
posterior_predict_lnr <- function(i, prep, ...) {
  # Get parameters
  mu <- brms::get_dpar(prep, "mu", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)  # Get minrt parameter

  # Calculate non-decision time
  ndt <- tau * minrt

  # Generate predictions
  n_draws <- length(mu)
  sim_data <- rlnr(n_draws, mu = mu, mudelta = mudelta,
                  sigmazero = sigmazero, sigmadelta = sigmadelta,
                  ndt = ndt)

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
       " and the user is better off computing 'predictions' rather than 'expectations' and", 
       " using the posterior draws to compute any quantities of interest.")
}






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
#   mu <- brms::get_dpar(prep, "mu")
#   mudelta <- brms::get_dpar(prep, "mudelta")
#   sigmazero <- brms::get_dpar(prep, "sigmazero")
#   sigmadelta <- brms::get_dpar(prep, "sigmadelta")
#   tau <- brms::get_dpar(prep, "tau")
#   minrt <- brms::get_dpar(prep, "minrt")

#   # Calculate ndt (matrix: draws x observations)
#   ndt <- tau * minrt

#   # Determine a reasonable upper integration limit.
#   # Using a fixed large value or Inf is okay, but might be slow.
#   # Alternatively, estimate from observed data quantiles if available,
#   # but that's tricky as it should reflect the *predicted* range.
#   # Let's use Inf for now, accepting potential slowness/integration issues.
#   upper_limit <- Inf

#   # Calculate P(response=1) for each draw and observation
#   # Apply the helper function. This will be slow.
#   # Need to handle the matrix structure appropriately.
#   # .posterior_epred_lnr_marginal_p expects vectors, so we iterate or use apply.
#   # Using apply might be cleaner but potentially not faster.

#   # Get dimensions
#   n_draws <- nrow(mu)
#   n_obs <- ncol(mu)
#   epred_matrix <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

#   # Loop through observations (columns) - often more efficient in R than rows
#   for (j in 1:n_obs) {
#       # Pass columns (vectors of parameters for this observation across draws)
#       epred_matrix[, j] <- .posterior_epred_lnr_marginal_p(
#           mu = mu[, j],
#           mudelta = mudelta[, j],
#           sigmazero = sigmazero[, j],
#           sigmadelta = sigmadelta[, j],
#           ndt = ndt[, j],
#           response = 1, # Calculate P(response=1)
#           upper_limit = upper_limit
#       )
#   }

#   # Return the matrix of probabilities (draws x observations)
#   epred_matrix
# }



# #' Calculate Marginal Probability for LNR (Internal Helper)
# #'
# #' Integrates the joint density dlnr to get the marginal probability of a response.
# #' NOTE: This involves numerical integration and can be slow.
# #' @keywords internal
# .posterior_epred_lnr_marginal_p <- function(mu, mudelta, sigmazero, sigmadelta, ndt, response, upper_limit = Inf, subdivisions = 100) {
#   # Vectorize integration over parameters
#   mapply(function(m, md, sz, sd, n) {
#     # Define the integrand for a single set of parameters
#     integrand <- function(rt_vec) {
#       dlnr(rt_vec, mu = m, mudelta = md, sigmazero = sz, sigmadelta = sd, ndt = n, response = response)
#     }
#     # Integrate from ndt to upper_limit
#     lower_bound <- n + .Machine$double.eps # Start integration just above ndt
#     if (lower_bound >= upper_limit) return(0) # Handle case where ndt >= upper_limit

#     prob <- tryCatch(
#       stats::integrate(integrand, lower = lower_bound, upper = upper_limit,
#                        subdivisions = subdivisions, stop.on.error = FALSE)$value,
#       error = function(e) {
#         warning("Integration failed in .marginal_p_lnr: ", e$message)
#         NA # Return NA on integration failure
#       }
#     )
#     # Clamp probability between 0 and 1 to handle potential integration errors
#     max(0, min(1, prob))
#   }, mu, mudelta, sigmazero, sigmadelta, ndt, SIMPLIFY = TRUE)
# }

