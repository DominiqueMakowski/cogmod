# Stanvars ----------------------------------------------------------------

# Stan Functions and Custom Family for brms

#' @keywords internal
.bext_lpdf <- function() {
"
// Log probability density function for the BEXT distribution
// Optimized for numerical stability
real bext_lpdf(real y, real mu, real phi, real pex, real bex) {
  real eps = 1e-9; // Tolerance for floating point comparisons near 0 and 1

  // --- Parameter Validation ---
  if (mu < 0.0 || mu > 1.0 || phi <= 0.0 || pex < 0.0 || pex > 1.0 || bex < 0.0 || bex > 1.0) {
    return negative_infinity();
  }
  // Validate outcome y
  if (y < 0.0 || y > 1.0) {
     return negative_infinity();
  }

  // Calculate internal precision
  real precision = phi * 2.0;

  // --- Edge Case: pex = 1 (Pure Bernoulli) ---
  // If pex is exactly or very close to 1
  if (abs(pex - 1.0) < eps) {
    if (abs(y - 0.0) < eps) { // y = 0 (using abs(y-0.0) for consistency)
      // log(1 - bex)
      if (bex >= 1.0 - eps) return negative_infinity(); // Avoid log(0) if bex=1
      return log1m(bex); // More stable than log(1.0 - bex)
    } else if (abs(y - 1.0) < eps) { // y = 1  <- CORRECTED CONDITION
      // log(bex)
      if (bex <= eps) return negative_infinity(); // Avoid log(0) if bex=0
      return log(bex);
    } else { // 0 < y < 1
      // Density is zero for continuous values in Bernoulli
      return negative_infinity();
    }
  }

  // --- Edge Case: pex = 0 (Pure Beta) ---
  // If pex is exactly or very close to 0
  if (abs(pex - 0.0) < eps) {
     // For pure Beta, density is defined only on (0, 1)
     if (y <= 0.0 + eps || y >= 1.0 - eps) {
         // Stan's beta_lpdf might handle this, but explicit check is safer
         return negative_infinity();
     } else {
         // Use Stan's built-in beta log PDF
         return beta_lpdf(y | mu * precision, (1.0 - mu) * precision);
     }
  }

  // --- Mixture Case (0 < pex < 1) ---
  // At this point, we know 0 < pex < 1

  if (abs(y - 0.0) < eps) { // y = 0
    // Log-probability is log(pex * (1 - bex)) = log(pex) + log(1 - bex)
    if (bex >= 1.0 - eps) return negative_infinity(); // Avoid log(1-bex) if bex=1
    // pex > 0 is guaranteed here
    return log(pex) + log1m(bex);
  } else if (abs(y - 1.0) < eps) { // y = 1
    // Log-probability is log(pex * bex) = log(pex) + log(bex)
    if (bex <= eps) return negative_infinity(); // Avoid log(bex) if bex=0
    // pex > 0 is guaranteed here
    return log(pex) + log(bex);
  } else { // 0 < y < 1
    // Log-density is log(1 - pex) + beta_lpdf(...)
    // pex < 1 is guaranteed here, so (1-pex) > 0
    return log1m(pex) + beta_lpdf(y | mu * precision, (1.0 - mu) * precision);
  }
}
"
}


#' @rdname rbext
#' @examples
#' # You can expose the lpdf function as follows:
#' # bext_lpdf <- bext_lpdf_expose()
#' # bext_lpdf(y = 0.5, mu = 0.6, phi = 10, pex = 0.2, bex = 0.5)
#'
#' @export
bext_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  
  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .bext_lpdf(),
    "\n}"
  )
  
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$bext_lpdf
}

#' @rdname rbext
#' @export
bext_stanvars <- function() {
  brms::stanvar(scode = .bext_lpdf(), block = "functions")
}


#' @rdname rbext
#' @param link_mu,link_phi,link_pex,link_bex Link functions for the parameters.
#' @export
#' @examples
#' \dontrun{
#' # Default usage:
#' family = bext()
#' # Custom link for phi:
#' family = bext(link_phi = "log")
#' }
bext <- function(link_mu = "logit", link_phi = "softplus", link_pex = "logit", link_bex = "logit") {
  brms::custom_family(
    name = "bext",
    dpars = c("mu", "phi", "pex", "bex"),
    links = c(link_mu, link_phi, link_pex, link_bex),
    lb = c(NA, 0, 0, 0),
    ub = c(NA, NA, 1, 1),  # Question: are lb and ub needed for pex and bex?
    type = "real"
  )
}

# brms --------------------------------------------------------------------

# Posterior Prediction

#' @rdname rbext
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_bext <- function(i, prep, ...) {
  # Extract draws for each parameter
  mu    <- brms::get_dpar(prep, "mu", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  pex   <- brms::get_dpar(prep, "pex", i = i)
  bex   <- brms::get_dpar(prep, "bex", i = i)

  n_draws <- length(mu) # Number of posterior draws

  # Use rbext()
  final_out <- rbext(n = n_draws, mu = mu, phi = phi, pex = pex, bex = bex)

  as.matrix(final_out)
}


#' @rdname rbext
#' @export
posterior_epred_bext <- function(prep) {
  # Extract draws for the necessary parameters
  # get_dpar returns draws x observations matrix
  mu  <- brms::get_dpar(prep, "mu")
  pex <- brms::get_dpar(prep, "pex")
  bex <- brms::get_dpar(prep, "bex")

  # Calculate the expectation (mean) for each draw and observation
  # The expectation of BEXT(mu, phi, pex, bex) is (1 - pex) * mu + pex * bex.
  # Formula: E[Y] = (1 - pex) * E[Beta] + pex * E[Bernoulli]
  # E[Beta] = mu
  # E[Bernoulli] = bex
  epred <- (1 - pex) * mu + pex * bex

  epred # Return the matrix of posterior expectations (draws x observations)
}



#' @rdname rbext
#' @export
log_lik_bext <- function(i, prep) {
  # Extract observed value for the i-th observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y_scalar <- prep$data$Y[i] # This is a single value

  # Extract model draws (vectors) for the i-th observation
  mu  <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  # Determine number of draws
  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0)) # Handle case with no draws

  # Replicate the scalar y to match the number of draws
  y_vec <- rep(y_scalar, length.out = n_draws)

  # Calculate log-likelihood using the vectorized dbext function
  # Now y_vec has the same length as the parameter vectors
  ll <- dbext(x = y_vec, mu = mu, phi = phi, pex = pex, bex = bex, log = TRUE)

  # Ensure no NaN/NA values (dbext should return -Inf for zero density)
  # This might be redundant if dbext is robust, but safe to keep.
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll # Return the vector of log-likelihoods for all draws
}