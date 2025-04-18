# Stan Functions and Custom Family for brms (Beta-Gate)

# Stanvars ----------------------------------------------------------------

#' @keywords internal
.betagate_lpdf <- function() {
"
// Log probability density function for the Beta-Gate distribution
real betagate_lpdf(real y, real mu, real phi, real pex, real bex) {
  // Tolerance for floating point comparisons near 0 and 1
  real eps = 1e-10;

  // --- Parameter Validation ---
  // Ensure parameters are within valid ranges. Note: brms often handles this
  // via link functions, but explicit checks add robustness.
  if (!(mu > 0.0 && mu < 1.0) || !(phi > 0.0) ||
      !(pex >= 0.0 && pex <= 1.0) || !(bex >= 0.0 && bex <= 1.0) ||
      !(y >= 0.0 && y <= 1.0)) {
    return negative_infinity();
  }

  // --- Calculate Cutpoints ---
  // Calculate probability-scale cutpoints and apply logit scale
  real cutzerolog = logit(pex * (1.0 - bex));
  real cutonelog = logit(1.0 - pex * bex);

  // --- Calculate Log Probability based on y ---
  // Location parameter on logit scale
  real mu_ql = logit(mu);

  if (abs(y - 0.0) < eps) { // Case: y = 0
    // Log probability P(latent <= cutzerolog)
    return log1m_inv_logit(mu_ql - cutzerolog);

  } else if (abs(y - 1.0) < eps) { // Case: y = 1
    // Log probability P(latent > cutonelog)
    return log_inv_logit(mu_ql - cutonelog);

  } else { // Case: 0 < y < 1
    // Log probability P(cutzerolog < latent <= cutonelog)
    real log_prob_middle = log_diff_exp(log_inv_logit(mu_ql - cutzerolog),
                                        log_inv_logit(mu_ql - cutonelog));

    // Beta distribution parameters
    real shape1 = mu * phi * 2.0;
    real shape2 = (1.0 - mu) * phi * 2.0;

    // Log Beta density for y
    real log_beta_dens = beta_lpdf(y | shape1, shape2);

    // Total log probability: log(P(middle)) + log(BetaPDF(y))
    return log_prob_middle + log_beta_dens;
  }
}
"
}


#' @rdname rbetagate
#' @examples
#' # You can expose the lpdf function as follows:
#' # betagate_lpdf <- betagate_lpdf_expose()
#' # betagate_lpdf(y = 0.5, mu = 0.6, phi = 10, pex = 0.2, bex = 0.5)
#'
#' @export
betagate_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .betagate_lpdf(),
    "\n}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$betagate_lpdf
}

#' @rdname rbetagate
#' @export
betagate_stanvars <- function() {
  brms::stanvar(scode = .betagate_lpdf(), block = "functions")
}


#' @rdname rbetagate
#' @param link_mu,link_phi,link_pex,link_bex Link functions for the parameters.
#' @export
#' @examples
#' \dontrun{
#' # Default usage:
#' family = betagate()
#' # Custom link for phi:
#' family = betagate(link_phi = "log")
#' }
betagate <- function(link_mu = "logit", link_phi = "softplus", link_pex = "logit", link_bex = "logit") {
  brms::custom_family(
    name = "betagate",
    dpars = c("mu", "phi", "pex", "bex"),
    links = c(link_mu, link_phi, link_pex, link_bex),
    lb = c(NA, 0, 0, 0), # Lower bounds: phi>0, pex>=0, bex>=0
    ub = c(NA, NA, 1, 1), # Upper bounds: pex<=1, bex<=1
    type = "real"       # Outcome variable type
  )
}

# brms Post-processing Functions -------------------------------------------

#' @rdname rbetagate
#' @export
log_lik_betagate <- function(i, prep) {
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

  # Calculate log-likelihood using the vectorized dbetagate function
  ll <- dbetagate(x = y_vec, mu = mu, phi = phi, pex = pex, bex = bex, log = TRUE)

  # Ensure no NaN/NA values (dbetagate should return -Inf for zero density)
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll # Return the vector of log-likelihoods for all draws
}


#' @rdname rbetagate
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_betagate <- function(i, prep, ...) {
  # Extract draws for each parameter for the i-th observation
  mu    <- brms::get_dpar(prep, "mu", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  pex   <- brms::get_dpar(prep, "pex", i = i)
  bex   <- brms::get_dpar(prep, "bex", i = i)

  n_draws <- length(mu) # Number of posterior draws
  if (n_draws == 0) return(matrix(numeric(0), ncol = 1)) # Handle case with no draws

  # Use rbetagate() to generate predictions for each draw
  final_out <- rbetagate(n = n_draws, mu = mu, phi = phi, pex = pex, bex = bex)

  # Return as a matrix with one column (standard for posterior_predict)
  as.matrix(final_out)
}


#' @rdname rbetagate
#' @export
posterior_epred_betagate <- function(prep) {
  # Extract draws for the necessary parameters (draws x observations matrices)
  mu  <- brms::get_dpar(prep, "mu")
  phi <- brms::get_dpar(prep, "phi") # phi is not directly needed for epred, but kept for consistency
  pex <- brms::get_dpar(prep, "pex")
  bex <- brms::get_dpar(prep, "bex")

  # Calculate the expectation (mean) based on the Ordered Beta logic
  # E[Y] = P(Y=1) * 1 + P(0 < Y < 1) * E[Y | 0 < Y < 1]
  # E[Y | 0 < Y < 1] is simply the mean of the Beta distribution, which is mu.
  # E[Y] = P(Y=1) + P(0 < Y < 1) * mu

  # --- Calculate Probabilities (Vectorized) ---
  # Cutpoints on probability scale
  cutzero_p <- pex * (1 - bex)
  cutone_p <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero_p)
  cutonelog <- stats::qlogis(cutone_p)

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # Probabilities for the three categories
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE) # P(eta <= cutzerolog) -> Y=0
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE) # P(eta > cutonelog) -> Y=1
  prob_mid <- 1 - prob_0 - prob_1 # P(cutzerolog < eta <= cutonelog) -> Y ~ Beta(mu, phi)
  # Ensure probabilities sum to 1 and handle potential floating point inaccuracies
  prob_mid <- pmax(0, prob_mid)

  # --- Calculate Expectation ---
  # E[Y] = P(Y=1) * 1 + P(0 < Y < 1) * E[Y | 0 < Y < 1]
  # E[Y] = prob_1 * 1 + prob_mid * mu
  epred <- prob_1 + prob_mid * mu

  # Handle edge cases where parameters might lead to NaN (e.g., invalid mu)
  epred[is.nan(epred) | is.na(epred)] <- NA # Or a default value if preferred

  epred # Return the matrix of posterior expectations (draws x observations)
}