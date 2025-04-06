# Stanvars ----------------------------------------------------------------


#' @rdname rlnr
#' @export
lnr_stanvars <- function() {
  brms::stanvar(scode = "
  // Log-likelihood for a single observation from the reparameterized Log-Normal Race model.
  // y: observed reaction time.
  // dec: decision indicator (0 or 1).
  // mu: baseline accumulator mean (in log-space) for choice 0.
  // mudelta: additive deviation for the mean of choice 1.
  // sigmazero: baseline accumulator standard deviation (in log-space) for choice 0.
  // sigmadelta: log-deviation for the standard deviation of choice 1.
  // tau: non-decision time (shift).
  real lnr_lpdf(real y, real mu, real mudelta, real sigmazero, real sigmadelta, real tau, int dec) {
    real eps = 1e-6;  // A small constant to prevent underflow
    real t_adj = y - tau;

    // If t_adj is too small, return a very low log probability
    if (t_adj < eps)
      return negative_infinity();

    // Convert 0-based decision to 1-based indexing for Stan.
    int dec_stan = dec + 1;

    // Construct means and standard deviations for the two accumulators.
    vector[2] nu;
    vector[2] sigma;
    nu[1] = mu;
    nu[2] = mu + mudelta;
    sigma[1] = fmax(sigmazero, eps);
    sigma[2] = fmax(sigmazero * exp(sigmadelta), eps);

    real lp = 0;
    // Sum contributions across both accumulators.
    for (i in 1:2) {
      if (i == dec_stan)
        lp += lognormal_lpdf(t_adj | nu[i], sigma[i]);
      else {
        // Use lognormal_lcdf, then compute the log survival probability safely.
        real log_cdf = lognormal_lcdf(t_adj | nu[i], sigma[i]);
        lp += log1m_exp(log_cdf);
      }
    }
    return lp;
  }
  ", block = "functions")
}





#' @rdname rlnr
#' @param link_mu Link function for the `mu` parameter in the custom family.
#'   Determines how `mu` is transformed in the model. Default: "identity".
#' @param link_mudelta Link function for the `mudelta` parameter in the custom family.
#'   Determines how `mudelta` is transformed in the model. Default: "identity".
#' @param link_sigmazero Link function for the `sigmazero` parameter in the custom family.
#'   Ensures `sigmazero` remains positive. Default: "softplus".
#' @param link_sigmadelta Link function for the `sigmadelta` parameter in the custom family.
#'   Determines how `sigmadelta` is transformed in the model. Default: "identity".
#' @param link_tau Link function for the `tau` parameter in the custom family.
#'   Ensures `tau` remains non-negative. Default: "softplus".
#' @export
lnr <- function(link_mu = "identity", link_mudelta = "identity",
                link_sigmazero = "softplus", link_sigmadelta = "identity",
                link_tau = "softplus") {
  brms::custom_family(
    name = "lnr",
    dpars = c("mu", "mudelta", "sigmazero", "sigmadelta", "tau"),  # Distributional parameters
    links = c(link_mu, link_mudelta, link_sigmazero, link_sigmadelta, link_tau),  # Link functions
    lb = c(NA, NA, 0, NA, 0),
    vars = "dec[n]"  # Additional variable for decision
  )
}

# brms --------------------------------------------------------------------


#' @rdname rlnr
#' @inheritParams choco
#' @export
posterior_predict_lnr <- function(i, prep, ...) {
  # Extract distributional parameters for draw i
  mu <- brms::get_dpar(prep, "mu", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)

  # Number of draws (here, each draw produces one simulated trial)
  n_draws <- length(tau)

  # Generate all predictions at once using the vectorized simulation function rlnr.
  sim_data <- rlnr(n_draws, mu = mu, mudelta = mudelta, sigmazero = sigmazero,
                   sigmadelta = sigmadelta, tau = tau)

  # Convert to matrix
  as.matrix(sim_data)
}




#' @rdname rlnr
log_lik_lnr <- function(i, prep) {
  # Extract the i-th observed reaction time
  y <- prep$data$Y[i]
  
  # Check if this is a valid observation (not NA)
  if (is.na(y)) {
    return(NA)
  }
  
  # Extract parameters
  mu        <- brms::get_dpar(prep, "mu", i = i)
  mudelta   <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  tau       <- brms::get_dpar(prep, "tau", i = i)
  
  # Extract the decision indicator (should be a scalar, 0 or 1)
  response <- prep$data[["dec"]][i]
  
  # Initialize result vector
  ll <- rep(NA_real_, length(mu))
  
  # Validate response
  if (!response %in% c(0, 1)) {
    warning("Response must be 0 or 1. Got: ", response)
    return(rep(-Inf, length(mu)))
  }
  
  # Basic parameter validation
  valid_idx <- which(sigmazero > 0 & tau >= 0 & y > tau)
  
  if (length(valid_idx) == 0) {
    return(rep(-Inf, length(mu)))
  }
  
  # For valid parameter sets, compute log-likelihood
  if (length(valid_idx) < length(mu)) {
    # Some invalid parameter sets exist
    ll[-valid_idx] <- -Inf
  }
  
  # Compute log likelihoods for valid parameter sets
  ll[valid_idx] <- dlnr(y = y, 
                      mu = mu[valid_idx], 
                      mudelta = mudelta[valid_idx], 
                      sigmazero = sigmazero[valid_idx], 
                      sigmadelta = sigmadelta[valid_idx], 
                      tau = tau[valid_idx], 
                      response = response, 
                      log = TRUE)
  
  ll
}