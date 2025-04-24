# Stanvars ----------------------------------------------------------------

#' @keywords internal
.rt_lba_lpdf <- function() {
"
// Log density function for the one-accumulator LBA model.
// Y: observed reaction time.
// mu: mean drift rate of the accumulator (named 'mu' because of brms requirements).
// sigma: standard deviation of the drift.
// sigmabias: the range of the starting point (A = sigmabias).
// bs: additive component so that the decision threshold is b = sigmabias + bs.
// tau: scale factor for non-decision time (scaled by minrt).
// minrt: minimum possible reaction time (used to scale tau).
real rt_lba_lpdf(real Y, real mu, real sigma, real sigmabias, real bs, real tau, real minrt) {
  
  // Derived parameters.
  real A = sigmabias;         // Uniform start-point range: U(0, A)
  real b = A + bs;            // Decision threshold.
  real ndt = tau * minrt;  // Non-decision time.
  real eps = 1e-10; // Small value for numerical stability.
  
  // If the observed RT is less than or equal to non-decision time, no decision is made.
  if (Y <= ndt) return negative_infinity();
  
  // Adjusted (decision) time.
  real T = Y - ndt;
  
  // Compute st = sigma * T; use a small floor for stability.
  real st = sigma * T;
  st = fmax(st, eps);
  
  // Compute z-values.
  // z1 corresponds to lower limit (based on threshold b - A),
  // z2 corresponds to upper limit (based on threshold b).
  real z1 = (b - A - mu * T) / st;
  real z2 = (b - mu * T) / st;
  
  // Compute the standard normal PDF at the z-values.
  // (Stan does not have a built-in normal_pdf; this is done via exp(normal_lpdf(...)).)
  real phi_z1 = exp(normal_lpdf(z1 | 0, 1));
  real phi_z2 = exp(normal_lpdf(z2 | 0, 1));
  
  // Compute the defective density using the conventional formula:
  real f_val = (1 / A) * ( mu * (Phi(z2) - Phi(z1))
                           + sigma * (phi_z1 - phi_z2) );
  
  // Normalize by the probability that the drift is positive.
  real prob_positive_drift = 1 - Phi((0 - mu) / sigma);
  f_val = f_val / prob_positive_drift;
  
  // Floor the density for numerical stability.
  f_val = fmax(f_val, eps);
  
  return log(f_val);
}
"
}



#' @rdname rrt_lba
#' @export
rt_lba_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  
  # This wraps the (new) Stan function rt_lba_lpdf using our new parametrization.
  stancode <- paste0(
"functions {
", .rt_lba_lpdf(), "
}")
  
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_lba_lpdf
}

#' @rdname rrt_lba
#' @export
rt_lba_stanvars <- function() {
  brms::stanvar(scode = .rt_lba_lpdf(), block = "functions")
}


#' @rdname rrt_lba
#' @param link_mu,link_sigma,link_sigmabias,link_bs,link_tau,link_minrt Link functions for the parameters.
#' @export
rt_lba <- function(link_mu = "identity",
                   link_sigma = "softplus",
                   link_sigmabias = "softplus",
                   link_bs = "softplus",
                   link_tau = "logit",
                   link_minrt = "identity") {
  brms::custom_family(
    name = "rt_lba",
    dpars = c("mu", "sigma", "sigmabias", "bs", "tau", "minrt"),
    links = c(link_mu, link_sigma, link_sigmabias, link_bs, link_tau, link_minrt),
    lb = c(NA, 0, 0, 0, 0, 0),
    ub = c(NA, NA, NA, NA, 1, NA),
    type = "real"
  )
}

# brms Post-processing Functions ------------------------------------------

#' @rdname rrt_lba
#' @inheritParams rlnr
#' @export
log_lik_rt_lba <- function(i, prep) {
  # Extract observed reaction time for observation i
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_) # Handle missing RTs

  # Get the parameter draws for observation i:
  driftzero <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)  # Minimum possible reaction time

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Get decision indicator for observation i (should be 0 or 1).
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning("Response (dec) must be 0 or 1. Got: ", response, " for observation ", i)
    return(rep(-Inf, length(driftzero)))
  }

  # Compute log-likelihood using our R density function 'dlba'
  ll <- dlba(x = y, response = response,
             driftzero = driftzero, driftone = driftone,
             sigmazero = sigmazero, sigmaone = sigmaone,
             sigmabias = sigmabias, bs = bs, ndt = ndt,
             log = TRUE)
  ll  # Return vector of log-likelihoods (one per draw)
}

#' @rdname rrt_lba
#' @export
posterior_predict_rt_lba <- function(i, prep, ...) {
  # Get the parameter draws for observation i:
  driftzero <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt
  
  # Number of posterior draws:
  n_draws <- length(driftzero)
  
  # Generate predictions using the R simulation function 'rlba'
  # It is assumed that 'rlba' is vectorized over parameters.
  sim_data <- rlba(n = n_draws,
                   driftzero = driftzero,
                   driftone = driftone,
                   sigmazero = sigmazero,
                   sigmaone = sigmaone,
                   sigmabias = sigmabias,
                   bs = bs,
                   ndt = ndt)
  as.matrix(sim_data)  
}

#' @rdname rrt_lba
#' @export
posterior_epred_rt_lba <- function(prep) {
  stop(
    "Calculating the posterior expected prediction (epred) for the LBA model ",
    "is computationally prohibitive within this framework.\n",
    "Please use `posterior_predict()` to obtain draws from the posterior ",
    "predictive distribution and calculate summaries manually if needed."
  )
}