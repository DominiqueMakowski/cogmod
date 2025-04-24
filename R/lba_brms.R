# Stanvars ----------------------------------------------------------------

#' @keywords internal
.lba_lpdf <- function() {
"
// Log probability density function for the LBA model using the new parametrization.
// Y: observed reaction time.
// mu: mean drift for accumulator 0.
// driftone: mean drift for accumulator 1.
// sigmazero, sigmaone: standard deviations of drift for accumulators 0 and 1.
// sigmabias: starting-point range, A = sigmabias.
// bs: threshold offset such that b = sigmabias + bs.
// tau: scale factor for non-decision time (0-1, scaled by minrt).
// minrt: minimum possible reaction time (used to scale tau).
// dec: decision indicator (0 if accumulator 0 wins; 1 if accumulator 1 wins).
real lba_lpdf(real Y,
              real mu,
              real driftone,
              real sigmazero,
              real sigmaone,
              real sigmabias,
              real bs,
              real tau,
              real minrt,
              int dec) {
  // Derived parameters.
  real A = sigmabias;             // Starting-point range.
  real b = sigmabias + bs;        // Decision threshold.
  real ndt = tau * minrt;         // Non-decision time.
  
  // Check if the observed RT is less than or equal to ndt.
  if (Y <= ndt) return negative_infinity();
  
  // Now compute the adjusted decision time.
  real t_adj = Y - ndt;
  
  // --- Early Parameter Validity Checks ---
  if (sigmazero <= 0 || sigmaone <= 0) return negative_infinity();
  if (dec != 0 && dec != 1) return negative_infinity();

  // --- Standard Case: Both drifts are finite.
  // Define winning and losing accumulator parameters based on dec.
  real v_win;
  real s_win;
  real v_loss;
  real s_loss;
  if (dec == 0) {
    v_win  = mu;
    s_win  = sigmazero;
    v_loss = driftone;
    s_loss = sigmaone;
  } else {  // dec == 1.
    v_win  = driftone;
    s_win  = sigmaone;
    v_loss = mu;
    s_loss = sigmazero;
  }
  
  // --- Compute Defective Density for the Winning Accumulator ---
  real st_win = s_win * t_adj;
  st_win = fmax(st_win, 1e-10);
  real z1_win = (b - A - v_win * t_adj) / st_win;
  real z2_win = (b - v_win * t_adj) / st_win;
  real phi_z1_win = exp(-0.5 * square(z1_win)) / sqrt(2 * pi());
  real phi_z2_win = exp(-0.5 * square(z2_win)) / sqrt(2 * pi());
  real cdf_z1_win = Phi(z1_win);
  real cdf_z2_win = Phi(z2_win);
  real f_win = (1 / A) * ( - v_win * cdf_z1_win + s_win * phi_z1_win
                           + v_win * cdf_z2_win - s_win * phi_z2_win );
  f_win = fmax(f_win, 1e-10);
  real log_f_win = log(f_win);
  
  // --- Compute the Cumulative Distribution for the Losing Accumulator ---
  real st_loss = s_loss * t_adj;
  st_loss = fmax(st_loss, 1e-10);
  real z1_loss = (b - A - v_loss * t_adj) / st_loss;
  real z2_loss = (b - v_loss * t_adj) / st_loss;
  real cdf_z1_loss = Phi(z1_loss);
  real cdf_z2_loss = Phi(z2_loss);
  real phi_z1_loss = exp(-0.5 * square(z1_loss)) / sqrt(2 * pi());
  real phi_z2_loss = exp(-0.5 * square(z2_loss)) / sqrt(2 * pi());
  real term_loss = ((b - A - v_loss * t_adj) / A) * cdf_z1_loss
                   - ((b - v_loss * t_adj) / A) * cdf_z2_loss
                   + (st_loss / A) * (phi_z1_loss - phi_z2_loss);
  real F_loss = 1 + term_loss;
  F_loss = fmin(fmax(F_loss, 0), 1 - 1e-10);
  real log_surv = log1m(F_loss);  // log(1 - F_loss)
  
  // --- Return the Log Likelihood ---
  real log_lik = log_f_win + log_surv;
  if (is_nan(log_lik) || is_inf(log_lik)) return negative_infinity();
  return log_lik;
}
"
}



#' @rdname rlba
#' @examples
#' # You can expose the lpdf function as follows:
#' # lba_lpdf <- lba_lpdf_expose()
#' # lba_lpdf(...)
#'
#' @export
lba_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  
  # This wraps the (new) Stan function lba_lpdf using our new parametrization.
  stancode <- paste0(
"functions {
", .lba_lpdf(), "
}")
  
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$lba_lpdf
}

#' @rdname rlba
#' @export
lba_stanvars <- function() {
  brms::stanvar(scode = .lba_lpdf(), block = "functions")
}


#' @rdname rlba
#' @param link_mu,link_driftone,link_sigmazero,link_sigmaone,link_sigmabias,link_bs,link_tau,link_minrt Link functions for the parameters.
#' @export
lba <- function(link_mu = "identity",
                link_driftone = "identity",
                link_sigmazero = "softplus",
                link_sigmaone = "softplus",
                link_sigmabias = "softplus",
                link_bs = "softplus",
                link_tau = "logit",
                link_minrt = "identity") {
  brms::custom_family(
    name = "lba",
    dpars = c("mu", "driftone", "sigmazero", "sigmaone", "sigmabias", "bs", "tau", "minrt"),
    links = c(link_mu, link_driftone, link_sigmazero, link_sigmaone, link_sigmabias, link_bs, link_tau, link_minrt),
    lb = c(NA, NA, 0, 0, 0, 0, 0, 0),
    ub = c(NA, NA, NA, NA, NA, NA, 1, NA),
    type = "real",
    vars = c("dec[n]")
  )
}

# brms Post-processing Functions ------------------------------------------

#' @rdname rlba
#' @inheritParams rlnr
#' @export
log_lik_lba <- function(i, prep) {
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

#' @rdname rlba
#' @export
posterior_predict_lba <- function(i, prep, ...) {
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

#' @rdname rlba
#' @export
posterior_epred_lba <- function(prep) {
  stop(
    "Calculating the posterior expected prediction (epred) for the LBA model ",
    "is computationally prohibitive within this framework.\n",
    "Please use `posterior_predict()` to obtain draws from the posterior ",
    "predictive distribution and calculate summaries manually if needed."
  )
}