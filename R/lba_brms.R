# Stanvars ----------------------------------------------------------------

#' @keywords internal
.lba_lpdf <- function() {
"
// Log probability density function for the Linear Ballistic Accumulator (LBA) model.
// This version uses a relative parameterization for two choices and tau*minrt for ndt.
//
// Y Observed reaction time.
// mu Mean drift rate for accumulator 0.
// vdelta Additive deviation for the mean drift rate of accumulator 1 (v1 = vzero + vdelta).
// sigmazero Standard deviation of the drift rate for accumulator 0 (must be positive).
// sigmadelta Log-deviation for the standard deviation of accumulator 1 (sigma1 = sigmazero * exp(sigmadelta)).
// A Maximum start point (must be positive).
// k Threshold offset (b = A + k, k must be positive).
// tau Proportion of minrt for non-decision time (0 <= tau <= 1).
// minrt Minimum possible reaction time (must be non-negative).
// dec Decision indicator (0 or 1).
real lba_lpdf(real Y, real mu, real vdelta, real sigmazero, real sigmadelta,
              real A, real k, real tau, real minrt, int dec) {
  real eps = 1e-10; // Small epsilon for numerical stability
  
  // --- 1. Input validation ---
  if (sigmazero <= 0 || A <= 0 || k <= 0 || tau < 0 || tau > 1 || minrt < 0) {
    return negative_infinity();
  }
  if (dec != 0 && dec != 1) return negative_infinity();
  
  // --- 2. Derived Parameters ---
  real b = A + k;            // Decision threshold
  real v0 = mu;              // Drift rate for accumulator 0
  real v1 = mu + vdelta;     // Drift rate for accumulator 1
  real s0 = fmax(sigmazero, eps);  
  real s1 = fmax(s0 * exp(sigmadelta), eps);
  real ndt = tau * minrt;    // Non-decision time (derived)
  real t_adj = Y - ndt;
  if (t_adj < eps) return negative_infinity();
  
  // --- 3. Determine Winner/Loser Parameters ---
  real v_win = dec == 0 ? v0 : v1;
  real s_win = dec == 0 ? s0 : s1;
  real v_loss = dec == 0 ? v1 : v0;
  real s_loss = dec == 0 ? s1 : s0;
  
  // --- 4. pnegative Correction ---
  // Instead of computing 1 - normal_cdf(...), use the log-complement functions.
  // For v > 0, 1 - Phi(v) = exp(std_normal_lccdf(v)), and for v <= 0, Phi(-v) = exp(std_normal_lcdf(-v)).
  real pneg0 = v0 > 0 ? exp(std_normal_lccdf(v0 / s0)) : exp(std_normal_lcdf(-v0 / s0));
  real pneg1 = v1 > 0 ? exp(std_normal_lccdf(v1 / s1)) : exp(std_normal_lcdf(-v1 / s1));
  real pnegative = pneg0 * pneg1;
  real log_correction = -log(fmax(1 - pnegative, eps));
  
  // --- 5. Winner PDF in log-space ---
  real st_win = s_win * t_adj;
  st_win = fmax(st_win, eps);
  real z1_win = (b - A - v_win * t_adj) / st_win;
  real z2_win = (b - v_win * t_adj) / st_win;
  
  // Use the std_normal functions to compute differences in CDFs and PDFs.
  real diff_cdf = std_normal_cdf(z2_win) - std_normal_cdf(z1_win);
  real diff_pdf = exp(std_normal_lpdf(z1_win)) - exp(std_normal_lpdf(z2_win));
  real winner_dens = (v_win * diff_cdf + s_win * diff_pdf) / A;
  winner_dens = fmax(winner_dens, eps);
  real log_winner = log(winner_dens);
  
  // --- 6. Loser Survival Function ---
  real st_loss = s_loss * t_adj;
  st_loss = fmax(st_loss, eps);
  real z1_loss = (b - A - v_loss * t_adj) / st_loss;
  real z2_loss = (b - v_loss * t_adj) / st_loss;
  // Compute term for loser CDF using std_normal functions.
  real term_loss = (b - A - v_loss * t_adj) * std_normal_cdf(z1_loss)
                    + st_loss * exp(std_normal_lpdf(z1_loss))
                    - (b - v_loss * t_adj) * std_normal_cdf(z2_loss)
                    - st_loss * exp(std_normal_lpdf(z2_loss));
  real cdf_loss = 1 + (1 / A) * term_loss;
  cdf_loss = fmin(fmax(cdf_loss, 0), 1 - eps);
  real log_surv = log1m(cdf_loss);
  
  // --- 7. Combine Components ---
  real log_lik = log_winner + log_surv + log_correction;
  return (is_nan(log_lik) || is_inf(log_lik)) ? negative_infinity() : log_lik;
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

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
"functions {
", .lba_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$lba_lpdf
}



#' @rdname lba_brms
#' @export
lba_stanvars <- function() {
  brms::stanvar(scode = .lba_lpdf(), block = "functions")
}



#' @rdname lba_brms
#' @export
lba <- function(link_mu = "identity", link_vdelta = "identity",
                link_sigmazero = "softplus", link_sigmadelta = "identity",
                link_A = "softplus", link_k = "softplus",
                link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "lba",
    dpars = c("mu", "vdelta", "sigmazero", "sigmadelta", "A", "k", "tau", "minrt"), # Changed ndt to tau/minrt
    links = c(link_mu, link_vdelta, link_sigmazero, link_sigmadelta, link_A, link_k, link_tau, link_minrt), # Changed ndt to tau/minrt
    lb = c(NA, NA, 0, NA, 0, 0, 0, 0), # Lower bounds: tau>=0, minrt>=0
    ub = c(NA, NA, NA, NA, NA, NA, 1, NA), # Upper bound: tau<=1
    type = "real", # Response variable type
    vars = c("dec[n]") # Additional data variables needed
  )
}

# brms Post-processing Functions ------------------------------------------

#' @rdname lba_brms
#' @export
log_lik_lba <- function(i, prep) {
  # Extract observation index i
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_) # Handle missing RTs

  # Get parameters for observation i across all posterior draws
  vzero <- brms::get_dpar(prep, "mu", i = i)
  vdelta <- brms::get_dpar(prep, "vdelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  A <- brms::get_dpar(prep, "A", i = i)
  k <- brms::get_dpar(prep, "k", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i) # Get tau
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate ndt
  ndt <- tau * minrt

  # Get decision for observation i (should be constant across draws)
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning("Response (dec) must be 0 or 1. Got: ", response, " for observation ", i)
    return(rep(-Inf, length(vzero))) # Return -Inf for all draws if response is invalid
  }

  # Compute log-likelihood using the R density function
  # dlba should be vectorized over its parameters, expects ndt
  ll <- dlba(x = y, response = response,
             vzero = vzero, vdelta = vdelta,
             sigmazero = sigmazero, sigmadelta = sigmadelta,
             A = A, k = k, ndt = ndt, log = TRUE) # Pass calculated ndt

  ll # Return vector of log-likelihoods (one per draw)
}

#' @rdname lba_brms
#' @export
posterior_predict_lba <- function(i, prep, ...) {
  # Get parameters for observation i across all posterior draws
  vzero <- brms::get_dpar(prep, "mu", i = i)
  vdelta <- brms::get_dpar(prep, "vdelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  A <- brms::get_dpar(prep, "A", i = i)
  k <- brms::get_dpar(prep, "k", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i) # Get tau
  # Get minrt from the data (constant for this observation across draws)
  minrt <- prep$data[["minrt"]][i]
  if (is.na(minrt)) {
      warning("minrt data is NA for observation ", i, ". Cannot generate predictions.")
      return(rep(NA_real_, length(vzero)))
  }

  # Calculate ndt
  ndt <- tau * minrt

  # Number of posterior draws
  n_draws <- length(vzero) # Length of any parameter vector

  # Generate predictions using the R simulation function
  # rlba should be vectorized over its parameters, expects ndt
  sim_data <- rlba(n = n_draws, vzero = vzero, vdelta = vdelta,
                   sigmazero = sigmazero, sigmadelta = sigmadelta,
                   A = A, k = k, ndt = ndt) # Pass calculated ndt

  # Return the simulated reaction times (vector, one per draw)
  sim_data$rt
}

#' @rdname lba_brms
#' @export
posterior_epred_lba <- function(prep) {
  # Calculating the expected value (mean RT) for the LBA model is complex,
  # often requiring numerical integration or extensive simulation for each draw.
  # This is generally too computationally intensive for standard brms post-processing.
  # It's recommended to use posterior_predict_lba to get draws from the
  # predictive distribution and calculate expectations or other summaries manually.
  stop(
    "Calculating the posterior expected prediction (epred) for the LBA model ",
    "is computationally prohibitive within this framework.\n",
    "Please use `posterior_predict()` to obtain draws from the posterior ",
    "predictive distribution and calculate summaries manually if needed."
  )
}
