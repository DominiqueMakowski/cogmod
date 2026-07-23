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
  real ndt = tau * minrt;         // Non-decision time.

  // Check if the observed RT is less than or equal to ndt.
  if (Y <= ndt) return negative_infinity();

  // Now compute the adjusted decision time.
  real t = Y - ndt;

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

  real inv_A = inv(A);
  real INV_SQRT_2PI = 0.3989422804014327;

  // --- Winning accumulator: defective density ---
  real st_win = fmax(s_win * t, 1e-10);
  real inv_st_win = inv(st_win);
  real z1w = (bs - v_win * t) * inv_st_win;        // = (b - A - v_win*t) / st_win
  real z2w = z1w + A * inv_st_win;                 // = (b - v_win*t) / st_win
  real Phi1w = Phi(z1w);
  real Phi2w = Phi(z2w);
  real phi1w = exp(-0.5 * square(z1w)) * INV_SQRT_2PI;
  real phi2w = exp(-0.5 * square(z2w)) * INV_SQRT_2PI;
  real f_win = inv_A * (v_win * (Phi2w - Phi1w) + s_win * (phi1w - phi2w));

  // --- Losing accumulator: survival probability ---
  real st_loss = fmax(s_loss * t, 1e-10);
  real inv_st_loss = inv(st_loss);
  real z1l = (bs - v_loss * t) * inv_st_loss;      // = (b - A - v_loss*t) / st_loss
  real z2l = z1l + A * inv_st_loss;                // = (b - v_loss*t) / st_loss
  real Phi1l = Phi(z1l);
  real Phi2l = Phi(z2l);
  real phi1l = exp(-0.5 * square(z1l)) * INV_SQRT_2PI;
  real phi2l = exp(-0.5 * square(z2l)) * INV_SQRT_2PI;
  real g1 = z1l * Phi1l + phi1l;
  real g2 = z2l * Phi2l + phi2l;
  real S_loss = (st_loss * inv_A) * (g2 - g1);     // = 1 - CDF_loss(t)

  // --- Combine into the joint (defective) density ---
  // Floor the *joint* density (not each factor separately) at double
  // precision epsilon, mirroring the R reference implementation `dlba()`.
  // Flooring each factor individually (e.g. at 1e-10) is much more
  // aggressive and can bias the log-density in the tails relative to R.
  real EPS = 2.220446049250313e-16;  // .Machine$double.eps
  real dens = fmax(f_win, 0) * fmin(fmax(S_loss, 0), 1);
  real log_lik = log(fmax(dens, EPS));

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
",
    .lba_lpdf(),
    "
}"
  )

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
lba <- function(
  link_mu = "identity",
  link_driftone = "identity",
  link_sigmazero = "softplus",
  link_sigmaone = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
) {
  brms::custom_family(
    name = "lba",
    dpars = c(
      "mu",
      "driftone",
      "sigmazero",
      "sigmaone",
      "sigmabias",
      "bs",
      "tau",
      "minrt"
    ),
    links = c(
      link_mu,
      link_driftone,
      link_sigmazero,
      link_sigmaone,
      link_sigmabias,
      link_bs,
      link_tau,
      link_minrt
    ),
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
  if (is.na(y)) {
    return(NA_real_)
  } # Handle missing RTs

  # Get the parameter draws for observation i:
  driftzero <- brms::get_dpar(prep, "mu", i = i)
  driftone <- brms::get_dpar(prep, "driftone", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmaone <- brms::get_dpar(prep, "sigmaone", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Minimum possible reaction time

  # Calculate non-decision time from tau and minrt
  ndt <- tau * minrt

  # Get decision indicator for observation i (should be 0 or 1).
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning(
      "Response (dec) must be 0 or 1. Got: ",
      response,
      " for observation ",
      i
    )
    return(rep(-Inf, length(driftzero)))
  }

  # Compute log-likelihood using our R density function 'dlba'
  ll <- dlba(
    x = y,
    response = response,
    driftzero = driftzero,
    driftone = driftone,
    sigmazero = sigmazero,
    sigmaone = sigmaone,
    sigmabias = sigmabias,
    bs = bs,
    ndt = ndt,
    log = TRUE
  )
  ll # Return vector of log-likelihoods (one per draw)
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
  sim_data <- rlba(
    n = n_draws,
    driftzero = driftzero,
    driftone = driftone,
    sigmazero = sigmazero,
    sigmaone = sigmaone,
    sigmabias = sigmabias,
    bs = bs,
    ndt = ndt
  )
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
