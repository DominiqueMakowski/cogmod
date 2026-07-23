# Stanvars ----------------------------------------------------------------

#' @keywords internal
.rt_lba_lpdf <- function() {
  "
// Log density function for the one-accumulator LBA model.
real rt_lba_lpdf(
  real Y,
  real mu,         // drift rate (named 'mu' to match Stan's convention)
  real sigma,      // standard deviation of drift
  real sigmabias,  // A = sigmabias
  real bs,         // so b = A + bs
  real tau,        // scales minrt to yield ndt = tau*minrt
  real minrt       // minimum possible RT
) {
  // Derived params
  real A   = sigmabias;      // start uniform [0,A]
  real ndt = tau * minrt;    // non-decision time

  // Must exceed non-dec time, and A/sigma must be valid
  if (Y <= ndt) return negative_infinity();
  if (sigma <= 0 || A <= 0) return negative_infinity();

  real t = Y - ndt;

  real inv_A = inv(A);
  real INV_SQRT_2PI = 0.3989422804014327;

  // Standardize (st = sigma * t, clamped away from 0)
  real st = fmax(sigma * t, 1e-10);
  real inv_st = inv(st);

  real z  = mu / sigma;
  real z1 = (bs - mu * t) * inv_st;  // = (b - A - mu*t) / st
  real z2 = z1 + A * inv_st;         // = (b - mu*t) / st

  // Standard normal CDFs
  real Phi1 = Phi(z1);
  real Phi2 = Phi(z2);

  // Standard normal PDFs
  real phi1 = exp(-0.5 * square(z1)) * INV_SQRT_2PI;
  real phi2 = exp(-0.5 * square(z2)) * INV_SQRT_2PI;

  // Defective density:
  // f = (1/A) * [mu*(Phi(z2)-Phi(z1)) + sigma*(phi(z1)-phi(z2))]
  real f = inv_A * (mu * (Phi2 - Phi1) + sigma * (phi1 - phi2));

  // log normalization = log[1 - Phi(-mu/sigma)] = log[1 - Phi(-z)]
  // = log1m_exp(log Phi(-z)) because log Phi(-z) = std_normal_lcdf(-z)
  real log_denom = log1m_exp(std_normal_lcdf(-z));

  // Final log density
  return log(fmax(f, 1e-10)) - log_denom;
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
",
    .rt_lba_lpdf(),
    "
}"
  )

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
rt_lba <- function(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_sigmabias = "softplus",
  link_bs = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
) {
  brms::custom_family(
    name = "rt_lba",
    dpars = c("mu", "sigma", "sigmabias", "bs", "tau", "minrt"),
    links = c(
      link_mu,
      link_sigma,
      link_sigmabias,
      link_bs,
      link_tau,
      link_minrt
    ),
    lb = c(NA, 0, 0, 0, 0, 0),
    ub = c(NA, NA, NA, NA, 1, NA),
    type = "real"
  )
}

# brms Post-processing Functions ------------------------------------------

#' @rdname rrt_lba
#' @export
log_lik_rt_lba <- function(i, prep) {
  # Extract the observed reaction time for observation i
  y <- prep$data$Y[i]
  if (is.na(y)) {
    return(NA_real_)
  } # Handle missing RTs

  # Extract posterior draws for observation i.
  # (Note: Under the new single-accumulator parametrization, we have one drift parameter.)
  drift <- brms::get_dpar(prep, "mu", i = i) # drift parameter
  sigma <- brms::get_dpar(prep, "sigma", i = i) # drift rate SD
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Minimum possible RT

  # Calculate non-decision time from tau and minrt.
  ndt <- tau * minrt

  # For a single-accumulator LBA, there is no response indicator,
  # so we simply compute the log density via our R density function `dlba`
  # (which here expects: x (RT), drift, sigma, sigmabias, bs, and ndt).
  ll <- drt_lba(
    x = y,
    drift = drift,
    sigma = sigma,
    sigmabias = sigmabias,
    bs = bs,
    ndt = ndt,
    log = TRUE
  )
  ll # Return the vector of log-likelihoods (one per posterior draw)
}

#' @rdname rrt_lba
#' @inheritParams rlnr
#' @export
posterior_predict_rt_lba <- function(i, prep, ...) {
  # Extract the posterior draws (using the new parametrization)
  drift <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time.
  ndt <- tau * minrt

  # Number of posterior draws.
  n_draws <- length(drift)

  # Generate predictions via the simulation function `rrt_lba`
  # (which now is for a single accumulator with parameters: drift, sigma, sigmabias, bs, ndt)
  sim_data <- rrt_lba(
    n = n_draws,
    drift = drift,
    sigma = sigma,
    sigmabias = sigmabias,
    bs = bs,
    ndt = ndt
  )
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
