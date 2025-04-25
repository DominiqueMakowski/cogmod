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
  real A    = sigmabias;      // start uniform [0,A]
  real b    = A + bs;         // threshold
  real ndt  = tau * minrt;    // non-decision time
  real eps  = 1e-10;

  // Must exceed non-dec time
  if (Y <= ndt) return negative_infinity();
  real T = Y - ndt;
  if (T <= 0) return negative_infinity();

  // st = sigma * T
  real st = sigma * T;
  st = fmax(st, eps);

  // Standardize
  real z = mu / sigma;
  real z1    = (b - A - mu * T) / st; // = (b - A)/(sigma * T) - z
  real z2    = (b - mu * T)   / st;   // = b/(sigma * T) - z

  // Phi(z)
  real phi1 = exp(std_normal_lpdf(z1));
  real phi2 = exp(std_normal_lpdf(z2));

  // Phi(z) - In Stan's log-CDF world (un-log at the end)
  real lcdf1 = std_normal_lcdf(z1);  // log phi(z1)
  real lcdf2 = std_normal_lcdf(z2);  // log phi(z2)
  real cdf1  = exp(lcdf1);
  real cdf2  = exp(lcdf2);

  // Defective density:
  // comp = mu*(phi(z2)-phi(z1)) + sigma*(phi(z1)-phi(z2))
  real comp = mu * (cdf2 - cdf1)
            + sigma * (phi1  - phi2);
  comp = fmax(comp, eps);

  // log normalization = log[1 - phi(-mu/sigma)] = log[1 - phi(-z)]
  // = log1m_exp(log phi(-z)) because log phi(-z) = std_normal_lcdf(-z)
  real log_denom = log1m_exp(std_normal_lcdf(-z));

  // Final log density
  // = -log(A) + log(comp) - log(1 - phi(-z))
  return -log(A) + log(comp) - log_denom;
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
rt_lba <- function(link_mu = "softplus",
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
#' @export
log_lik_rt_lba <- function(i, prep) {
  # Extract the observed reaction time for observation i
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)  # Handle missing RTs

  # Extract posterior draws for observation i.
  # (Note: Under the new single-accumulator parametrization, we have one drift parameter.)
  drift     <- brms::get_dpar(prep, "mu", i = i)      # drift parameter
  sigma  <- brms::get_dpar(prep, "sigma", i = i)   # drift rate SD
  sigmabias  <- brms::get_dpar(prep, "sigmabias", i = i)
  bs         <- brms::get_dpar(prep, "bs", i = i)
  tau        <- brms::get_dpar(prep, "tau", i = i)
  minrt      <- brms::get_dpar(prep, "minrt", i = i)     # Minimum possible RT

  # Calculate non-decision time from tau and minrt.
  ndt <- tau * minrt

  # For a single-accumulator LBA, there is no response indicator,
  # so we simply compute the log density via our R density function `dlba`
  # (which here expects: x (RT), drift, sigma, sigmabias, bs, and ndt).
  ll <- drt_lba(x = y,
             drift = drift,
             sigma = sigma,
             sigmabias = sigmabias,
             bs = bs,
             ndt = ndt,
             log = TRUE)
  ll  # Return the vector of log-likelihoods (one per posterior draw)
}

#' @rdname rrt_lba
#' @inheritParams rlnr
#' @export
posterior_predict_rt_lba <- function(i, prep, ...) {
  # Extract the posterior draws (using the new parametrization)
  drift    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  bs        <- brms::get_dpar(prep, "bs", i = i)
  tau       <- brms::get_dpar(prep, "tau", i = i)
  minrt     <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time.
  ndt <- tau * minrt

  # Number of posterior draws.
  n_draws <- length(drift)

  # Generate predictions via the simulation function `rrt_lba`
  # (which now is for a single accumulator with parameters: drift, sigma, sigmabias, bs, ndt)
  sim_data <- rrt_lba(n = n_draws,
                      drift   = drift,
                      sigma   = sigma,
                      sigmabias = sigmabias,
                      bs      = bs,
                      ndt     = ndt)
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
