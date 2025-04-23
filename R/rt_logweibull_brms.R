#' Shifted Log-Weibull Model (Gumbel Distribution)
#'
#' Provides the necessary functions to use a shifted Log-Weibull distribution
#' as a custom family in `brms`. This version is reparameterized to use
#' `tau` (proportion of non-decision time relative to minimum RT) and `minrt`
#' (minimum possible RT), where the non-decision time `ndt = tau * minrt`.
#' The distribution is parameterized by the shape (`mu`), scale (`sigma`),
#' as well as `tau`, and `minrt` (which is typically fixed).
#'
#' @param link_mu,link_sigma,link_tau,link_minrt Link functions for the parameters.
#'
#' @export
rt_logweibull <- function(link_mu = "softplus", link_sigma = "softplus",
                          link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "rt_logweibull", 
    dpars = c("mu", "sigma", "tau", "minrt"), 
    links = c(link_mu, link_sigma, link_tau, link_minrt),
    lb = c(0, 0, 0, 0), # Lower bounds: mu > 0, sigma > 0, tau >= 0, minrt > 0
    ub = c(NA, NA, 1, NA), # Upper bound: tau <= 1
    type = "real" # Continuous outcome variable (RT)
  )
}

#' @keywords internal
.rt_logweibull_lpdf <- function() {
"
// Log-likelihood for a single observation from the Shifted Log-Weibull distribution.
// Y: observed reaction time.
// mu: shape parameter of the Log-Weibull distribution (> 0).
// sigma: scale parameter of the Log-Weibull distribution (> 0).
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
// minrt: Minimum possible reaction time (> 0).
real rt_logweibull_lpdf(real Y, real mu, real sigma, real tau, real minrt) {
    // Parameter checks
    if (mu <= 0 || sigma <= 0 || tau < 0 || tau > 1 || minrt <= 0) return negative_infinity();

    // Compute non-decision time and adjusted time
    real ndt   = tau * minrt;
    real t_adj = Y - ndt;
    if (t_adj <= 0) return negative_infinity(); // Density is 0 if Y <= ndt

    // Use Stan's built-in Gumbel lpdf
    return gumbel_lpdf(log(t_adj) | mu, sigma) - log(t_adj);
}
"
}

#' @rdname rt_logweibull
#' @export
rt_logweibull_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block
  stancode <- paste0(
"functions {",
.rt_logweibull_lpdf(),
"}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$rt_logweibull_lpdf
}

#' @rdname rt_logweibull
#' @export
rt_logweibull_stanvars <- function() {
  brms::stanvar(scode = .rt_logweibull_lpdf(), block = "functions")
}

#' @rdname rt_logweibull
#' @export
log_lik_rt_logweibull <- function(i, prep) {
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)

  mu    <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0))

  y_vec <- rep(y, length.out = n_draws)
  ndt <- tau * minrt
  t_adj <- y_vec - ndt

  ll <- numeric(length(t_adj))
  valid_idx <- t_adj > 0

  if (any(valid_idx)) {
    log_t_adj <- log(t_adj[valid_idx])
    z <- (log_t_adj - mu[valid_idx]) / sigma[valid_idx]
    ll[valid_idx] <- -log(sigma[valid_idx]) - z - exp(-z) - log_t_adj
  }

  ll[!valid_idx] <- -Inf
  ll[is.nan(ll) | is.na(ll)] <- -Inf

  ll
}



#' @rdname rt_logweibull
#' @inheritParams rt_invgaussian
#' @export
posterior_predict_rt_logweibull <- function(i, prep, ...) {
  # Get parameters for observation i across all draws
  mu    <- brms::get_dpar(prep, "mu", i = i)    # Gumbel location
  sigma <- brms::get_dpar(prep, "sigma", i = i) # Gumbel scale
  tau   <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  n_draws <- length(mu)

  # Check parameter validity
  if (any(sigma <= 0)) stop("Scale parameter must be greater than 0.")

  # Calculate non-decision time
  ndt <- tau * minrt

  # Generate Gumbel samples for log(t_adj)
  u <- stats::runif(n_draws)
  log_t_adj <- mu - sigma * log(-log(u))
  t_adj <- exp(log_t_adj)

  # Add non-decision time
  final_out <- t_adj + ndt

  as.matrix(final_out)
}

#' @rdname rt_logweibull
#' @export
posterior_epred_rt_logweibull <- function(prep) {
  mu    <- brms::get_dpar(prep, "mu")
  sigma <- brms::get_dpar(prep, "sigma")
  tau   <- brms::get_dpar(prep, "tau")
  minrt <- brms::get_dpar(prep, "minrt")

  ndt <- tau * minrt
  epred <- exp(mu + sigma * 0.5772156649) + ndt

  epred
}