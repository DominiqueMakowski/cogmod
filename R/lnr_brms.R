# Stanvars ----------------------------------------------------------------

#' @keywords internal
.lnr_lpdf <- function() {
"
// Log-likelihood for a single observation from the reparameterized Log-Normal Race model.
// Y: observed reaction time.
// dec: decision indicator (0 or 1).
// mu: baseline accumulator mean (in log-space) for choice 0.
// mudelta: log-deviation for the mean of choice 1.
// sigmazero: baseline accumulator standard deviation (in log-space) for choice 0.
// sigmadelta: log-deviation for the standard deviation of choice 1.
// tau: Scale factor for non-decision time (0-1, scaled by minimum RT).
real lnr_lpdf(real Y, real mu, real mudelta, real sigmazero, real sigmadelta, real tau, real minrt, int dec) {
  real eps = 1e-6;

  // Ensure positive standard deviations with minimal clipping
  real sig0 = fmax(sigmazero, eps);
  real sig1 = fmax(sig0 * exp(sigmadelta), eps);

  // Calculate non-decision time with extra safeguards
  real ndt = tau * minrt;

  // More robust RT adjustment - ensure it's positive
  real ndt_adj = Y - ndt;
  if (ndt_adj < eps)
    return negative_infinity();

  // Convert 0-based decision to 1-based indexing for Stan
  int dec_stan = dec + 1;

  // Define vectors in a more explicit, safer way
  vector[2] nu;
  vector[2] sigma;
  nu[1] = mu;
  nu[2] = mu + mudelta;
  sigma[1] = sig0;
  sigma[2] = sig1;

  real lp = 0;
  for (i in 1:2) {
    if (i == dec_stan) {
      // The chosen accumulator: add log density with bounds check
      lp += lognormal_lpdf(ndt_adj | nu[i], sigma[i]);
    } else {
      // Non-chosen accumulator: add log survival function with better numerical handling
      real log_cdf = lognormal_lcdf(ndt_adj | nu[i], sigma[i]);

      // More forgiving check for extreme CDF values
      if (log_cdf > -1e-10)  // Near-1 CDF is problematic
        return negative_infinity();

      // More stable survival calculation
      if (log_cdf < -20)  // If survival is basically 1
        lp += 0;  // log(1) = 0
      else
        lp += log1m_exp(log_cdf);
    }
  }
  return lp;
}
"
}

#' @rdname rlnr
#' @examples
#' # You can expose the lpdf function as follows:
#' # lnr_lpdf <- lnr_lpdf_expose()
#' # lnr_lpdf(Y = 0.5, mu = 0, mudelta = 0, sigmazero = 1, sigmadelta = 0.5,
#' #          tau = 0.1, minrt = 0.2, dec = 1)
#'
#' @export
lnr_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
"functions {
", .lnr_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$lnr_lpdf
}




#' @rdname rlnr
#' @export
lnr_stanvars <- function() {
  brms::stanvar(scode = .lnr_lpdf(), block = "functions")
}


#' @rdname rlnr
#' @param link_mu Link function for the `mu` parameter in the custom family.
#'   Determines how `mu` is transformed in the model. Default: "identity".
#' @param link_mudelta Link function for the `mudelta` parameter in the custom family.
#'   Determines how `mudelta` is transformed in the model. Default: "identity".
#' @param link_sigmazero Link function for the `sigmazero` parameter in the custom family.
#'   Ensures `sigmazero` remains positive. Default: "identity".
#' @param link_sigmadelta Link function for the `sigmadelta` parameter in the custom family.
#'   Determines how `sigmadelta` is transformed in the model. Default: "identity".
#' @param link_tau Link function for the `tau` parameter in the custom family.
#' @param link_minrt Link function for the `minrt` parameter in the custom family.
#' @export
lnr <- function(link_mu = "identity", link_mudelta = "identity",
                link_sigmazero = "softplus", link_sigmadelta = "identity",
                link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "lnr",
    dpars = c("mu", "mudelta", "sigmazero", "sigmadelta", "tau", "minrt"),
    links = c(link_mu, link_mudelta, link_sigmazero, link_sigmadelta, link_tau, link_minrt),
    lb = c(NA, NA, 0, NA, 0, 0),
    ub = c(NA, NA, NA, NA, 1, NA),    # Ensure tau stays <= 1
    vars = "dec[n]"
  )
}



# # brms --------------------------------------------------------------------

#' @rdname rlnr
#' @inheritParams choco
#' @export
posterior_predict_lnr <- function(i, prep, ...) {
  # Get parameters
  mu <- brms::get_dpar(prep, "mu", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)  # Get minrt parameter

  # Calculate non-decision time
  ndt <- tau * minrt

  # Generate predictions
  n_draws <- length(mu)
  sim_data <- rlnr(n_draws, mu = mu, mudelta = mudelta,
                  sigmazero = sigmazero, sigmadelta = sigmadelta,
                  ndt = ndt)

  as.matrix(sim_data)
}

#' @rdname rlnr
log_lik_lnr <- function(i, prep) {
  # Extract observation
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA)

  # Get parameters
  mu <- brms::get_dpar(prep, "mu", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  sigmazero <- brms::get_dpar(prep, "sigmazero", i = i)
  sigmadelta <- brms::get_dpar(prep, "sigmadelta", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Calculate non-decision time
  ndt <- tau * minrt

  # Get decision
  response <- prep$data[["dec"]][i]
  if (!response %in% c(0, 1)) {
    warning("Response must be 0 or 1. Got: ", response)
    return(rep(-Inf, length(mu)))
  }

  # Compute log-likelihood
  ll <- dlnr(x = y, mu = mu, mudelta = mudelta,
           sigmazero = sigmazero, sigmadelta = sigmadelta,
           ndt = ndt, response = response, log = TRUE)

  ll
}
