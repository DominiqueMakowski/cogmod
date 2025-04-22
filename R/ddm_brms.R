#' @keywords internal
.ddm_lpdf <- function() {
"
// Reparameterized Wiener diffusion log-PDF for a single response
// mu: drift rate parameter
// threshold: boundary separation parameter > 0
// bias: initial bias parameter in [0, 1]
// tau: non-decision time proportion in [0, 1]
// minrt: minimum reaction time > 0
real ddm_lpdf(real y, real mu, real bs, real bias,
              real tau, real minrt, int dec) {

    // derive non-decision time
    real ndt = tau * minrt;

    // Compute log-likelihood using Stan's wiener_lpdf
    if (dec == 1) { // Upper boundary response
        return wiener_lpdf(y | bs, ndt, bias, mu);
    } else { // Lower boundary response (dec == 0)
        // For lower boundary, flip bias (bias -> 1-bias) and drift (mu -> -mu)
        // when calling the function that calculates the density for the *upper* boundary.
        return wiener_lpdf(y | bs, ndt, 1 - bias, -mu);
    }
}
"
}



#' @rdname rddm
#' @export
ddm_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
"functions {
", .ddm_lpdf(), "
}")

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$ddm_lpdf
}


#' @rdname rddm
#' @export
ddm_stanvars <- function() {
  brms::stanvar(scode = .ddm_lpdf(), block = "functions")
}


#' @rdname rddm
#' @param link_mu,link_bs,link_bias,link_tau,link_minrt Link functions for the parameters.
#' @export
ddm <- function(link_mu = "identity", link_bs = "softplus", link_bias = "logit",
                link_tau = "logit", link_minrt = "identity") {
  brms::custom_family(
    name = "ddm",
    dpars = c("mu", "bs", "bias", "tau", "minrt"),
    links = c(link_mu, link_bs, link_bias, link_tau, link_minrt),
    lb = c(NA, 0, 0, 0, 0),  # Lower bounds for bs > 0, tau >= 0, minrt >= 0, bias >= 0
    ub = c(NA, NA, 1, 1, NA), # Upper bounds for bias <= 1, tau <= 1
    type = "real",            # Response variable type
    vars = "dec[n]"           # Required additional data variable for the decision
  )
}


# brms --------------------------------------------------------------------

#' @rdname rddm
#' @inheritParams rlnr
#' @export
log_lik_ddm <- function(i, prep, ...) { 
  # Extract parameters for observation i across all posterior draws (vectors)
  drift <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Vector if varying, scalar if fixed

  # Extract observed data for observation i (scalars)
  y <- prep$data$Y[i]      # Observed RT
  dec <- prep$data$dec[i]  # Observed decision

  # Handle missing data for the observation
  if (is.na(y) || is.na(dec)) {
    # Return vector of NAs, one for each draw
    return(rep(NA_real_, prep$ndraws))
  }

  # Basic check for valid response coding (0 or 1)
  if (!dec %in% c(0, 1)) {
    warning("Response ('dec') must be 0 or 1. Found: ", dec, " for observation ", i, ". Returning -Inf log-likelihood.")
    # Return -Inf for all draws for this observation
    return(rep(-Inf, prep$ndraws))
  }

  # Calculate ndt (vector, one element per draw)
  # Element-wise multiplication handles if minrt is scalar or vector
  ndt <- tau * minrt

  # Compute log-likelihood using brms::dwiener
  # dwiener is vectorized over parameters when x and resp are scalar
  ll <- dddm(x = y, drift = drift, bs = bs, bias = bias, ndt = ndt, response = dec, log = TRUE, ...)

  # Return vector of log-likelihoods (one per posterior draw)
  ll 
}



#' @rdname rddm
#' @inheritParams rlnr
#' @export
posterior_epred_ddm <- function(prep) {
  # Extract parameters
  mu <- prep$dpars$mu
  bs <- prep$dpars$bs
  bias <- prep$dpars$bias
  tau <- prep$dpars$tau
  minrt <- prep$dpars$minrt

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Compute expected reaction time (E[RT])
  # Formula adapted from https://doi.org/10.1016/j.jmp.2009.01.006
  epred <- ndt - bias / mu + bs / mu * 
    (exp(-2 * mu * bias) - 1) / (exp(-2 * mu * bs) - 1)

  # Return expected reaction time
  epred
}


#' @rdname rddm
#' @inheritParams rlnr
#' @export
posterior_predict_ddm <- function(i, prep, ...) {
  # Extract parameters for observation i
  mu <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Simulate reaction times and responses using rddm()
  sim_data <- rddm(n = prep$ndraws, drift = mu, bs = bs, bias = bias, ndt = ndt)

  # Return simulated reaction times
  as.matrix(sim_data)
}