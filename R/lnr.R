#' @title Log-Normal Race (LNR) Model
#'
#' @description
#' The Log-Normal Race (LNR) model is useful for modeling reaction times and choices in decision-making tasks.
#' The model assumes that each choice option (accumulator) draws a processing time from a LogNormal distribution.
#' The winning accumulator (minimum draw) determines the observed reaction time and choice.
#' The observed RT includes a non-decision time component (tau).
#'
#' Functions:
#' - `rlnr()`: Simulates random draws from the LNR model.
#' - `dlnr()`: Computes the likelihood/density of observed reaction times under the LNR model.
#' - `lnr()`: Creates a custom family to be used with `brms`.
#' - `lnr_stanvars()`: For `brms`, generates a `stanvars` object to pass to `brm()` when fitting the model.
#' - `posterior_predict_lnr()`: For `brms`, simulates predicted outcomes using sampled parameters.
#' - `log_lik_lnr()`: For `brms`, computes the log-likelihood of observed data.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param mu The log-space mean parameter for the baseline accumulator (choice 0).
#'   Controls the central tendency of the reaction time for choice 0.
#'   Can take any real value, with larger values leading to longer RTs. Range: (-Inf, Inf).
#' @param mudelta The additive deviation (in log-space) for the mean of accumulator 1 (choice 1).
#'   Positive values make choice 1 slower on average, while negative values make it faster.
#'   Can take any real value. Range: (-Inf, Inf).
#' @param sigmazero The log-space standard deviation for the baseline accumulator (choice 0).
#'   Controls the variability of reaction times for choice 0. Must be positive.
#'   Larger values increase variability. Range: (0, Inf).
#' @param sigmadelta The log-deviation for the standard deviation of accumulator 1.
#'   If positive, accumulator 1 has more variability; if negative, it has less variability
#'   compared to accumulator 0. Can take any real value. Range: (-Inf, Inf).
#' @param ndt Non-decision time (shift parameter). Represents the time taken for processes
#'   unrelated to the decision (e.g., encoding, motor response).
#'   Must be non-negative. Range: [0, Inf).
#' 
#' @details 
#' The LNR model conceptualizes decision-making as a race between multiple independent accumulators, 
#' each corresponding to a potential choice. Each accumulator gathers evidence at a fixed rate, 
#' and the first to reach a threshold determines the decision. In contrast, Wiener (DDM) models 
#' decision-making as a continuous stochastic process where a single decision variable accumulates 
#' noisy evidence over time until it reaches one of two decision boundaries.
#'
#' @examples
#' # Simulate data
#' data <- rlnr(1000, mu = 0, mudelta = 0.5, sigmazero = 1, sigmadelta = -0.5, ndt = 0.2)
#' hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' @export
rlnr <- function(n, mu = 1, mudelta = 0, sigmazero = 1, sigmadelta = 0, ndt = 0.2) {
  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) stop("n must be a positive integer.")
  if (any(sigmazero <= 0)) stop("sigmazero must be positive.")
  if (any(ndt < 0)) stop("ndt must be non-negative.")

  # --- Vectorization ---
  # Determine output length based on n and parameter vector lengths
  param_lengths <- c(length(mu), length(mudelta), length(sigmazero),
                     length(sigmadelta), length(ndt))
  n_params <- max(param_lengths)
  n_out <- max(n, n_params)

  # Recycle parameters to match output length
  if (n_out > 1) {
      mu <- rep(mu, length.out = n_out)
      mudelta <- rep(mudelta, length.out = n_out)
      sigmazero <- rep(sigmazero, length.out = n_out)
      sigmadelta <- rep(sigmadelta, length.out = n_out)
      ndt <- rep(ndt, length.out = n_out)
  }

  # --- Simulation ---
  # Compute the means and standard deviations for both accumulators (vectorized)
  nu0 <- mu
  nu1 <- mu + mudelta
  sigma0 <- sigmazero
  sigma1 <- sigmazero * exp(sigmadelta)

  # Generate log-normal draws for both accumulators (vectorized)
  # Each column corresponds to an accumulator, each row to a draw/observation
  draws0 <- stats::rlnorm(n_out, meanlog = nu0, sdlog = sigma0)
  draws1 <- stats::rlnorm(n_out, meanlog = nu1, sdlog = sigma1)

  # Determine responses and reaction times (vectorized)
  response <- ifelse(draws0 < draws1, 0, 1)
  rt_decision <- pmin(draws0, draws1) # Minimum of the two decision times
  rt <- rt_decision + ndt # Add non-decision time

  # Return data frame matching the determined output length
  data.frame(rt = rt, response = response)
}




#' @rdname rlnr
#' @param x The observed reaction time (RT). Must be greater than `ndt`.
#' @param response The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.
#' @param log Logical; if TRUE, returns the log-density. Default: FALSE.
#' @export
dlnr <- function(x, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = FALSE) {
  eps <- 1e-9 # Tolerance for checking RT > ndt

  # --- Input Validation ---
  if (any(sigmazero <= 0)) {
      warning("sigmazero must be positive. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(ndt < 0)) {
      warning("ndt must be non-negative. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  # For invalid response, just return 0/-Inf without warning (treat as impossible data)
  if (any(!response %in% c(0, 1))) {
      # stop("response must be 0 or 1") # Changed from stop
      return(ifelse(log, -Inf, 0))
  }


  # --- Vectorization ---
  # Recycle vectors to the same length
  n <- max(length(x), length(mu), length(mudelta), length(sigmazero),
           length(sigmadelta), length(ndt), length(response))

  x <- rep_len(x, n)
  mu <- rep_len(mu, n)
  mudelta <- rep_len(mudelta, n)
  sigmazero <- rep_len(sigmazero, n)
  sigmadelta <- rep_len(sigmadelta, n)
  ndt <- rep_len(ndt, n)
  response <- rep_len(response, n)

  # --- Calculation ---
  # Compute adjusted reaction times
  t_adj <- x - ndt

  # Initialize log-density vector with -Inf for invalid RTs (t_adj < eps)
  log_density <- ifelse(t_adj < eps, -Inf, 0) # Use 0 for now, will add log-densities

  # Identify valid RTs for calculation
  valid_idx <- which(t_adj >= eps)
  if (length(valid_idx) == 0) { # No valid RTs
      return(ifelse(log, log_density, exp(log_density)))
  }

  # Subset parameters and adjusted times for valid indices
  t_adj_valid <- t_adj[valid_idx]
  mu_valid <- mu[valid_idx]
  mudelta_valid <- mudelta[valid_idx]
  sigmazero_valid <- sigmazero[valid_idx]
  sigmadelta_valid <- sigmadelta[valid_idx]
  response_valid <- response[valid_idx]

  # Precompute accumulator parameters for valid indices
  nu0 <- mu_valid
  nu1 <- mu_valid + mudelta_valid
  sigma0 <- sigmazero_valid
  sigma1 <- sigmazero_valid * exp(sigmadelta_valid)

  # Identify winning and losing parameters based on response
  nu_win <- ifelse(response_valid == 0, nu0, nu1)
  sigma_win <- ifelse(response_valid == 0, sigma0, sigma1)
  nu_loss <- ifelse(response_valid == 0, nu1, nu0)
  sigma_loss <- ifelse(response_valid == 0, sigma1, sigma0)

  # Compute log-density for the winning accumulator (vectorized)
  log_pdf_win <- stats::dlnorm(t_adj_valid, meanlog = nu_win, sdlog = sigma_win, log = TRUE)

  # Compute log-survival probability for the losing accumulator (vectorized)
  log_cdf_loss <- stats::plnorm(t_adj_valid, meanlog = nu_loss, sdlog = sigma_loss, lower.tail = TRUE, log.p = TRUE)
  log_surv_loss <- .log1m_exp(log_cdf_loss) # Stable log(1 - p)

  # Combine log-density and log-survival probability for valid indices
  log_density[valid_idx] <- log_pdf_win + log_surv_loss

  # Handle potential -Inf results from dlnorm/plnorm if parameters were invalid
  # (although sigmazero check should prevent most)
  log_density[is.na(log_density) | !is.finite(log_density)] <- -Inf

  # Return result
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}