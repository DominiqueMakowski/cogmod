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
#' @examples
#' # Simulate data
#' data <- rlnr(1000, mu = 0, mudelta = 0.5, sigmazero = 1, sigmadelta = -0.5, ndt = 0.2)
#' hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' @export
rlnr <- function(n, mu = 1, mudelta = 0, sigmazero = 1, sigmadelta = 0, ndt = 0.2) {
  # Compute the means and standard deviations for both accumulators
  nu <- c(mu, mu + mudelta)
  sigma <- c(sigmazero, sigmazero * exp(sigmadelta))

  # Generate log-normal draws for both accumulators across all trials
  draws <- matrix(stats::rlnorm(2 * n, meanlog = rep(nu, each = n), sdlog = rep(sigma, each = n)), nrow = n, ncol = 2) + ndt

  # Determine responses and reaction times
  response <- apply(draws, 1, which.min) - 1  # 0-based index
  rt <- draws[cbind(seq_len(n), response + 1)]

  data.frame(rt = rt, response = response)
}



# dlnr: computes the log-density for one observation from the LNR model.
#' @rdname rlnr
#' @param x The observed reaction time (RT). Must be greater than `ndt`.
#' @param response The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.
#' @param log Logical; if TRUE, returns the log-density. Default: FALSE.
dlnr <- function(x, mu, mudelta, sigmazero, sigmadelta, ndt, response, log = FALSE) {
  eps <- 1e-6
  
  # Input validation
  if (!all(response %in% c(0, 1))) {
    stop("response must be 0 or 1")
  }
  
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
  
  # Compute adjusted reaction times
  t_adj <- x - ndt
  
  # Initialize log-density vector
  log_density <- numeric(n)
  
  # Handle each entry
  for (i in 1:n) {
    # Short-circuit for invalid RTs
    if (t_adj[i] < eps) {
      log_density[i] <- -Inf
      next
    }
    
    # Precompute accumulator parameters
    nu0 <- mu[i]
    nu1 <- mu[i] + mudelta[i]
    sigma0 <- sigmazero[i]
    sigma1 <- sigmazero[i] * exp(sigmadelta[i])
    
    # Compute log-density for the winning accumulator
    log_pdf_win <- if (response[i] == 0) {
      stats::dlnorm(t_adj[i], meanlog = nu0, sdlog = sigma0, log = TRUE)
    } else {
      stats::dlnorm(t_adj[i], meanlog = nu1, sdlog = sigma1, log = TRUE)
    }
    
    # Compute log-survival probability for the losing accumulator
    log_cdf_loss <- if (response[i] == 0) {
      stats::plnorm(t_adj[i], meanlog = nu1, sdlog = sigma1, lower.tail = TRUE, log.p = TRUE)
    } else {
      stats::plnorm(t_adj[i], meanlog = nu0, sdlog = sigma0, lower.tail = TRUE, log.p = TRUE)
    }
    
    # Use stable computation for log(1-p)
    log_surv_loss <- .log1m_exp(log_cdf_loss)
    
    # Combine log-density and log-survival probability
    log_density[i] <- log_pdf_win + log_surv_loss
  }
  
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}