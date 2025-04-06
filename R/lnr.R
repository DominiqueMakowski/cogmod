#' @title Log-Normal Race (LNR) Model
#'
#' @description
#' The Log-Normal Race (LNR) model is useful for modeling reaction times and errors in decision-making tasks.
#' The model assumes that each accumulator draws a value from a LogNormal distribution (shifted by a non-decision time τ).
#' The winning accumulator (minimum draw) determines the observed reaction time and choice.
#'
#' Functions:
#' - `rlnr()`: Simulates random draws from the LNR model.
#' - `dlnr()`: Computes the log-likelihood of observed reaction times under the LNR model.
#' - `lnr()`: Creates a custom family to be used with `brms`.
#' - `lnr_stanvars()`: For `brms`, generates a `stanvars` object to pass to `brm()` when fitting the model.
#' - `posterior_predict_lnr()`: For `brms`, simulates predicted outcomes using sampled parameters.
#' - `log_lik_lnr()`: For `brms`, computes the log-likelihood of observed data.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param mu The log-space mean (ν) for the baseline accumulator (choice 0).
#'   Determines the central tendency of the reaction time for choice 0.
#'   Plausible range: (-Inf, Inf).
#' @param mudelta The additive deviation (in log-space) from `mu` to obtain the mean for accumulator 1 (choice 1).
#'   Positive values make choice 1 faster on average, while negative values make it slower.
#'   Plausible range: (-Inf, Inf).
#' @param sigmazero The log-space standard deviation for the baseline accumulator (choice 0).
#'   Controls the variability of reaction times for choice 0. Must be positive.
#'   Plausible range: (0, Inf).
#' @param sigmadelta The log-deviation for the standard deviation so that the standard deviation for accumulator 1 is `sigmazero * exp(sigmadelta)`.
#'   Positive values increase variability for choice 1, while negative values decrease it.
#'   Plausible range: (-Inf, Inf).
#' @param tau A non-decision time (shift), τ. Represents the time taken for processes unrelated to the decision (e.g., motor response).
#'   Must be non-negative. Plausible range: [0, Inf).
#'
#' @examples
#' # Simulate data
#' data <- rlnr(1000, mu = 0, mudelta = 0.5, sigmazero = 1, sigmadelta = -0.5, tau = 0.2)
#' hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' @export
rlnr <- function(n, mu = 1, mudelta = 0, sigmazero = 1, sigmadelta = 0, tau = 0.2) {
  # Compute the means and standard deviations for both accumulators
  nu <- c(mu, mu + mudelta)
  sigma <- c(sigmazero, sigmazero * exp(sigmadelta))

  # Generate log-normal draws for both accumulators across all trials
  draws <- matrix(stats::rlnorm(2 * n, meanlog = rep(nu, each = n), sdlog = rep(sigma, each = n)), nrow = n, ncol = 2) + tau

  # Determine responses and reaction times
  response <- apply(draws, 1, which.min) - 1  # 0-based index
  rt <- draws[cbind(seq_len(n), response + 1)]

  data.frame(rt = rt, response = response)
}



# dlnr: computes the log-density for one observation from the LNR model.
#' @rdname rlnr
#' @param y The observed reaction time (RT).
#' @param response The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.
#' @param log Logical; if TRUE, returns the log-density. Default: FALSE.
dlnr <- function(y, mu, mudelta, sigmazero, sigmadelta, tau, response, log = FALSE) {
  eps <- 1e-6
  if (!response %in% c(0, 1)) {
    stop("response must be 0 or 1")
  }
  # Compute the adjusted reaction time for each draw
  t_adj <- y - tau
  if (t_adj < eps) {
    return(if (log) -Inf else 0)
  }

  # Precompute accumulator parameters
  nu0 <- mu
  nu1 <- mu + mudelta
  sigma0 <- sigmazero
  sigma1 <- sigmazero * exp(sigmadelta)

  # Compute log-density for the winning accumulator
  log_pdf_win <- if (response == 0) {
    stats::dlnorm(t_adj, meanlog = nu0, sdlog = sigma0, log = TRUE)
  } else {
    stats::dlnorm(t_adj, meanlog = nu1, sdlog = sigma1, log = TRUE)
  }

  # Compute log-survival probability for the losing accumulator
  log_cdf_loss <- if (response == 0) {
    stats::plnorm(t_adj, meanlog = nu1, sdlog = sigma1, lower.tail = TRUE, log.p = TRUE)
  } else {
    stats::plnorm(t_adj, meanlog = nu0, sdlog = sigma0, lower.tail = TRUE, log.p = TRUE)
  }

  # Use .log1m_exp for numerically stable computation
  log_surv_loss <- .log1m_exp(log_cdf_loss)

  # Combine log-density and log-survival probability
  log_density <- log_pdf_win + log_surv_loss

  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
