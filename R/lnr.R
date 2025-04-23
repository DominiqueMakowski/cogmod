#' @title Log-Normal Race (LNR) Model
#'
#' @description
#' The Log-Normal Race (LNR) model is useful for modeling reaction times and choices in decision-making tasks.
#' The model assumes that each choice option (accumulator) draws a processing time from a LogNormal distribution.
#' The winning accumulator (minimum draw) determines the observed reaction time and choice.
#' The observed RT includes a non-decision time component (ndt).
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
#' @param nuzero,nuone The (inverse of the) log-space mean parameter for both accumulators (choice 0 and 1).
#'   Controls the central tendency of the reaction time. Can take any real value (-Inf, Inf), with larger
#'   values leading to *faster* RTs. Named 'nu' (=-meanlog) for consistency with other race models.
#' @param sigmazero,sigmaone The log-space standard deviation for both accumulators (choice
#'   0 and 1) Controls the variability of reaction times. Must be positive (0, Inf).
#'   Larger values increase variability.
#' @param ndt Non-decision time (shift parameter). Represents the time taken for processes
#'   unrelated to the decision (e.g., encoding, motor response). Must be non-negative [0, Inf).
#'
#' @details
#' The LNR model conceptualizes decision-making as a race between two independent accumulators,
#' each corresponding to a potential choice. Each accumulator's finishing time is drawn from a
#' LogNormal distribution. The underlying `meanlog` parameter for the LogNormal distribution
#' for accumulator 0 is `-nuzero`, and for accumulator 1 is `-nuone`. The `sdlog` parameters
#' are `sigmazero` and `sigmaone` respectively.
#' The first accumulator to finish determines the choice and the decision time.
#' The observed reaction time (RT) is the decision time plus the non-decision time (ndt).
#' Higher values of `nuzero` or `nuone` correspond to faster processing speeds and thus
#' shorter reaction times for the respective accumulator.
#'
#' @references
#' - Rouder, J. N., Province, J. M., Morey, R. D., Gomez, P., & Heathcote, A. (2015).
#'   The lognormal race: A cognitive-process model of choice and latency with desirable
#'   psychometric properties. Psychometrika, 80(2), 491-513.
#'
#' @examples
#' # Simulate data
#' data <- rlnr(1000, nuzero = 1, nuone = 0.5, sigmazero = 1, sigmaone = 0.8, ndt = 0.2)
#' hist(data[data$response == 0, "rt"], breaks = 50, main = "Reaction Times", xlab = "RT")
#' hist(data[data$response == 1, "rt"], breaks = 50, add = TRUE, col = rgb(1, 0, 0, 0.5))
#'
#' @export
rlnr <- function(n, nuzero = 0, nuone = 0, sigmazero = 1, sigmaone = 1, ndt = 0.2) {
  # --- Prepare and Validate Parameters ---
  params <- .prepare_lnr(n = n, nuzero = nuzero, nuone = nuone,
                         sigmazero = sigmazero, sigmaone = sigmaone, ndt = ndt)
  n_out <- params$ndraws # Use the determined output length

  # --- Simulation ---
  # Generate log-normal draws for both accumulators (vectorized)
  # Use negative nu as meanlog for faster processing with higher nu
  draws0 <- stats::rlnorm(n_out, meanlog = -params$nuzero, sdlog = params$sigmazero)
  draws1 <- stats::rlnorm(n_out, meanlog = -params$nuone, sdlog = params$sigmaone)

  # Determine responses and reaction times (vectorized)
  response <- ifelse(draws0 < draws1, 0, 1)
  rt_decision <- pmin(draws0, draws1) # Minimum of the two decision times
  rt <- rt_decision + params$ndt # Add non-decision time

  # Return data frame matching the determined output length
  data.frame(rt = rt, response = response)
}




#' @rdname rlnr
#' @param x The observed reaction time (RT). Must be greater than `ndt`.
#' @param response The decision indicator (0 or 1). 0 for choice 0, 1 for choice 1.
#' @param log Logical; if TRUE, returns the log-density. Default: FALSE.
#' @export
dlnr <- function(x, nuzero, nuone, sigmazero, sigmaone, ndt, response, log = FALSE) {
  eps <- 1e-9 # Tolerance for checking RT > ndt

  # --- Prepare and Validate Parameters ---
  # Note: Basic parameter validation (positivity, etc.) happens in .prepare_lnr
  # Response validation (0/1) also happens there.
  # We still need to handle cases where parameters might lead to NA/Inf density later.
  params <- tryCatch(
      .prepare_lnr(x = x, nuzero = nuzero, nuone = nuone,
                   sigmazero = sigmazero, sigmaone = sigmaone,
                   ndt = ndt, response = response),
      error = function(e) {
          # If basic validation fails (e.g., negative sigma), return 0/-Inf density
          # Need to determine length based on x if possible, otherwise assume 1
          len <- if (!is.null(x)) length(x) else 1
          warning(conditionMessage(e), ". Returning 0 density / -Inf log-density.")
          return(list(ndraws = len, error = TRUE)) # Signal error state
      }
  )

  # If .prepare_lnr signaled an error, return appropriate value
  if (!is.null(params$error) && params$error) {
      return(rep(ifelse(log, -Inf, 0), params$ndraws))
  }

  n <- params$ndraws # Use the determined length

  # --- Calculation ---
  # Compute adjusted reaction times
  t_adj <- params$x - params$ndt

  # Initialize log-density vector with -Inf for invalid RTs (t_adj < eps)
  # or where parameters were invalid (e.g., sigma <= 0, caught by tryCatch earlier)
  log_density <- ifelse(t_adj < eps, -Inf, 0) # Use 0 for now, will add log-densities

  # Identify valid RTs for calculation
  valid_idx <- which(t_adj >= eps)
  if (length(valid_idx) == 0) { # No valid RTs
      return(ifelse(log, log_density, exp(log_density)))
  }

  # Subset parameters and adjusted times for valid indices
  t_adj_valid <- t_adj[valid_idx]
  # Use negative nu as meanlog
  meanlog0_valid <- -params$nuzero[valid_idx]
  meanlog1_valid <- -params$nuone[valid_idx]
  sigmazero_valid <- params$sigmazero[valid_idx]
  sigmaone_valid <- params$sigmaone[valid_idx]
  response_valid <- params$response[valid_idx]

  # Identify winning and losing parameters based on response
  meanlog_win <- ifelse(response_valid == 0, meanlog0_valid, meanlog1_valid)
  sigma_win <- ifelse(response_valid == 0, sigmazero_valid, sigmaone_valid)
  meanlog_loss <- ifelse(response_valid == 0, meanlog1_valid, meanlog0_valid)
  sigma_loss <- ifelse(response_valid == 0, sigmaone_valid, sigmazero_valid)

  # Compute log-density for the winning accumulator (vectorized)
  log_pdf_win <- stats::dlnorm(t_adj_valid, meanlog = meanlog_win, sdlog = sigma_win, log = TRUE)

  # Compute log-survival probability for the losing accumulator (vectorized)
  log_cdf_loss <- stats::plnorm(t_adj_valid, meanlog = meanlog_loss, sdlog = sigma_loss, lower.tail = TRUE, log.p = TRUE)
  log_surv_loss <- .log1m_exp(log_cdf_loss) # Stable log(1 - p)

  # Combine log-density and log-survival probability for valid indices
  log_density[valid_idx] <- log_pdf_win + log_surv_loss

  # Handle potential -Inf results from dlnorm/plnorm if parameters were invalid
  log_density[is.na(log_density) | !is.finite(log_density)] <- -Inf

  # Return result
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}


# Interals ---------------------------------------------------

#' @keywords internal
.prepare_lnr <- function(n = NULL, x = NULL, nuzero, nuone, sigmazero, sigmaone, ndt, response = NULL) {
  # --- Basic Validation ---
  if (any(sigmazero <= 0, na.rm = TRUE)) stop("sigmazero must be positive.")
  if (any(sigmaone <= 0, na.rm = TRUE)) stop("sigmaone must be positive.")
  if (any(ndt < 0, na.rm = TRUE)) stop("ndt must be non-negative.")
  # nu can be any real number, no check needed.

  # --- Determine Target Length ---
  if (!is.null(n)) {
    # For RNG (rlnr)
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n
    # Check lengths of parameters relative to n only if n > 1
    if (n > 1) {
        param_lengths <- c(length(nuzero), length(nuone), length(sigmazero),
                           length(sigmaone), length(ndt))
        m <- max(n, param_lengths)
    } else {
        m <- 1 # If n=1, output length is 1
    }

  } else if (!is.null(x)) {
    # For Density/Likelihood (dlnr)
    if (is.null(response)) stop("response must be provided for dlnr.")
    if (any(!response %in% c(0, 1), na.rm = TRUE)) stop("response must contain only 0 or 1.")

    param_lengths <- c(length(x), length(nuzero), length(nuone), length(sigmazero),
                       length(sigmaone), length(ndt), length(response))
    m <- max(param_lengths)
    if (m == 0) stop("At least one input vector (x, response, parameters) must have non-zero length.")

  } else {
    stop("Internal error: Either 'n' or 'x' must be provided.")
  }

  # --- Recycle Parameters ---
  params <- list(
    nuzero    = rep_len(nuzero, m),
    nuone     = rep_len(nuone, m),
    sigmazero = rep_len(sigmazero, m),
    sigmaone  = rep_len(sigmaone, m),
    ndt       = rep_len(ndt, m)
  )

  # --- Add x and response if provided ---
  if (!is.null(x)) {
    params$x <- rep_len(x, m)
  }
  if (!is.null(response)) {
    params$response <- rep_len(response, m)
  }

  params$ndraws <- m
  params
}
