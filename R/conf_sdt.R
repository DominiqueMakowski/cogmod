#' @title Confidence Signal Detection Theory Model (EXPERIMENTAL)
#'
#' @description
#' Functions for the Signal Detection Theory (SDT) model of confidence. This model
#' assumes that a single sample of sensory evidence is used to generate both a
#' discrimination response and a confidence rating.
#'
#' Functions:
#' - `rconf_sdt()`: Simulates random draws from the SDT model.
#' - `dconf_sdt()`: Computes the likelihood/density of observed responses and ratings.
#'
#' @param n Number of simulated trials.
#' @param dprime Sensitivity parameter(s) (d'). Can be a single value or a vector.
#' @param c Response bias parameter. Must be a single value.
#' @param thetazero A numeric vector of confidence criteria for response=0. Must be sorted in decreasing order.
#' @param thetaone A numeric vector of confidence criteria for response=1. Must be sorted in increasing order.
#' @param truth Stimulus value (0 or 1). For simulation, can be a vector of length `n`. If NULL, it's sampled. For density, it's the observed stimulus.
#' @param response Observed response (0 or 1). For density calculation.
#' @param confidence Observed confidence rating (integer >= 0). For density calculation.
#' @param log Logical; if TRUE, returns the log-density.
#'
#' @details
#' The SDT model assumes that for a given stimulus `truth` (0 or 1), sensory
#' evidence `x` is drawn from a normal distribution. When `truth` is 1, the
#' evidence is drawn from `N(dprime / 2, 1)`; when `truth` is 0, it's drawn
#' from `N(-dprime / 2, 1)`.
#' The `response` is 1 if `x > -c` and 0 otherwise. The decision criterion is thus `-c`.
#' The confidence rating is determined by which region the evidence `x` falls into,
#' defined by the confidence criteria `thetaone` and `thetazero`.
#' A rating of 0 corresponds to the lowest confidence (i.e., evidence closest to the decision criterion).
#'
#' **Parameter Interpretation:**
#' - `dprime`: The sensitivity (d'). It reflects the observer's ability to distinguish between the two stimulus categories. Higher values indicate better performance.
#' - `c`: The response bias. It reflects the observer's tendency to favor one response over the other, independent of the stimulus. A positive value indicates a bias towards response 1, and a negative value indicates a bias towards response 0. A value of 0 is unbiased.
#' - `thetazero` & `thetaone`: The confidence criteria. These thresholds partition the evidence space into different confidence levels for each response.
#'   - For `response = 1`, `thetaone` are the criteria. They must be greater than the decision criterion (`-c`). For a given response, higher confidence ratings correspond to values of `x` further away from the decision criterion. For example, `confidence = 0` if `-c < x < thetaone[1]`, `confidence = 1` if `thetaone[1] < x < thetaone[2]`, and so on.
#'   - For `response = 0`, `thetazero` are the criteria. They must be less than the decision criterion (`-c`). For example, `confidence = 0` if `thetazero[1] < x < -c`, `confidence = 1` if `thetazero[2] < x < thetazero[1]`, and so on.
#'
#' @references
#' - Green, D. M., & Swets, J. A. (1966). Signal detection theory and psychophysics. Wiley.
#'
#' @examples
#' # Simulate data with 2 criteria per response (3 confidence levels: 0, 1, 2)
#' sim_data <- rconf_sdt(
#'   n = 1000, dprime = 1, c = 0.2,
#'   thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
#' )
#' table(sim_data$response, sim_data$confidence)
#'
#' # Calculate density for confidence=2 (highest) for response=1
#' dconf_sdt(
#'   truth = 1, response = 1, confidence = 2, dprime = 1, c = 0.2,
#'   thetazero = c(-0.5, -1.5), thetaone = c(0.5, 1.5)
#' )
#'
#' @name conf_sdt
NULL

#' @rdname conf_sdt
#' @export
rconf_sdt <- function(n, dprime, c, thetazero, thetaone, truth = NULL) {
  # --- Prepare and Validate Parameters ---
  if (is.null(truth)) {
    truth <- sample(c(0, 1), n, replace = TRUE)
  }
  params <- .prepare_conf_sdt(n = n, dprime = dprime, c = c, truth = truth)
  if (length(c) > 1) stop("'c' must be a single value.")

  decision_criterion <- -c

  if (any(thetaone <= decision_criterion)) stop("All thetaone must be > -c.")
  if (any(thetazero >= decision_criterion)) stop("All thetazero must be < -c.")
  if (is.unsorted(thetaone)) stop("thetaone must be sorted in increasing order.")
  if (is.unsorted(rev(thetazero))) stop("thetazero must be sorted in decreasing order.")

  # --- Simulation ---
  # Generate sensory evidence
  mean_x <- ((params$truth * 2) - 1) * params$dprime / 2
  x <- stats::rnorm(params$ndraws, mean = mean_x, sd = 1)

  # Generate response (0 or 1)
  response <- ifelse(x > decision_criterion, 1, 0)

  # Generate confidence rating (0-indexed)
  confidence <- numeric(params$ndraws)

  idx_r1 <- which(response == 1)
  if (length(idx_r1) > 0) {
    # Confidence = 0 for -c < x < t1, 1 for t1 < x < t2, etc.
    # findInterval is 1-based, so we subtract 1 to get 0-based ratings.
    confidence[idx_r1] <- findInterval(x[idx_r1], vec = c(decision_criterion, thetaone)) - 1
  }

  idx_r0 <- which(response == 0)
  if (length(idx_r0) > 0) {
    # Confidence = 0 for tz1 < x < -c, 1 for tz2 < x < tz1, and so on.
    # This is equivalent to counting how many criteria in thetazero x is smaller than.
    confidence[idx_r0] <- rowSums(outer(x[idx_r0], thetazero, `<`))
  }

  # Return data frame
  data.frame(truth = params$truth, response = response, confidence = confidence)
}

#' @rdname conf_sdt
#' @export
dconf_sdt <- function(truth, response, confidence, dprime, c, thetazero, thetaone, log = FALSE) {
  # --- Prepare and Validate Parameters ---
  params <- .prepare_conf_sdt(truth = truth, response = response, confidence = confidence, dprime = dprime, c = c)
  if (length(c) > 1) stop("'c' must be a single value.")

  decision_criterion <- -c

  if (any(thetaone <= decision_criterion)) stop("All thetaone must be > -c.")
  if (any(thetazero >= decision_criterion)) stop("All thetazero must be < -c.")
  if (is.unsorted(thetaone)) stop("thetaone must be sorted in increasing order.")
  if (is.unsorted(rev(thetazero))) stop("thetazero must be sorted in decreasing order.")

  # --- Calculation ---
  lower_bound <- numeric(params$ndraws)
  upper_bound <- numeric(params$ndraws)

  # --- response == 1 ---
  idx_r1 <- which(params$response == 1)
  if (length(idx_r1) > 0) {
    conf_r1 <- params$confidence[idx_r1]
    if (any(conf_r1 > length(thetaone))) {
      stop("A confidence value for response=1 is higher than allowed by the number of thetaone criteria.")
    }
    boundaries <- c(decision_criterion, thetaone, Inf)
    lower_bound[idx_r1] <- boundaries[conf_r1 + 1]
    upper_bound[idx_r1] <- boundaries[conf_r1 + 2]
  }

  # --- response == 0 ---
  idx_r0 <- which(params$response == 0)
  if (length(idx_r0) > 0) {
    conf_r0 <- params$confidence[idx_r0]
    if (any(conf_r0 > length(thetazero))) {
      stop("A confidence value for response=0 is higher than allowed by the number of thetazero criteria.")
    }
    # Upper bound candidates: -c, tzero1, tzero2, ...
    upper_bound_cand <- c(decision_criterion, thetazero)
    upper_bound[idx_r0] <- upper_bound_cand[conf_r0 + 1]
    # Lower bound candidates: tzero1, tzero2, ..., -Inf
    lower_bound_cand <- c(thetazero, -Inf)
    lower_bound[idx_r0] <- lower_bound_cand[conf_r0 + 1]
  }

  # Calculate probability
  mean_x <- ((params$truth * 2) - 1) * params$dprime / 2
  prob <- stats::pnorm(upper_bound, mean = mean_x, sd = 1) - stats::pnorm(lower_bound, mean = mean_x, sd = 1)

  # Return
  if (log) {
    return(log(prob))
  } else {
    return(prob)
  }
}


#' @keywords internal
.prepare_conf_sdt <- function(n = NULL, response = NULL, truth = NULL, confidence = NULL, dprime, c) {
  # --- Determine Target Length ---
  if (!is.null(n)) {
    # For RNG (rconf_sdt)
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    param_lengths <- c(length(dprime), length(truth))
    m <- max(n, param_lengths)
  } else {
    # For Density/Likelihood (dconf_sdt)
    if (is.null(response) || is.null(confidence) || is.null(truth)) stop("response, truth and confidence must be provided for dconf_sdt.")
    if (any(!response %in% c(0, 1), na.rm = TRUE)) stop("response must contain only 0 or 1.")
    if (any(!truth %in% c(0, 1), na.rm = TRUE)) stop("truth must contain only 0 or 1.")
    if (any(confidence < 0, na.rm = TRUE)) stop("confidence must be >= 0.")

    param_lengths <- c(length(truth), length(response), length(confidence), length(dprime))
    m <- max(param_lengths)
    if (m == 0) stop("At least one input vector must have non-zero length.")
  }

  # --- Recycle Parameters ---
  params <- list(
    dprime = rep_len(dprime, m),
    c = c, # Keep as scalar, checked in main function
    truth = rep_len(truth, m)
  )

  # --- Add other variables if provided ---
  if (!is.null(response)) {
    params$response <- rep_len(response, m)
  }
  if (!is.null(confidence)) {
    params$confidence <- rep_len(confidence, m)
  }

  params$ndraws <- m
  params
}
