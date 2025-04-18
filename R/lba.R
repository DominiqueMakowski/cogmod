#' @title Linear Ballistic Accumulator (LBA) Model Simulation
#'
#' @description
#' Simulates random draws (reaction times and choices) from a two-choice Linear Ballistic Accumulator (LBA) model.
#' This version uses a parameterization where the second accumulator's parameters are defined relative to the first.
#' The LBA model assumes that evidence for each choice option accumulates linearly and independently
#' until one accumulator reaches a threshold. The start point of accumulation is variable, drawn
#' from a uniform distribution `[0, A]`. Drift rates are also variable, drawn from a normal distribution.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param vzero Mean drift rate for the first accumulator (accumulator 0). Range: (-Inf, Inf).
#' @param vdelta Additive deviation for the mean drift rate of accumulator 1 (`v1 = vzero + vdelta`). Range: (-Inf, Inf).
#' @param sigmazero Standard deviation of the drift rate for the first accumulator (accumulator 0). Must be positive. Range: (0, Inf).
#' @param sigmadelta Log-deviation for the standard deviation of accumulator 1 (`sigma1 = sigmazero * exp(sigmadelta)`). Range: (-Inf, Inf).
#' @param A Maximum start point for the uniform distribution of starting evidence `[0, A]`. Must be positive. Range: (0, Inf). Default: 0.8.
#' @param k Difference between the decision threshold `b` and the maximum start point `A` (i.e., `b = A + k`). Must be positive. Range: (0, Inf). Default: 0.2.
#' @param ndt Non-decision time (shift parameter). Represents time for processes like encoding and motor response. Must be non-negative. Range: [0, Inf). Default: 0.3.
#'
#' @return A data frame with `n` rows and columns:
#'   - `rt`: Simulated reaction time.
#'   - `response`: The index (0 or 1) of the winning accumulator/choice.
#'
#' @details
#' The LBA model assumes that evidence for each choice option accumulates linearly and independently
#' until one accumulator reaches a threshold. The start point of accumulation is variable, drawn
#' from a uniform distribution `[0, A]`. Drift rates are also variable, drawn from a normal distribution.
#'
#' The simulation follows the standard LBA process:
#' 1. For each trial:
#'    - Sample drift rates `rate_0` from `Normal(vzero, sigmazero)` and `rate_1` from `Normal(vzero + vdelta, sigmazero * exp(sigmadelta))`.
#'    - Sample start points `start_0` and `start_1` from `Uniform(0, A)`.
#' 2. Resample drift rates for a trial if *both* sampled rates are non-positive (<= 0).
#' 3. Calculate the time `dt_i` for each accumulator to reach the threshold `b = A + k`: `dt_i = (b - start_i) / rate_i`. If `rate_i <= 0`, `dt_i` is effectively infinity.
#' 4. The accumulator with the minimum positive `dt_i` determines the choice (`response`) and the decision time (`min(dt_i)`).
#' 5. The final reaction time is `rt = min(dt_i) + ndt`.
#'
#' **Psychological Interpretation:**
#' - **Drift Rate (`vzero`, `vdelta`)**: Reflects the rate at which evidence accumulates for each choice. Higher drift rates indicate faster 
#'   evidence accumulation and a higher likelihood of selecting the corresponding choice. Differences in drift rates between accumulators 
#'   (via `vdelta`) can represent differences in preference, difficulty, or bias between the two options.
#' - **Drift Rate Variability (`sigmazero`, `sigmadelta`)**: Captures trial-to-trial variability in the evidence accumulation process. 
#'   Higher variability indicates less consistent evidence accumulation, leading to greater variability in reaction times and choices.
#' - **Start Point Variability (`A`)**: Represents the range of initial evidence levels for each accumulator. Larger values of `A` introduce 
#'   more variability in reaction times, as the starting point can vary more widely between trials.
#' - **Threshold (`b = A + k`)**: Represents the amount of evidence required to make a decision. Higher thresholds lead to longer reaction times 
#'   but more accurate decisions, as more evidence is required before a choice is made.
#' - **Non-Decision Time (`ndt`)**: Accounts for processes unrelated to evidence accumulation, such as sensory encoding and motor response. 
#'   This parameter shifts all reaction times by a constant amount.
#'
#' **Special Cases:**
#' - When `vdelta = 0` and `sigmadelta = 0`, the two accumulators are symmetric, meaning both choices are equally likely (assuming no bias in the start points or thresholds).
#' - When `A` is small relative to `k`, the model behaves more deterministically, as the start point variability has less influence on reaction times.
#' - When `sigmazero` or `sigmadelta` are large, the model produces more variable reaction times and less predictable choices.
#'
#' @references
#' Brown, S. D., & Heathcote, A. (2008). The simplest complete model of choice response time: Linear ballistic accumulation. *Cognitive Psychology*, *57*(3), 153-178. \doi{10.1016/j.cogpsych.2007.12.002}
#'
#' @examples
#' lba_data_rel <- rlba(n = 1000, vzero = 3, vdelta = -0.5,
#'                      sigmazero = 0.5, sigmadelta = 0,
#'                      A = 0.5, k = 0.5, ndt = 0.3)
#' hist(lba_data_rel$rt[lba_data_rel$response == 0], breaks = 50,
#'      col = rgb(0,0,1,0.5), xlab = "RT", main = "LBA Simulation (Relative Params)")
#' hist(lba_data_rel$rt[lba_data_rel$response == 1], breaks = 50,
#'      col = rgb(1,0,0,0.5), add = TRUE)
#'
#' @export
rlba <- function(n, vzero = 3, vdelta = 0, sigmazero = 1, sigmadelta = 0, A = 0.5, k = 0.5, ndt = 0.3) {
  # --- Input Validation ---
  if (length(n) != 1 || n <= 0 || n != floor(n))
    stop("n must be a single positive integer.")
  if (A <= 0 || k <= 0 || sigmazero <= 0 || ndt < 0)
    stop("A, k, sigmazero must be positive, and ndt must be non-negative.")

  # --- Derived Parameters ---
  b <- A + k                              # Decision threshold
  v1 <- vzero + vdelta                    # Mean drift rate for accumulator 1
  sigma1 <- sigmazero * exp(sigmadelta)     # SD of drift rate for accumulator 1

  # --- Drift Rate Sampling ---
  rates_matrix <- matrix(NA_real_, nrow = n, ncol = 2)
  valid_trials <- rep(FALSE, n)
  while (any(!valid_trials)) {
    remaining <- which(!valid_trials)
    n_remaining <- length(remaining)
    # Sample drift rates for each accumulator
    rates_0 <- stats::rnorm(n_remaining, mean = vzero, sd = sigmazero)
    rates_1 <- stats::rnorm(n_remaining, mean = v1, sd = sigma1)
    rates_matrix[remaining, ] <- cbind(rates_0, rates_1)
    # Valid if at least one accumulator has a positive rate
    valid_trials[remaining] <- rowSums(rates_matrix[remaining, , drop = FALSE] > 0) > 0
  }

  # --- Starting Points ---
  start_points <- matrix(stats::runif(n * 2, min = 0, max = A), nrow = n, ncol = 2)

  # --- Decision Times Calculation ---
  decision_times <- (b - start_points) / rates_matrix
  # Set any trial with non-positive rate (in a given accumulator) to Inf for that accumulator.
  decision_times[rates_matrix <= 0] <- Inf

  # --- Determine Winner ---
  # which.min returns the index (1 or 2), convert to 0-based for consistency with other functions.
  choices <- apply(decision_times, 1, which.min) - 1
  min_decision_times <- apply(decision_times, 1, min)

  # --- Reaction Times ---
  rts <- ndt + min_decision_times

  # --- Return as Data Frame ---
  data.frame(rt = rts, response = choices)
}




#' The density function `dlba` calculates the likelihood of observing a specific
#' reaction time `x` and response `response`, given the LBA parameters. It is
#' based on the formulation by Brown & Heathcote (2008), where the likelihood
#' is the product of the probability density of the winning accumulator finishing
#' at time `t = x - ndt` and the probability (survival function) that the losing
#' accumulator has not finished by time `t`. The density is normalized by
#' `(1 - pnegative)`, where `pnegative` is the probability that both drift rates
#' are non-positive, to account for the resampling process in `rlba`.
#' @rdname rlba
#' @inheritParams rlnr
#' @export
dlba <- function(x, response,
                 vzero = 3, vdelta = 0,
                 sigmazero = 1, sigmadelta = 0,
                 A = 0.5, k = 0.5, ndt = 0.3, log = FALSE) {
  eps <- 1e-10
  
  # --- Input Validation ---
  if (any(sigmazero <= 0) || any(A <= 0) || any(k <= 0) || any(ndt < 0)) {
    warning("sigmazero, A, k must be positive and ndt non-negative. Returning 0 density / -Inf log-density.")
    return(if (log) -Inf else 0)
  }

  # --- Vectorization ---
  n <- max(length(x), length(response), length(vzero), length(vdelta),
           length(sigmazero), length(sigmadelta), length(A), length(k), length(ndt))
  x <- rep_len(x, n)
  response <- rep_len(response, n)
  vzero <- rep_len(vzero, n)
  vdelta <- rep_len(vdelta, n)
  sigmazero <- rep_len(sigmazero, n)
  sigmadelta <- rep_len(sigmadelta, n)
  A <- rep_len(A, n)
  k <- rep_len(k, n)
  ndt <- rep_len(ndt, n)
  
  # --- Adjusted Decision Time ---
  t_adj <- x - ndt
  
  # Initialize log-density vector
  log_density <- rep(-Inf, n)
  
  # Identify valid cases
  valid_rt_idx <- which(t_adj > eps)
  valid_resp_idx <- which(response %in% c(0, 1))
  valid_params_idx <- which(
    is.finite(vzero) & is.finite(vdelta) &
    is.finite(sigmazero) & is.finite(sigmadelta) &
    is.finite(A) & is.finite(k) & is.finite(ndt)
  )
  valid_idx <- Reduce(intersect, list(valid_rt_idx, valid_resp_idx, valid_params_idx))
  
  if (length(valid_idx) == 0) {
    return(if (log) log_density else exp(log_density))
  }
  
  # --- Extract parameters for valid indices ---
  t_adj_v   <- t_adj[valid_idx]
  response_v <- response[valid_idx]
  v0        <- vzero[valid_idx]
  v1        <- vzero[valid_idx] + vdelta[valid_idx]
  s0        <- sigmazero[valid_idx]
  s1        <- sigmazero[valid_idx] * exp(sigmadelta[valid_idx])
  A_v       <- A[valid_idx]
  k_v       <- k[valid_idx]
  b_v       <- A_v + k_v
  
  # Determine winning and losing accumulators
  v_win <- ifelse(response_v == 0, v0, v1)
  s_win <- ifelse(response_v == 0, s0, s1)
  v_loss <- ifelse(response_v == 0, v1, v0)
  s_loss <- ifelse(response_v == 0, s1, s0)
  
  # --- Winner PDF Calculation ---
  st_win <- s_win * t_adj_v
  st_win[st_win < eps] <- eps
  z1_win <- (b_v - A_v - v_win * t_adj_v) / st_win
  z2_win <- (b_v - v_win * t_adj_v) / st_win
  
  # Compute PDF using the standard LBA density formula.
  pdf_win <- (v_win * (stats::pnorm(z2_win) - stats::pnorm(z1_win)) +
              s_win * (stats::dnorm(z1_win) - stats::dnorm(z2_win))) / A_v
  pdf_win <- pmax(pdf_win, eps)
  
  # --- Loser Survival Function Calculation ---
  st_loss <- s_loss * t_adj_v
  st_loss[st_loss < eps] <- eps
  z1_loss <- (b_v - A_v - v_loss * t_adj_v) / st_loss
  z2_loss <- (b_v - v_loss * t_adj_v) / st_loss
  
  cdf_loss <- 1 + (1 / A_v) * (
    (b_v - A_v - v_loss * t_adj_v) * stats::pnorm(z1_loss) +
      st_loss * stats::dnorm(z1_loss) -
      (b_v - v_loss * t_adj_v) * stats::pnorm(z2_loss) -
      st_loss * stats::dnorm(z2_loss)
  )
  # Clamp CDF to [0,1]
  cdf_loss <- pmin(pmax(cdf_loss, 0), 1)
  surv_loss <- pmax(1 - cdf_loss, eps)
  
  # --- pnegative Correction ---
  pneg0 <- stats::pnorm(-v0 / s0)
  pneg1 <- stats::pnorm(-v1 / s1)
  pnegative <- pneg0 * pneg1
  correction <- -log(pmax(1 - pnegative, eps))
  
  # --- Total Log-Density ---
  log_dens_valid <- log(pdf_win) + log(surv_loss) + correction
  log_density[valid_idx] <- log_dens_valid
  
  return(if (log) log_density else exp(log_density))
}