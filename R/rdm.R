#' @title Simulate from the Two-Accumulator Racing Diffusion Model (RDM)
#'
#' @description
#' Simulates choice reaction times from a two-accumulator Racing Diffusion Model (RDM).
#' This is a specialized version where exactly two accumulators race towards a
#' common threshold. The model assumes variability in the starting point of the
#' diffusion process, drawn from a uniform distribution. This version is optimized
#' for performance using vectorized operations and allows one (but not both)
#' drift rate to be zero.
#'
#' @details
#' The RDM implemented here follows the formulation where the two accumulators
#' have drift rates `vzero` and `vone`. The diffusion process starts at a point `z` drawn
#' from `Uniform(0, A)`. The process terminates when either accumulator reaches
#' a threshold `b`. The parameter `k` is defined as `k = b - A`, representing
#' the distance from the maximum starting point `A` to the threshold `b`.
#' Therefore, the effective distance to threshold for a given trial is
#' `alpha = b - z = k + A - z`.
#'
#' The finishing time for a single accumulator, given its drift rate `v`, `k`, `A`, and `ndt`,
#' is simulated by:
#' 1. Sampling a starting point `z ~ Uniform(0, A)`.
#' 2. Calculating the distance `alpha = k + A - z`.
#' 3. If `v > 0`, simulating the time to reach `alpha` from an Inverse Gaussian distribution
#'    with mean `alpha / v` and shape `alpha^2`. This simulation uses an
#'    internal implementation based on Michael et al. (1976).
#' 4. If `v = 0`, the finishing time is considered infinite (`Inf`).
#' 5. Adding the non-decision time `ndt` to finite finishing times.
#'
#' The function simulates this process for both accumulators using vectorized operations.
#' The accumulator that finishes first determines the choice (1 for `vzero`, 2 for `vone`)
#' and the reaction time (RT) for that trial. If one drift rate is zero, the
#' accumulator with the positive drift rate will always win.
#'
#' This implementation is based on the description and parameters used in:
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#' without random between-trial variability: The racing diffusion model of
#' speeded decision making. *Psychonomic Bulletin & Review*, *27*, 911-936.
#' \doi{10.3758/s13423-020-01738-8} (specifically matching the `WaldA` component
#' used within their RDM simulation).
#'
#' @param n Number of trials to simulate. Must be a positive integer.
#' @param vzero Drift rate for the first accumulator. Must be a single non-negative number.
#' @param vone Drift rate for the second accumulator. Must be a single non-negative number.
#'   At least one of `vzero` or `vone` must be positive.
#' @param k Threshold parameter, defined as `k = b - A`, where `b` is the
#'   decision threshold and `A` is the maximum starting point. Must be a single
#'   positive number.
#' @param A Maximum starting point parameter. The starting point `z` for each
#'   accumulator on each trial is drawn from `Uniform(0, A)`. Must be a single
#'   positive number.
#' @param ndt Non-decision time (encoding and motor time offset). Must be a
#'   single non-negative number.
#'
#' @return A data frame with `n` rows and two columns:
#'   \item{rt}{The simulated reaction time (minimum finishing time across the two accumulators).}
#'   \item{choice}{The index of the winning accumulator (1 for `vzero`, 2 for `vone`).}
#'
#' @references
#' - Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating Random Variates Using
#'     Transformations with Multiple Roots. *The American Statistician*, *30*(2), 88–90. \doi{10.2307/2683801}
#' - Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#'     without random between-trial variability: The racing diffusion model of
#'     speeded decision making. *Psychonomic Bulletin & Review*, *27*, 911-936.
#'     \doi{10.3758/s13423-020-01738-8}
#' - Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian distribution and its
#'     statistical application—a review. *Journal of the Royal Statistical Society Series B:
#'     Statistical Methodology*, *40*(3), 263-275.
#'
#' @seealso `rshifted_wald`
#'
#' @examples
#' # Basic simulation with two positive drifts
#' sim_data_pos <- rrdm(n = 1000, vzero = 0.8, vone = 0.6, k = 0.5, A = 0.2, ndt = 0.15)
#' head(sim_data_pos)
#' mean(sim_data_pos$choice == 1) # Proportion of choices for accumulator 1
#' hist(sim_data_pos$rt, breaks = 30, main = "Simulated RDM RTs (v>0)")
#'
#' # Simulation with one drift rate = 0
#' sim_data_zero <- rrdm(n = 1000, vzero = 0.8, vone = 0, k = 0.5, A = 0.2, ndt = 0.15)
#' head(sim_data_zero)
#' mean(sim_data_zero$choice == 1) # Should be 1 (or very close due to floating point)
#' hist(sim_data_zero$rt, breaks = 30, main = "Simulated RDM RTs (vone=0)")
#' # Note: RT distribution differs from shifted_wald due to starting point variability
#'
#' @export
rrdm <- function(n, vzero, vone, k, A, ndt) {
  # Prepare and validate parameters using the helper function
  # This ensures all parameters are recycled to length 'n' (stored in params$n)
  params <- .prepare_rrdm(n = n, vzero = vzero, vone = vone, k = k, A = A, ndt = ndt)
  nobs <- params$ndraws

  # Handle special cases where one drift rate is zero across all trials
  # Note: This simplifies by simulating only the non-zero accumulator with a fixed threshold k+A
  if (all(params$vzero == 0)) {
    # Calculate alpha based on k and A (no z variability in this simplified case)
    alpha_fixed <- params$k + params$A
    return(data.frame(
      rt = rshifted_wald(n = nobs, nu = params$vone, alpha = alpha_fixed, ndt = params$ndt),
      choice = rep(2L, nobs)
    ))
  }
  if (all(params$vone == 0)) {
    # Calculate alpha based on k and A (no z variability in this simplified case)
    alpha_fixed <- params$k + params$A
    return(data.frame(
      rt = rshifted_wald(n = nobs, nu = params$vzero, alpha = alpha_fixed, ndt = params$ndt),
      choice = rep(1L, nobs)
    ))
  }

  # --- Simulation for cases where both drifts can be non-zero ---
  num_accumulators_total <- nobs * 2

  # Simulate starting points (z) for all trials and accumulators
  # params$A is already recycled to length nobs
  # We need 2*nobs starting points, one for each accumulator in each trial
  z_vec <- stats::runif(num_accumulators_total, min = 0, max = rep(params$A, each = 2))

  # Calculate distance to threshold (alpha = b - z = k + A - z)
  # Interleave k and A vectors to match the structure (trial1_acc1, trial1_acc2, ...)
  k_interleaved <- rep(params$k, each = 2)
  A_interleaved <- rep(params$A, each = 2)
  alpha_vec <- k_interleaved + A_interleaved - z_vec # Length = nobs * 2

  # Prepare drift rate vector (interleaving vzero_i, vone_i for each trial i)
  v_matrix <- cbind(params$vzero, params$vone) # nobs x 2 matrix
  v_vec <- as.vector(t(v_matrix)) # Length = nobs * 2, ordered (v0_1, v1_1, v0_2, v1_2, ...)

  # Prepare non-decision time vector (interleaving ndt_i, ndt_i for each trial i)
  ndt_vec_interleaved <- rep(params$ndt, each = 2) # Length = nobs * 2

  # Simulate finishing times using rshifted_wald
  # Pass the correctly sized and structured vectors
  finish_time_vec <- rshifted_wald(n = num_accumulators_total,
                                   nu = v_vec,
                                   alpha = alpha_vec,
                                   ndt = ndt_vec_interleaved)

  # Ensure the length matches before reshaping (should pass now)
  if (length(finish_time_vec) != num_accumulators_total) {
    stop("Internal error: Length of finish_time_vec does not match expected size after rshifted_wald.")
  }

  # Reshape into a matrix: nobs rows, 2 columns (acc1, acc2)
  finish_time_matrix <- matrix(finish_time_vec, nrow = nobs, ncol = 2, byrow = TRUE)

  # Extract times for each accumulator
  ft_acc1 <- finish_time_matrix[, 1]
  ft_acc2 <- finish_time_matrix[, 2]

  # Find minimum RT using pmin
  rt <- pmin(ft_acc1, ft_acc2)

  # Determine choice
  choice <- ifelse(ft_acc1 < ft_acc2, 1L, 2L)

  # --- Return Results ---
  data.frame(rt = rt, choice = choice)
}



#' @keywords internal
.prepare_rrdm <- function(n = NULL, vzero, vone, k, A, ndt) {
  # Validate parameters
  if (any(vzero < 0)) stop("Drift rate 'vzero' must be non-negative.")
  if (any(vone < 0)) stop("Drift rate 'vone' must be non-negative.")
  if (all(vzero == 0 & vone == 0)) stop("At least one drift rate (vzero or vone) must be positive.")
  if (any(k <= 0)) stop("Threshold parameter 'k' must be positive.")
  if (any(A <= 0)) stop("Maximum starting point 'A' must be positive.")
  if (any(ndt < 0)) stop("Non-decision time 'ndt' must be non-negative.")

  # Determine target length
  if (!is.null(n)) {
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n
  } else {
    stop("Parameter 'n' must be provided for simulation.")
  }

  # Recycle parameters to length m
  params <- list(
    vzero = rep_len(vzero, m),
    vone  = rep_len(vone,  m),
    k     = rep_len(k,     m),
    A     = rep_len(A,     m),
    ndt   = rep_len(ndt,   m)
  )

  params$ndraws <- m
  params
}