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
#' from `Uniform(0, bias)`. The process terminates when either accumulator reaches
#' a threshold `b`. The parameter `bs` is defined as `bs = b - bias`, representing
#' the distance from the maximum starting point `bias` to the threshold `b`.
#' Therefore, the effective distance to threshold for a given trial is
#' `bs = b - z = bs + bias - z`.
#'
#' The finishing time for a single accumulator, given its drift rate `v`, `bs`, `bias`, and `ndt`,
#' is simulated by:
#' 1. Sampling a starting point `z ~ Uniform(0, bias)`.
#' 2. Calculating the distance `bs = bs + bias - z`.
#' 3. If `v > 0`, simulating the time to reach `bs` from an Inverse Gaussian distribution
#'    with mean `bs / v` and shape `bs^2`. This simulation uses an
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
#' @param bs Threshold parameter, defined as `bs = b - bias`, where `b` is the
#'   decision threshold and `bias` is the maximum starting point. Must be a single
#'   positive number.
#' @param bias Maximum starting point parameter. The starting point for each
#'   accumulator on each trial is drawn from `Uniform(0, bias)`. Must be a single
#'   positive number.
#' @param ndt Non-decision time (encoding and motor time offset). Must be a
#'   single non-negative number.
#'
#' @return bias data frame with `n` rows and two columns:
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
#' @seealso `rrt_invgaussian`
#' @examples
#' rdm_pos <- rrdm(n = 1000, vzero = 0.8, vone = 0.6, bs = 0.5, bias = 0.2, ndt = 0.15)
#'
#' @export
rrdm <- function(n, vzero, vone, bs, bias, ndt) {
  # Prepare and validate parameters using the helper function
  # This ensures all parameters are recycled to length 'n' (stored in params$n)
  params <- .prepare_rrdm(n = n, vzero = vzero, vone = vone, bs = bs, bias = bias, ndt = ndt)
  nobs <- params$ndraws

  # Handle special cases where one drift rate is zero across all trials
  # Note: This simplifies by simulating only the non-zero accumulator with a fixed threshold bs+bias
  if (all(params$vzero == 0)) {
    # Calculate alpha based on bs and bias (no z variability in this simplified case)
    alpha_fixed <- params$bs + params$bias
    return(data.frame(
      rt = rrt_invgaussian(n = nobs, drift = params$vone, bs = alpha_fixed, ndt = params$ndt),
      choice = rep(2L, nobs)
    ))
  }
  if (all(params$vone == 0)) {
    # Calculate alpha based on bs and bias (no z variability in this simplified case)
    alpha_fixed <- params$bs + params$bias
    return(data.frame(
      rt = rrt_invgaussian(n = nobs, drift = params$vzero, bs = alpha_fixed, ndt = params$ndt),
      choice = rep(1L, nobs)
    ))
  }

  # --- Simulation for cases where both drifts can be non-zero ---
  num_accumulators_total <- nobs * 2

  # Simulate starting points (z) for all trials and accumulators
  # params$bias is already recycled to length nobs
  # We need 2*nobs starting points, one for each accumulator in each trial
  z_vec <- stats::runif(num_accumulators_total, min = 0, max = rep(params$bias, each = 2))

  # Calculate distance to threshold (alpha = b - z = bs + bias - z)
  # Interleave bs and bias vectors to match the structure (trial1_acc1, trial1_acc2, ...)
  k_interleaved <- rep(params$bs, each = 2)
  A_interleaved <- rep(params$bias, each = 2)
  alpha_vec <- k_interleaved + A_interleaved - z_vec # Length = nobs * 2

  # Prepare drift rate vector (interleaving vzero_i, vone_i for each trial i)
  v_matrix <- cbind(params$vzero, params$vone) # nobs x 2 matrix
  v_vec <- as.vector(t(v_matrix)) # Length = nobs * 2, ordered (v0_1, v1_1, v0_2, v1_2, ...)

  # Prepare non-decision time vector (interleaving ndt_i, ndt_i for each trial i)
  ndt_vec_interleaved <- rep(params$ndt, each = 2) # Length = nobs * 2

  # Simulate finishing times using rrt_invgaussian
  # Pass the correctly sized and structured vectors
  finish_time_vec <- rrt_invgaussian(n = num_accumulators_total,
                                   drift = v_vec,
                                   bs = alpha_vec,
                                   ndt = ndt_vec_interleaved)

  # Ensure the length matches before reshaping (should pass now)
  if (length(finish_time_vec) != num_accumulators_total) {
    stop("Internal error: Length of finish_time_vec does not match expected size after rrt_invgaussian.")
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



#' @rdname rrdm
#' @inheritParams rrt_invgaussian
#' @export
drdm <- function(x, vzero, vone, bs, bias, ndt, log = FALSE) {  
  # Prepare and validate parameters
  n    <- length(x)
  params <- .prepare_rrdm(n = n, vzero = vzero,
                            vone = vone, bs = bs,
                            bias = bias, ndt = ndt)
  # One-accumulator cases: delegate to .dwald
  if (all(params$vzero == 0)) {
    return(.dwald(x, drift = params$vone,
                  bs = params$bs, bias = params$bias,
                  ndt = params$ndt, log = log))
  }
  if (all(params$vone == 0)) {
    return(.dwald(x, drift = params$vzero,
                  bs = params$bs, bias = params$bias,
                  ndt = params$ndt, log = log))
  }
  # Two-accumulator case
  # Expand x vector
  x_vec <- rep(x, length.out = params$ndraws)
  # Log-densities and log-survivals
  logf1 <- .dwald(x_vec, params$vzero, params$bs, params$bias,
                  params$ndt, log = TRUE)
  logf2 <- .dwald(x_vec, params$vone,  params$bs, params$bias,
                  params$ndt, log = TRUE)
  logS1 <- .pwald(x_vec, params$vzero, params$bs, params$bias,
                  params$ndt, lower.tail = FALSE, log.p = TRUE)
  logS2 <- .pwald(x_vec, params$vone,  params$bs, params$bias,
                  params$ndt, lower.tail = FALSE, log.p = TRUE)
  # Stable log-sum-exp combine
  a <- logf1 + logS2
  b <- logf2 + logS1
  both_inf <- is.infinite(a) & is.infinite(b) & a < 0 & b < 0
  dens_log <- numeric(length(a))
  dens_log[both_inf] <- -Inf
  idx <- !both_inf
  m <- pmax(a[idx], b[idx])
  dens_log[idx] <- m + log(exp(a[idx] - m) + exp(b[idx] - m))
  # Zero density for x <= ndt
  is_zero <- x <= params$ndt
  if (log) {
    dens_log[is_zero] <- -Inf
    # Underflow correction: ensure finite log-density for x > ndt
    dens_log[!is_zero & is.infinite(dens_log)] <- log(.Machine$double.xmin)
    return(dens_log)
  } else {
    dens <- exp(dens_log)
    dens[is_zero] <- 0
    # Underflow correction: ensure positive density for x > ndt
    zero_idx <- (!is_zero & dens == 0)
    dens[zero_idx] <- .Machine$double.xmin
    return(dens)
  }
}


# Internals ---------------------------------------------------------------


#' Analytical PDF for Wald with Uniform Start Point Variability
#'
#' Calculates the probability density function (PDF) for the time it takes a
#' diffusion process with drift `drift`, starting point `z ~ U(0, bias)`, to reach
#' threshold `b = bs + bias`, shifted by `ndt`. Based on Tillman et al. (2020).
#' This function calculates the density for the *unadjusted* time `x`.
#'
#' @keywords internal
.dwald <- function(x, drift, bs, bias, ndt, log = FALSE) {
  # 1. Quick length check
  n <- length(x)
  if (!all(lengths(list(drift, bs, bias, ndt)) == n)) {
    stop("All inputs must be numeric vectors of the same length.")
  }

  dens    <- rep(NA_real_, n)
  # 2. Identify invalid parameters (NA or non-positive bounds)
  invalid <- is.na(x) | is.na(drift) | is.na(bs) | is.na(bias) | is.na(ndt) |
             (bs <= 0) | (bias <= 0) | (ndt < 0)
  dens[invalid] <- NA_real_

  # 3. Compute adjusted times and baseline zeros
  ok       <- !invalid
  t_adj    <- x[ok] - ndt[ok]
  zero_pdf <- (t_adj <= 0) | (drift[ok] < 0)
  dens_ok  <- numeric(sum(ok))
  dens_ok[zero_pdf] <- 0

  # 4. drift == 0 case (Eq.6: no drift)
  is0 <- !zero_pdf & (drift[ok] == 0)
  if (any(is0)) {
    t0      <- t_adj[is0]
    sqrt_t0 <- sqrt(t0)
    k0      <- bs[ok][is0];  A0 <- bias[ok][is0]
    alpha   <- k0 / sqrt_t0         # (bs - 0*t)/sqrt(t)
    beta    <- (k0 + A0) / sqrt_t0  # (bs+bias - 0*t)/sqrt(t)
    # g0(t) = (phi(alpha) - phi(beta)) / (bias*sqrt(t))
    dens_ok[is0] <- (stats::dnorm(alpha) - stats::dnorm(beta)) / (A0 * sqrt_t0)
  }

  # 5. drift > 0 case (Eq.5)
  pos <- !zero_pdf & (drift[ok] > 0)
  if (any(pos)) {
    t1      <- t_adj[pos];    sqrt_t1 <- sqrt(t1)
    nu1     <- drift[ok][pos];   k1 <- bs[ok][pos];   A1 <- bias[ok][pos]
    alpha   <- (k1 - nu1 * t1) / sqrt_t1       # see alpha definition
    beta    <- (k1 + A1 - nu1 * t1) / sqrt_t1  # see beta definition
    Phi_a   <- stats::pnorm(alpha);   Phi_b <- stats::pnorm(beta)
    phi_a   <- stats::dnorm(alpha);   phi_b <- stats::dnorm(beta)
    invA    <- 1 / A1
    inv_st  <- 1 / sqrt_t1
    dens_ok[pos] <- invA * (
      - nu1 * Phi_a + inv_st * phi_a +
        nu1 * Phi_b - inv_st * phi_b
    )
  }

  # 6. Assemble densities, clamp, and optional log
  dens[ok] <- pmax(0, dens_ok)
  if (log) dens <- log(dens)
  dens
}


#' Analytical CDF for Wald with Uniform Start Point Variability
#'
#' Calculates the cumulative distribution function (CDF) for the time it takes a
#' diffusion process with drift `drift`, starting point `z ~ U(0, bias)`, to reach
#' threshold `b = bs + bias`, shifted by `ndt`. Based on Tillman et al. (2020).
#' This function calculates the CDF for the *unadjusted* time `x`.
#' 
#' @keywords internal
.pwald <- function(x, drift, bs, bias, ndt,
                  lower.tail = TRUE, log.p = FALSE) {
  # 1. Quick length check
  n <- length(x)
  if (!all(lengths(list(drift, bs, bias, ndt)) == n)) {
    stop("All inputs must be numeric vectors of the same length.")
  }

  C       <- rep(NA_real_, n)
  invalid <- is.na(x) | is.na(drift) | is.na(bs) | is.na(bias) | is.na(ndt) |
             (bs <= 0) | (bias <= 0) | (ndt < 0)
  C[invalid] <- NA_real_

  ok        <- !invalid
  t_adj     <- x[ok] - ndt[ok]
  zero_cdf  <- (t_adj <= 0) | (drift[ok] < 0)
  C_ok      <- numeric(sum(ok))
  C_ok[zero_cdf] <- 0

  # 2. drift == 0 branch (Appendix bias, v=0)
  is0 <- !zero_cdf & (drift[ok] == 0)
  if (any(is0)) {
    t0      <- t_adj[is0]; sqrt_t0 <- sqrt(t0)
    k0      <- bs[ok][is0]; A0 <- bias[ok][is0]
    b0      <- k0 + A0
    a1      <- -b0 / sqrt_t0  # alpha1
    a2      <- -k0 / sqrt_t0  # alpha2
    P1      <- stats::pnorm(a1); P2 <- stats::pnorm(a2)
    p1      <- stats::dnorm(a1); p2 <- stats::dnorm(a2)
    # C0 = (sqrt(t)/bias)[a2*P2 - a1*P1 + p2 - p1]
    C_ok[is0] <- (sqrt_t0/A0) * (a2*P2 - a1*P1 + p2 - p1)
  }

  # 3. drift > 0 branch (Appendix bias full formula)
  pos <- !zero_cdf & (drift[ok] > 0)
  if (any(pos)) {
    t1 <- t_adj[pos]; sqrt_t1 <- sqrt(t1)
    nu1 <- drift[ok][pos]; k1 <- bs[ok][pos]; A1 <- bias[ok][pos]
    b1  <- k1 + A1
    a1  <- (nu1*t1 - b1) / sqrt_t1
    a2  <- (nu1*t1 - k1) / sqrt_t1
    bth1<- (-nu1*t1 - b1)/sqrt_t1
    bth2<- (-nu1*t1 - k1)/sqrt_t1
    P_a1<- stats::pnorm(a1); P_a2<- stats::pnorm(a2)
    P_b1<- stats::pnorm(bth1);P_b2<- stats::pnorm(bth2)
    p_a1<- stats::dnorm(a1); p_a2<- stats::dnorm(a2)
    inv2vA <- 1/(2*nu1*A1)
    stA    <- sqrt_t1/A1
    C_ok[pos] <- (
      inv2vA*(P_a2 - P_a1) +
      stA*(a2*P_a2 - a1*P_a1) -
      inv2vA*(exp(2*nu1*k1)*P_b2 - exp(2*nu1*b1)*P_b1) +
      stA*(p_a2 - p_a1)
    )
  }

  # 4. Clamp to [0,1]; tail and log options
  C[ok] <- pmin(1, pmax(0, C_ok))
  if (!lower.tail) C <- 1 - C
  if (log.p)       C <- log(C)
  C
}



#' @keywords internal
.prepare_rrdm <- function(n = NULL, vzero, vone, bs, bias, ndt) {
  # Validate parameters
  if (any(vzero < 0)) stop("Drift rate 'vzero' must be non-negative.")
  if (any(vone < 0)) stop("Drift rate 'vone' must be non-negative.")
  if (all(vzero == 0 & vone == 0)) stop("At least one drift rate (vzero or vone) must be positive.")
  if (any(bs <= 0)) stop("Threshold parameter 'bs' must be positive.")
  if (any(bias <= 0)) stop("Maximum starting point 'bias' must be positive.")
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
    bs     = rep_len(bs,     m),
    bias     = rep_len(bias,     m),
    ndt   = rep_len(ndt,   m)
  )

  params$ndraws <- m
  params
}