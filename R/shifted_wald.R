#' @title Shifted Wald (Shifted Inverse Gaussian) Model Functions
#'
#' @description
#' Density, distribution function, and random generation for the Shifted Wald
#' distribution, also known as the Shifted Inverse Gaussian distribution. This
#' distribution is commonly used in modeling reaction times in cognitive tasks.
#' It is characterized by a drift rate (`nu`), a decision threshold (`alpha`),
#' and a non-decision time (`ndt`).
#'
#' Functions:
#' - `rshifted_wald()`: Simulates random draws from the Shifted Wald model.
#' - `dshifted_wald()`: Computes the density (likelihood) of the Shifted Wald distribution.
#' - `pshifted_wald()`: Computes the cumulative distribution function (CDF).
#' - `qshifted_wald()`: Computes the quantile function (inverse CDF).
#'
#' @details
#' The Shifted Wald distribution describes the time it takes for a Wiener diffusion
#' process starting at 0 to reach a threshold `alpha` > 0, given a positive drift
#' rate `nu` > 0. The resulting time is then shifted by a non-decision time `ndt` >= 0.
#'
#' The distribution is mathematically equivalent to shifting an Inverse Gaussian
#' distribution with mean `mu = alpha / nu` and shape parameter `lambda = alpha^2`.
#' That is, `ShiftedWald(nu, alpha, ndt) = InverseGaussian(mean = alpha/nu, shape = alpha^2) + ndt`.
#'
#' The random generation algorithm implemented here is based on the method described
#' by Michael, Schucany, and Haas (1976), as used in the `statmod` package.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param nu Drift rate. Must be positive. Represents the average speed of evidence accumulation.
#'   Range: (0, Inf).
#' @param alpha Decision threshold. Must be positive. Represents the amount of evidence needed
#'   to make a decision. Range: (0, Inf).
#' @param ndt Non-decision time (shift parameter). Must be non-negative. Represents time for
#'   processes like stimulus encoding and response execution. Range: [0, Inf).
#'
#' @return `rshifted_wald` returns a vector of random deviates.
#'
#' @references
#' - Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating Random Variates Using
#'     Transformations with Multiple Roots. *The American Statistician*, *30*(2), 88–90. \doi{10.2307/2683801}
#' - Anders, R., Alario, F., & Van Maanen, L. (2016). The shifted Wald distribution for
#'     response time data analysis. *Psychological Methods*, *21*(3), 309–327. \doi{10.1037/met0000063}
#' - Matzke, D., & Wagenmakers, E. J. (2009). Psychological interpretation of the ex-Gaussian
#'     and shifted Wald parameters: A diffusion model analysis. *Psychonomic Bulletin & Review*,
#'     *16*(5), 798–817. \doi{10.3758/PBR.16.5.798}
#' - Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian distribution and its
#'     statistical application—a review. *Journal of the Royal Statistical Society Series B:
#'     Statistical Methodology*, *40*(3), 263-275.
#'
#' @seealso `dshifted_wald`, `pshifted_wald`, `qshifted_wald`
#'
#' @examples
#' # Simulate 1000 RTs
#' rts <- rshifted_wald(1000, nu = 3, alpha = 0.5, ndt = 0.2)
#' hist(rts, breaks = 50, main = "Simulated Shifted Wald RTs", xlab = "Reaction Time")
#'
#' @export
rshifted_wald <- function(n, nu = 3, alpha = 0.5, ndt = 0.2) {
  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) stop("n must be a positive integer.")
  if (any(nu <= 0)) stop("Drift rate 'nu' must be positive.")
  if (any(alpha <= 0)) stop("Threshold 'alpha' must be positive.")
  if (any(ndt < 0)) stop("Non-decision time 'ndt' must be non-negative.")

  # --- Vectorization ---
  param_lengths <- c(length(nu), length(alpha), length(ndt))
  n_params <- max(param_lengths)
  n_out <- max(n, n_params)
  if (n_out == 0) return(numeric(0))
  if (n_out > 1) {
      nu <- rep(nu, length.out = n_out)
      alpha <- rep(alpha, length.out = n_out)
      ndt <- rep(ndt, length.out = n_out)
  }

  # --- Simulation (Inverse Gaussian part, adapted from statmod::rinvgauss) ---
  ig_mean <- alpha / nu # mu parameter for IG
  ig_lambda <- alpha^2  # lambda (shape) parameter for IG
  ig_phi <- 1 / ig_lambda # dispersion = 1 / lambda

  # Handle cases where parameters might lead to non-finite results early
  valid_params <- is.finite(ig_mean) & ig_mean > 0 & is.finite(ig_phi) & ig_phi > 0
  shifted_draws <- numeric(n_out)
  shifted_draws[!valid_params] <- NA_real_

  if (!any(valid_params)) return(shifted_draws)

  n_valid <- sum(valid_params)
  mu_valid <- ig_mean[valid_params]
  phi_valid <- ig_phi[valid_params]

  # Generate chisquare(1) deviates
  y <- stats::rnorm(n_valid)^2

  # Calculate intermediate term y * phi * mu
  y_phi_mu <- y * phi_valid * mu_valid

  # Calculate x1/mu (denoted X1 in statmod code)
  x1_over_mu <- numeric(n_valid)

  # Check for large intermediate values (Taylor series approximation)
  # Threshold from statmod
  taylor_threshold <- 5e5
  use_taylor <- y_phi_mu > taylor_threshold
  use_standard <- !use_taylor

  # Standard calculation
  idx_standard <- which(use_standard) # Get indices for standard calculation
  if (length(idx_standard) > 0) {
      y_phi_mu_std <- y_phi_mu[idx_standard]
      # Handle y_phi_mu == 0 case within standard calculation
      non_zero_std <- y_phi_mu_std > 0
      zero_std <- !non_zero_std
      if (any(non_zero_std)) {
          y_phi_mu_std_nz <- y_phi_mu_std[non_zero_std]
          x1_over_mu[idx_standard[non_zero_std]] <- 1 + y_phi_mu_std_nz / 2 * (1 - sqrt(1 + 4 / y_phi_mu_std_nz))
      }
      if (any(zero_std)) {
          x1_over_mu[idx_standard[zero_std]] <- 1 # Limit as y_phi_mu -> 0
      }
  }

  # Taylor series approximation for large y_phi_mu
  idx_taylor <- which(use_taylor) # Get indices for Taylor calculation
  if (length(idx_taylor) > 0) {
      # Approximation: x1/mu ≈ 1 / y_phi_mu for large y_phi_mu
      # This comes from 1 + z/2 * (1 - sqrt(1+4/z)) ≈ 1 + z/2 * (1 - (1 + 2/z - 2/z^2)) = 1 + z/2 * (-2/z + 2/z^2) = 1 - 1 + 1/z = 1/z
      x1_over_mu[idx_taylor] <- 1 / y_phi_mu[idx_taylor]
  }

  # Generate uniform random numbers for root selection
  u <- stats::runif(n_valid)

  # Select the root based on the uniform draw
  ig_draws_over_mu <- ifelse(u < 1 / (1 + x1_over_mu), x1_over_mu, 1 / x1_over_mu)

  # Multiply by mu to get the IG draw
  ig_draws_valid <- mu_valid * ig_draws_over_mu

  # Add shift and place results back into the full vector
  shifted_draws[valid_params] <- ig_draws_valid + ndt[valid_params]

  # --- Return ---
  if (n_out == n && n > 0) {
      return(shifted_draws[1:n])
  } else if (n_out > n && n == 1) {
      return(shifted_draws)
  } else {
      return(shifted_draws)
  }
}




#' @rdname rshifted_wald
#' @param x Vector of quantiles (observed reaction times).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @return `dshifted_wald` returns the density (PDF).
#' @export
dshifted_wald <- function(x, nu = 3, alpha = 0.5, ndt = 0.2, log = FALSE) {
  # --- Input Validation ---
  if (any(nu <= 0)) stop("Drift rate 'nu' must be positive.")
  if (any(alpha <= 0)) stop("Threshold 'alpha' must be positive.")
  if (any(ndt < 0)) stop("Non-decision time 'ndt' must be non-negative.")

  # --- Vectorization ---
  # Determine output length based on x and parameter vector lengths
  param_lengths <- c(length(x), length(nu), length(alpha), length(ndt))
  n_out <- max(param_lengths)

  # Handle zero-length output
  if (n_out == 0) return(numeric(0))

  # Recycle parameters to match output length
  if (n_out > 1) {
      x <- rep(x, length.out = n_out)
      nu <- rep(nu, length.out = n_out)
      alpha <- rep(alpha, length.out = n_out)
      ndt <- rep(ndt, length.out = n_out)
  }

  # --- Density Calculation ---
  # Initialize result vector
  log_density <- numeric(n_out)

  # Calculate adjusted x (time relative to non-decision time)
  x_adj <- x - ndt

  # Identify valid cases for calculation (x > ndt and finite parameters)
  valid_params <- is.finite(nu) & is.finite(alpha) & is.finite(ndt)
  valid_x <- is.finite(x_adj) & x_adj > 0
  calculate_idx <- valid_params & valid_x

  # Set density to 0 (log-density to -Inf) for x <= ndt
  log_density[!valid_x & valid_params] <- -Inf

  # Set density to NA for invalid parameters or NA x
  log_density[!valid_params | is.na(x)] <- NA_real_

  # Proceed with calculation only for valid indices
  if (any(calculate_idx)) {
    x_adj_calc <- x_adj[calculate_idx]
    nu_calc <- nu[calculate_idx]
    alpha_calc <- alpha[calculate_idx]
    # ndt_calc <- ndt[calculate_idx] # Not needed directly in formula

    # Calculate Inverse Gaussian parameters for valid cases
    mu_ig <- alpha_calc / nu_calc
    lambda_ig <- alpha_calc^2

    # Calculate log-density using the standard Inverse Gaussian PDF formula
    # log(f(x; mu, lambda)) = 0.5 * (log(lambda) - log(2*pi) - 3*log(x)) - (lambda * (x - mu)^2) / (2 * mu^2 * x)
    term1 <- 0.5 * (log(lambda_ig) - log(2 * pi) - 3 * log(x_adj_calc))
    term2_num <- lambda_ig * (x_adj_calc - mu_ig)^2
    term2_den <- 2 * mu_ig^2 * x_adj_calc
    # Avoid division by zero if mu_ig or x_adj_calc are somehow zero (shouldn't happen with validation)
    term2 <- ifelse(term2_den > 0, term2_num / term2_den, Inf) # Assign Inf if denominator is zero

    log_density[calculate_idx] <- term1 - term2
  }

  # Return log or standard density
  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}
