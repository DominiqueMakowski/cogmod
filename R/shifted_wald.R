#' @title Shifted Wald Model (Inverse Gaussian)
#'
#' @description
#' Density, distribution function, and random generation for the Shifted Wald
#' distribution, also known as the Shifted Inverse Gaussian distribution. This
#' distribution is commonly used in modeling reaction times in cognitive tasks.
#' It is characterized by a drift rate (`drift`), a decision threshold (`bs`),
#' and a non-decision time (`ndt`).
#'
#' Functions:
#' - `rshifted_wald()`: Simulates random draws from the Shifted Wald model.
#' - `dshifted_wald()`: Computes the density (likelihood) of the Shifted Wald distribution.
#' - `pshifted_wald()`: Computes the cumulative distribution function (CDF).
#'
#' @details
#' The Shifted Wald distribution describes the time it takes for a Wiener diffusion
#' process starting at 0 to reach a threshold `bs` > 0, given a positive drift
#' rate `drift` > 0. The resulting time is then shifted by a non-decision time `ndt` >= 0.
#'
#' The distribution is mathematically equivalent to shifting an Inverse Gaussian
#' distribution with mean `mu = bs / drift` and shape parameter `lambda = bs^2`.
#' That is, `ShiftedWald(drift, bs, ndt) = InverseGaussian(mean = bs/drift, shape = bs^2) + ndt`.
#'
#' The random generation algorithm implemented here is based on the method described
#' by Michael, Schucany, and Haas (1976), as used in the `statmod` package.
#'
#' @param n Number of observations. If `length(n) > 1`, the length is taken to be the number required.
#' @param drift Drift rate. Must be positive. Represents the average speed of evidence accumulation.
#'   Range: (0, Inf).
#' @param bs Decision threshold (boundary separation). Must be positive. Represents the amount of evidence needed
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
#' @examples
#' # Simulate 1000 RTs
#' rts <- rshifted_wald(1000, drift = 3, bs = 0.5, ndt = 0.2)
#' hist(rts, breaks = 50, main = "Simulated Shifted Wald RTs", xlab = "Reaction Time")
#'
#' @export
rshifted_wald <- function(n, drift = 3, bs = 0.5, ndt = 0.2) {
  # Prepare and validate all inputs for RNG
  params <- .prepare_shifted_wald(x = NULL, n = n, drift = drift, bs = bs, ndt = ndt)

  # Parameters for inverse Gaussian
  mu     <- params$bs / params$drift
  lambda <- params$bs^2

  # Generate IG draws via the two-root method
  y <- stats::rnorm(params$ndraws)^2
  z <- y * (mu / lambda)
  # Stable root ratio formula
  x1_over_mu <- 1 + z/2 * (1 - sqrt(1 + 4/z))

  # Root selection via uniform draws
  u <- stats::runif(params$ndraws)
  root_ratio <- ifelse(u < 1/(1 + x1_over_mu), x1_over_mu, 1/x1_over_mu)

  # Final draws: scale by mu and add shift
  ig_draws <- mu * root_ratio
  ig_draws + params$ndt
}



#' @rdname rshifted_wald
#' @param x Vector of quantiles (observed reaction times).
#' @param log Logical; if TRUE, probabilities p are given as log(p).
#' @export
dshifted_wald <- function(x, drift = 3, bs = 0.5, ndt = 0.2, log = FALSE) {
  # Prepare and validate inputs for density
  params <- .prepare_shifted_wald(x = x, n = NULL, drift = drift, bs = bs, ndt = ndt)

  # Time relative to non-decision component
  x_adj <- params$x - params$ndt
  log_density <- rep(NA_real_, params$ndraws)

  ## Valid indices where x > ndt
  valid <- x_adj > 0
  if (any(valid)) {
    mu    <- params$bs[valid] / params$drift[valid]
    lambda<- params$bs[valid]^2
    xa    <- x_adj[valid]

    # Standard inverse Gaussian log-PDF
    term1 <- 0.5 * (log(lambda) - log(2 * pi) - 3 * log(xa))
    term2 <- lambda * (xa - mu)^2 / (2 * mu^2 * xa)
    log_density[valid] <- term1 - term2
  }

  # For x <= ndt, density = 0 => log = -Inf
  log_density[!valid] <- -Inf

  if (log) {
    return(log_density)
  } else {
    return(exp(log_density))
  }
}


#' @rdname rshifted_wald
#' @param q Vector of quantiles (observed reaction times).
#' @param lower.tail Logical; if TRUE (default), probabilities are `P[X <= x]`, otherwise, `P[X > x]`.
#' @param log.p Logical; if TRUE, probabilities p are given as log(p). Defaults to FALSE.
#' @export
pshifted_wald <- function(q, drift = 3, bs = 0.5, ndt = 0.2, lower.tail = TRUE, log.p = FALSE) {
  params <- .prepare_shifted_wald(x = q, n = NULL, drift = drift, bs = bs, ndt = ndt)

  t    <- params$x - params$ndt
  cdf  <- rep(NA_real_, params$ndraws)

  # Cases where q <= ndt: CDF = 0
  cdf[t <= 0] <- 0

  # Cases where q = Inf: CDF = 1
  inf_idx <- is.infinite(t)
  cdf[inf_idx] <- 1

  # Remaining finite t > 0
  calc_idx <- (t > 0) & !inf_idx
  if (any(calc_idx)) {
    mu     <- params$bs[calc_idx] / params$drift[calc_idx]
    lambda <- params$bs[calc_idx]^2
    tv     <- t[calc_idx]

    sqrt_l_t <- sqrt(lambda / tv)
    t_mu     <- tv / mu

    term1 <- stats::pnorm(sqrt_l_t * (t_mu - 1))
    term2 <- exp(2 * lambda / mu) * stats::pnorm(-sqrt_l_t * (t_mu + 1))

    cdf_calc <- term1 + term2
    cdf_calc[cdf_calc < 0] <- 0
    cdf_calc[cdf_calc > 1] <- 1
    cdf[calc_idx] <- cdf_calc
  }

  if (!lower.tail) cdf <- 1 - cdf
  if (log.p) cdf <- log(cdf)
  cdf
}



# Internals ---------------------------------------------------------------




#' @keywords internal
.prepare_shifted_wald <- function(x = NULL, n = NULL, drift, bs, ndt) {
  # Validate parameters once
  if (any(drift <= 0))    stop("Drift rate 'drift' must be positive.")
  if (any(bs <= 0)) stop("Threshold 'bs' must be positive.")
  if (any(ndt < 0))    stop("Non-decision time 'ndt' must be non-negative.")

  # Determine target length:
  if (!is.null(x)) {
    # For PDF/CDF: based on x length
    m <- length(x)
  } else if (!is.null(n)) {
    # For RNG: based on n draw count
    if (length(n) != 1 || n < 0 || n != floor(n)) {
      stop("n must be a single non-negative integer.")
    }
    m <- n
  } else {
    stop("Either 'x' or 'n' must be provided.")
  }

  # Recycle vectors to length m
  params <- list(
    drift    = rep(drift,    length.out = m),
    bs = rep(bs, length.out = m),
    ndt   = rep(ndt,   length.out = m)
  )

  # Incorporate x if provided
  if (!is.null(x)) {
    params$x <- rep(x, length.out = m)
  }
  params$ndraws <- m
  params
}

