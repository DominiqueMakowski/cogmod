# Based on the ordbetareg package by Robert Kubinec (2023).

#' @title Beta-Extreme (BEXT) Model
#'
#' @description
#' The ordered beta model ([Kubinec, 2023](https://doi.org/10.1017/pan.2022.20)) was
#' introduced as an appropriate and parsimonious way of describing data commonly
#' observed in psychological science (such as from slider scales). It is defined
#' with a Beta distribution on the interval 0-1 with additional point masses at 0 and 1.
#'
#' The BeXt model corresponds to a reparametrized ordered beta model.
#' Instead of defining the left and right cutpoints, the BeXt parametrization
#' uses the likelihood of extreme values (0 and 1) and their balance (i.e., the
#' relative proportion of zeros and ones).
#'
#' It differs from the Zero-One-Inflated Beta (ZOIB) model in that the ZOIB model has `zoi`
#' and `coi` parameters, directly controlling the likelihood of extreme values. Instead,
#' BeXt uses `pex` and `bex` to define "cutpoints" after which extreme values become likely.
#'
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param mu Mean of the continuous component (between 0 and 1).
#' @param phi Precision parameter (positive).
#' @param pex Likelihood of extreme values (zeros and ones). Must be in the range 0-1.
#' @param bex Balance of extreme values. Represents the proportion of ones relative to zeros. Default: 0.5.
#'
#' @return A vector of simulated outcomes in the range 0-1.
#'
#' @examples
#' # Simulate data with different parameterizations
#' x <- rbext(10000, mu = 0.5, phi = 3, pex = 0.05, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#' @export
rbext <- function(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5) {
  # Validate inputs
  if (pex < 0 || pex > 1)
    stop("pex must be between 0 and 1")
  if (bex < 0 || bex > 1)
    stop("bex must be between 0 and 1")
  if (mu < 0 || mu > 1)
    stop("mu must be between 0 and 1")
  if (phi <= 0)
    stop("phi must be positive")

  # Compute kleft and kright from pex and bex
  kleft <- pex * (1 - bex)   # threshold below which outcomes are set to 0 (zeros)
  kright <- 1 - (pex * bex)  # threshold above which outcomes are set to 1 (ones)

  # Special case: pex = 1
  if (pex == 1) {
    # When pex = 1, all mass is at extremes, distributed according to bex
    return(sample(c(0, 1), size = n, replace = TRUE, prob = c(1-bex, bex)))
  }

  # Validate kleft and kright for cases where pex < 1
  if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1 || kleft >= kright)
    stop("Invalid parameterization: ensure 0 <= kleft < kright <= 1")

  # Draw uniform random numbers to determine the outcome category
  u <- stats::runif(n)

  outcomes <- numeric(n)
  outcomes[u < kleft] <- 0            # Category 1: outcome = 0
  outcomes[u >= kright] <- 1          # Category 3: outcome = 1

  # Category 2: continuous outcomes from the Beta distribution
  cont_idx <- which(u >= kleft & u < kright)
  if (length(cont_idx) > 0) {
    outcomes[cont_idx] <- stats::rbeta(length(cont_idx),
                                       shape1 = mu * phi,
                                       shape2 = (1 - mu) * phi)
  }

  outcomes
}

#' @rdname rbext
#' @param x A vector of values for which to compute the density.
#' @param log Logical; if TRUE, returns the log-density.
#' @examples
#' x <- seq(0, 1, length.out = 1001)
#' densities <- dbext(x, mu = 0.5, phi = 5, pex = 0.2, bex = 0.5)
#' plot(x, densities, type = "l", main = "Density Function", xlab = "y", ylab = "Density")
#' @export
dbext <- function(x, mu = 0.5, phi = 1, pex = 0.1, bex = 0.5, log = FALSE) {
  # Validate inputs
  if (!all(mu > 0 & mu < 1))
    stop("mu must be between 0 and 1")
  if (!all(phi > 0))
    stop("phi must be positive")
  if (!(length(mu) %in% c(1, length(x))))
    stop("Please pass a vector for mu that is either length 1 or the same length as x.")
  if (!(length(phi) %in% c(1, length(x))))
    stop("Please pass a vector for phi that is either length 1 or the same length as x.")
  if (pex < 0 || pex > 1)
    stop("pex must be between 0 and 1.")
  if (bex < 0 || bex > 1)
    stop("bex must be between 0 and 1.")

  # Initialize density vector
  dens <- numeric(length(x))

  # Special case: pex = 1 (discrete Bernoulli distribution)
  if (pex == 1) {
    # Handle values outside [0,1]
    dens[x < 0 | x > 1] <- 0

    # Handle exact zeros (point mass at 0)
    dens[x == 0] <- 1 - bex

    # Handle exact ones (point mass at 1)
    dens[x == 1] <- bex

    # All other values have zero density
    dens[x > 0 & x < 1] <- 0
  } else {
    # Regular case: mixture of continuous and discrete
    # Compute kleft and kright from pex and bex
    kleft <- pex * (1 - bex)   # probability mass at 0
    kright <- 1 - (pex * bex)  # threshold above which outcomes are set to 1

    # Validate parameters
    if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1 || kleft >= kright)
      stop("Invalid parameterization: ensure 0 <= kleft < kright <= 1")

    # Handle values outside [0,1]
    dens[x < 0 | x > 1] <- 0

    # Handle exact zeros (point mass at 0)
    dens[x == 0] <- kleft

    # Handle exact ones (point mass at 1)
    dens[x == 1] <- 1 - kright

    # Handle values in the continuous range (0,1)
    cont_idx <- which(x > 0 & x < 1)
    if (length(cont_idx) > 0) {
      # Get mu and phi for these indices
      mu_use <- if (length(mu) == 1) mu else mu[cont_idx]
      phi_use <- if (length(phi) == 1) phi else phi[cont_idx]

      # Calculate beta density scaled by the probability mass in the continuous component
      dens[cont_idx] <- (kright - kleft) * stats::dbeta(x[cont_idx],
                                                      shape1 = mu_use * phi_use,
                                                      shape2 = (1 - mu_use) * phi_use)
    }
  }

  # Return log density if requested
  if (log) {
    # Avoid log(0)
    log_dens <- rep(-Inf, length(dens))
    log_dens[dens > 0] <- log(dens[dens > 0])
    return(log_dens)
  } else {
    return(dens)
  }
}
