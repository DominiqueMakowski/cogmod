#' @title Choice-Confidence (CHOCO) Model
#'
#' @description
#' The Choice-Confidence (CHOCO) model is useful to model data from subjective ratings,
#' such as Likert-type or analog scales, in which the left and the right side correspond
#' to different processes or higher order categorical responses (e.g., "disagree" vs.
#' "agree", "true" vs. "false"). They can be used to jointly model choice (left or right)
#' and confidence (the degree of left or right).
#'
#' To better represent potentially bimodal bounded data with zeros and ones, CHOCO models
#' consist of two separate ordered beta distributions (Kubinec, 2023) reparametrized as
#' BetaXtreme to model the left and the right hand side of the scale.
#'
#' The left-hand side is defined by `muleft` and `phileft`, while the right-hand side is defined
#' as the mirrored value of `muleft` around the threshold. A general parameter `p` controls the
#' proportion of data on the right-hand side versus the left-hand side of the threshold.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param p General proportion of data on the right-hand side of the threshold. Must be in the range `[0, 1]`.
#'   If `p = 0.3`, 30% of the data will be on the right-hand side.
#' @param muleft Mean of the Beta distribution for the left-hand side (range `[0, threshold]`).
#'   Higher values push the mass toward `threshold`.
#' @param mudelta Deviation in mean for the right-hand side relative to the left-hand side.
#' @param phileft Precision parameter of the Beta distribution for the left-hand side (positive).
#' @param phidelta Deviation in precision for the right-hand side relative to the left-hand side,
#'   expressed on the log scale.
#' @param pex Likelihood of extreme values (zeros and ones). Must be in the range `[0, 1]`.
#' @param bex Balance of extreme values. Represents the proportion of ones relative to zeros. Default: 0.5.
#' @param threshold The threshold separating the two Beta distributions. Default: 0.5.
#'
#' @return A vector of simulated outcomes in the range `[0, 1]`.
#'
#' @examples
#' # Simulate data with different parameterizations
#' x <- rchoco(50000, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5,
#'   phidelta = 0, pex = 0.1, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#'
#' x <- rchoco(50000, p = 0.7, muleft = 0.3, mudelta = 0, phileft = 5,
#'   phidelta = 0, pex = 0.1, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#' @export
rchoco <- function(n, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 4, phidelta = 0, pex = 0.1, bex = 0.5, threshold = 0.5) {
  # Validate inputs
  if (p < 0 || p > 1)
    stop("p must be between 0 and 1")
  if (pex < 0 || pex > 1)
    stop("pex must be between 0 and 1")
  if (bex < 0 || bex > 1)
    stop("bex must be between 0 and 1")
  if (muleft <= 0 || muleft >= 1)
    stop("muleft must be between 0 and 1 (exclusive)")
  if (phileft <= 0)
    stop("phileft must be positive")
  if (threshold <= 0 || threshold >= 1)
    stop("threshold must be between 0 and 1")

  # Special case: pex = 1 (all mass at extremes)
  if (pex == 1) {
    # When pex = 1, all outcomes are either 0 or 1, distributed according to p and bex
    # p controls overall left/right proportion, bex controls extreme value balance
    prob0 <- (1 - p) * (1 - bex) / ((1 - p) * (1 - bex) + p * bex)
    return(sample(c(0, 1), size = n, replace = TRUE, prob = c(prob0, 1 - prob0)))
  }

  # Transform parameters for the right-side distribution
  # Helper function to compute right-side parameters from left-side
  compute_right_params <- function(muleft, mudelta, phileft, phidelta) {
    # Mirror muleft in logit space and add mudelta
    logit_muleft <- log(muleft / (1 - muleft))
    logit_muright <- -logit_muleft + mudelta
    muright <- exp(logit_muright) / (1 + exp(logit_muright))

    # Adjust precision parameter
    phiright <- phileft * exp(phidelta)

    return(list(muright = muright, phiright = phiright))
  }

  # Get right-side parameters
  right_params <- compute_right_params(muleft, mudelta, phileft, phidelta)
  muright <- right_params$muright
  phiright <- right_params$phiright

  # Validate derived parameters
  if (muright <= 0 || muright >= 1)
    stop("Invalid parameterization: muright = ", round(muright, 4),
         ". Try adjusting muleft and mudelta.")

  # Compute probabilities for zeros and ones
  zero_prob <- pex * (1 - bex) * (1 - p)  # Probability of observing a zero
  one_prob <- pex * bex * p               # Probability of observing a one

  # Validate extreme value probabilities
  if (zero_prob + one_prob > 1)
    stop("Invalid parameterization: sum of zero and one probabilities exceeds 1. ",
         "Try reducing pex, or adjusting p and bex.")

  # Initialize outcomes
  outcomes <- numeric(n)

  # Draw uniform random numbers for classification
  u <- stats::runif(n)

  # Assign values based on the drawn uniform numbers
  zero_idx <- u < zero_prob
  one_idx <- u >= (1 - one_prob)
  cont_idx <- !(zero_idx | one_idx)

  # Assign extreme values
  outcomes[zero_idx] <- 0
  outcomes[one_idx] <- 1

  # Process continuous outcomes
  if (any(cont_idx)) {
    # Rescale uniform values to determine left/right assignment
    # (more efficient than drawing new random numbers)
    scaled_u <- (u[cont_idx] - zero_prob) / (1 - zero_prob - one_prob)
    left_side <- scaled_u < (1 - p)

    # Generate Beta variates for left side
    if (any(left_side)) {
      outcomes[cont_idx][left_side] <- threshold * stats::rbeta(
        sum(left_side),
        shape1 = muleft * phileft,
        shape2 = (1 - muleft) * phileft
      )
    }

    # Generate Beta variates for right side
    if (any(!left_side)) {
      outcomes[cont_idx][!left_side] <- threshold + (1 - threshold) * stats::rbeta(
        sum(!left_side),
        shape1 = muright * phiright,
        shape2 = (1 - muright) * phiright
      )
    }
  }

  outcomes
}


#' @rdname rchoco
#' @examples
#' x <- seq(0, 1, length.out = 101)
#' dens <- dchoco(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0, pex = 0.1, bex = 0.5)
#' plot(x, dens, type = "l", main = "Density of CHOCO Model", xlab = "x", ylab = "Density")
#' @export
dchoco <- function(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 4, phidelta = 0,
                   pex = 0.1, bex = 0.5, threshold = 0.5, log = FALSE) {
  # Validate inputs
  if (p < 0 || p > 1)
    stop("p must be between 0 and 1")
  if (pex < 0 || pex > 1)
    stop("pex must be between 0 and 1")
  if (bex < 0 || bex > 1)
    stop("bex must be between 0 and 1")
  if (muleft <= 0 || muleft >= 1)
    stop("muleft must be between 0 and 1 (exclusive)")
  if (phileft <= 0)
    stop("phileft must be positive")
  if (threshold <= 0 || threshold >= 1)
    stop("threshold must be between 0 and 1")
    
  # Transform parameters for the right-side distribution
  # Helper function to compute right-side parameters from left-side
  compute_right_params <- function(muleft, mudelta, phileft, phidelta) {
    # Mirror muleft in logit space and add mudelta
    logit_muleft <- log(muleft / (1 - muleft))
    logit_muright <- -logit_muleft + mudelta
    muright <- exp(logit_muright) / (1 + exp(logit_muright))
    
    # Adjust precision parameter
    phiright <- phileft * exp(phidelta)
    
    return(list(muright = muright, phiright = phiright))
  }
  
  # Get right-side parameters
  right_params <- compute_right_params(muleft, mudelta, phileft, phidelta)
  muright <- right_params$muright
  phiright <- right_params$phiright
  
  # Validate derived parameters
  if (muright <= 0 || muright >= 1)
    stop("Invalid parameterization: muright = ", round(muright, 4),
         ". Try adjusting muleft and mudelta.")
  
  # Initialize density vector
  dens <- numeric(length(x))
  
  # Special case: pex = 1 (all mass at extremes)
  if (pex == 1) {
    # When pex = 1, all outcomes are either 0 or 1, distributed according to p and bex
    prob0 <- (1 - p) * (1 - bex) / ((1 - p) * (1 - bex) + p * bex)
    prob1 <- 1 - prob0
    
    # Handle values outside [0,1]
    dens[x < 0 | x > 1] <- 0
    
    # Handle exact zeros and ones
    dens[x == 0] <- prob0
    dens[x == 1] <- prob1
    
    # All other values have zero density
    dens[x > 0 & x < 1] <- 0
  } else {
    # Compute probabilities for zeros and ones
    zero_prob <- pex * (1 - bex) * (1 - p)  # Probability of observing a zero
    one_prob <- pex * bex * p               # Probability of observing a one
    
    # Validate extreme value probabilities
    if (zero_prob + one_prob > 1)
      stop("Invalid parameterization: sum of zero and one probabilities exceeds 1. ",
           "Try reducing pex, or adjusting p and bex.")
    
    # Handle values outside [0,1]
    dens[x < 0 | x > 1] <- 0
    
    # Handle exact zeros (point mass at 0)
    dens[x == 0] <- zero_prob
    
    # Handle exact ones (point mass at 1)
    dens[x == 1] <- one_prob
    
    # Handle values in the continuous range
    cont_mass <- 1 - zero_prob - one_prob

    # Scale Beta density by p
    left_scale <- cont_mass * (1 - p)
    right_scale <- cont_mass * p

    # Left side of threshold (only if there's any mass on the left)
    if (p < 1) {
      left_idx <- which(x > 0 & x < threshold)
      if (length(left_idx) > 0) {
        # Scale x to [0,1] range for Beta
        scaled_x <- x[left_idx] / threshold
        
        # Calculate scaled beta density
        dens[left_idx] <- left_scale * stats::dbeta(
          scaled_x,
          shape1 = muleft * phileft,
          shape2 = (1 - muleft) * phileft
        ) / threshold  # Additional scaling factor due to change of variables
      }
    }
    
    # Right side of threshold (only if there's any mass on the right)
    if (p > 0) {
      right_idx <- which(x > threshold & x < 1)
      if (length(right_idx) > 0) {
        # Scale x to [0,1] range for Beta
        scaled_x <- (x[right_idx] - threshold) / (1 - threshold)
        
        # Calculate scaled beta density
        dens[right_idx] <- right_scale * stats::dbeta(
          scaled_x,
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright
        ) / (1 - threshold)  # Additional scaling factor due to change of variables
      }
    }
    
    # Special handling for exactly at threshold
    thresh_idx <- which(x == threshold)
    if (length(thresh_idx) > 0) {
      # Use points very close to boundaries instead of exactly at them
      overlap <- 1e-6
      
      if (p == 0) {
        # When p = 0, only use left limit
        left_limit <- left_scale * stats::dbeta(
          1.0 - overlap,  # Almost at x=1 but not quite
          shape1 = muleft * phileft,
          shape2 = (1 - muleft) * phileft
        ) / threshold
        dens[thresh_idx] <- left_limit
      } else if (p == 1) {
        # When p = 1, only use right limit
        right_limit <- right_scale * stats::dbeta(
          0.0 + overlap,  # Almost at x=0 but not quite
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright
        ) / (1 - threshold)
        dens[thresh_idx] <- right_limit
      } else {
        # Normal case: use weighted average
        left_limit <- left_scale * stats::dbeta(
          1.0 - overlap,  # Almost at x=1 but not quite
          shape1 = muleft * phileft,
          shape2 = (1 - muleft) * phileft
        ) / threshold
        
        right_limit <- right_scale * stats::dbeta(
          0.0 + overlap,  # Almost at x=0 but not quite
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright
        ) / (1 - threshold)
        
        dens[thresh_idx] <- (1 - p) * left_limit + p * right_limit
      }
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