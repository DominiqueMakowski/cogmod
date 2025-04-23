#' @title Beta-Gate Model
#'
#' @description
#' The Beta-Gate model represents subjective ratings as a mixture of a continuous Beta distribution
#' with additional point masses at the extremes (0 and 1). This structure effectively captures
#' common patterns in subjective rating data where respondents often select extreme values
#' at higher rates than would be expected from a Beta distribution alone.
#'
#' The Beta-Gate model corresponds to a reparametrized ordered beta model ([Kubinec, 2023](https://doi.org/10.1017/pan.2022.20)).
#' In the ordered Beta model, the extreme values (0 and 1) arise from censoring an underlying
#' latent process based on cutpoints ("gates"). Values falling past the gates are considered extremes
#' (zeros and ones). The difference from the Ordered Beta is the way the cutpoints are defined,
#' as well as the scale of the precision parameter phi.
#'
#' It differs from the Zero-One-Inflated Beta (ZOIB) model in that the ZOIB model has `zoi`
#' and `coi` parameters, directly controlling the likelihood of extreme values. Instead,
#' Beta-Gate uses `pex` and `bex` to define "cutpoints" after which extreme values become likely.
#' In an ordered beta framework, the boundary probabilities arise through a single underlying
#' ordering process (the location of the cutpoints on the latent scale). In a ZOIB framework,
#' the boundaries are more like additional mass points inserted into a beta distribution.
#' In Beta-gate models, extreme values arise naturally from thresholding a single latent process.
#'
#' @param n Number of simulated values.
#' @param mu Mean of the underlying Beta distribution (0 < mu < 1).
#' @param phi Precision parameter of the underlying Beta distribution (phi > 0).
#'   Note: `precision = phi * 2`. `phi = 1` corresponds to uniform when `mu = 0.5`.
#' @param pex Controls the location of the lower and upper boundary gates (`0 <= pex <= 1`). It defines 
#'   the total probability mass allocated to the extremes (0 or 1).  Higher `pex` increases the probability 
#'   of extreme values (0 or 1).
#' @param bex Balances the extreme probability mass `pex` between 0 and 1 (`0 <= bex <= 1`). A balance
#'   of `0.5` means that the 'gates' are symmetrically placed around the center of the distribution, and 
#'   values higher or lower than `0.5` will shift the relative "ease" of crossing the gates towards 1 
#'   or 0, respectively.
#'
#' @details
#' **Special cases:**
#' - When `pex = 0`: Pure Beta distribution with mean `mu` and precision `phi * 2`.
#' - When `pex = 1`: Pure Bernoulli distribution with `P(1) = bex`, `P(0) = 1-bex`.
#' - When `bex = 0` and `pex = 1`: All mass at 0.
#' - When `bex = 1` and `pex = 1`: All mass at 1.
#'
#' **Psychological Interpretation:**
#' - `mu`: Can be interpreted as the underlying average tendency or preference strength,
#'   disregarding extreme "all-or-nothing" responses.
#' - `phi`: Reflects the certainty or consistency of the non-extreme responses. Higher `phi`
#'   indicates responses tightly clustered around `mu` (more certainty), while lower `phi`
#'   (especially `phi = 1`) suggests more uniform or uncertain responses.
#' - `pex`: Represents the overall tendency towards extreme responding (choosing 0 or 1).
#'   This could reflect individual response styles (e.g., acquiescence, yea-saying/nay-saying)
#'   or properties of the item itself (e.g., polarizing questions).
#' - `bex`: Indicates the *direction* of the extreme response bias. `bex > 0.5` suggests a bias
#'   for producing ones more easily, while `bex < 0.5` suggests a bias towards zero.
#'
#' @return A vector of simulated outcomes in the range 0-1.
#'
#' @references
#' - Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with
#'     lower and upper bounds. Political Analysis, 31(4), 519-536.
#'
#' @examples
#' # Symmetric gates (c0=0.05, c1=0.95), pex=0.1, bex=0.5
#' x1 <- rbetagate(10000, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5)
#' # hist(x1, breaks=50, main="rbetagate: Symmetric Cutpoints (pex=0.1)")
#'
#' # Asymmetric gates (c0=0.15, c1=0.95), pex=0.2, bex=0.25
#' x2 <- rbetagate(10000, mu = 0.5, phi = 3, pex = 0.2, bex = 0.25)
#' # hist(x2, breaks=50, main="rbetagate: Asymmetric Cutpoints (pex=0.2, bex=0.25)")
#'
#' # No gating (pure Beta)
#' x3 <- rbetagate(10000, mu = 0.7, phi = 5, pex = 0, bex = 0.5)
#' # hist(x3, breaks=50, main="rbetagate: No Extreme Values (pex=0)")
#'
#' @export
rbetagate <- function(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5) {
  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) stop("n must be a positive integer.")
  if (any(mu <= 0 | mu >= 1)) stop("mu must be strictly between 0 and 1.")
  if (any(phi <= 0)) stop("phi must be positive.")
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1.")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1.")

  # --- Vectorization ---
  n_out <- max(n, length(mu), length(phi), length(pex), length(bex))
  if (n_out > 1 || n > 1) { # Ensure vectorization if n>1 even if params are scalar
      mu <- rep(mu, length.out = n_out)
      phi <- rep(phi, length.out = n_out)
      pex <- rep(pex, length.out = n_out)
      bex <- rep(bex, length.out = n_out)
  } else {
      n_out <- n # Case where n=1 and params are scalar
  }

  # --- Parameter Calculation ---
  eps <- .Machine$double.eps # Smallest representable positive number

  # Cutpoints on probability scale
  cutzero <- pex * (1 - bex)
  cutone <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero)
  cutonelog <- stats::qlogis(cutone)

  # Beta distribution parameters
  shape1 <- mu * phi * 2
  shape2 <- (1 - mu) * phi * 2

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # --- Probabilities for Outcome Categories ---
  # P(outcome = 0) = P(eta < cutzerolog) = P(logistic(mu_ql) < cutzero)
  # P(outcome = 1) = P(eta > cutonelog) = P(logistic(mu_ql) > cutone)
  # P(0 < outcome < 1) = P(cutzerolog <= eta <= cutonelog)
  # Using plogis(q, location) = P(X <= q) where X ~ Logistic(location, scale=1)
  # P(eta < cutzerolog) = plogis(cutzerolog, location = mu_ql) = 1 - plogis(mu_ql - cutzerolog)
  # P(eta > cutonelog) = 1 - plogis(cutonelog, location = mu_ql) = plogis(mu_ql - cutonelog)
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE) # P(eta <= cutzerolog)
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE) # P(eta > cutonelog)
  prob_mid <- 1 - prob_0 - prob_1
  # Ensure probabilities sum to 1 and handle potential floating point inaccuracies
  prob_mid <- pmax(0, prob_mid)
  probs_matrix <- cbind(prob_0, prob_mid, prob_1)
  probs_matrix <- probs_matrix / rowSums(probs_matrix) # Normalize row-wise

  # --- Simulation ---
  # Sample outcome category (0, middle, 1) for each observation
  # Uses Gumbel-max trick for efficient multinomial sampling
  gumbel_noise <- -log(-log(matrix(stats::runif(n_out * 3), nrow = n_out, ncol = 3)))
  outcome_category <- max.col(log(probs_matrix) + gumbel_noise) # 1 for 0, 2 for middle, 3 for 1

  # Generate underlying Beta draws for all (only used for middle category)
  beta_draws <- stats::rbeta(n = n_out, shape1 = shape1, shape2 = shape2)
  # Clamp beta draws slightly away from exact 0/1 for stability if needed downstream
  beta_draws <- pmax(eps, pmin(1 - eps, beta_draws))

  # --- Combine Results ---
  final_out <- numeric(n_out)
  final_out[outcome_category == 1] <- 0
  final_out[outcome_category == 3] <- 1
  mid_indices <- which(outcome_category == 2)
  if (length(mid_indices) > 0) {
      final_out[mid_indices] <- beta_draws[mid_indices]
  }

  # Ensure output length matches original n if n was the max dimension
  if (n_out == n && n > 0) {
      return(final_out[1:n])
  } else if (n_out > n && n == 1) {
      return(final_out) # Return vector if n=1 but params were vectors
  } else {
      return(final_out) # Default case
  }
}


#' @rdname rbetagate
#' @param x Vector of quantiles (values at which to evaluate the density). Must be between 0 and 1, inclusive.
#' @param log Logical; if TRUE, returns the log-density.
#' @examples
#' x <- seq(0, 1, length.out = 1001)
#' densities <- dbetagate(x, mu = 0.5, phi = 5, pex = 0.2, bex = 0.5)
#' plot(x, densities, type = "l", main = "Density Function", xlab = "y", ylab = "Density")
#' @export
dbetagate <- function(x, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5, log = FALSE) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) stop("mu must be strictly between 0 and 1.")
  if (any(phi <= 0)) stop("phi must be positive.")
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1.")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1.")
  # Allow x=0 and x=1
  if (any(x < 0 | x > 1)) warning("x must be between 0 and 1.")

  # --- Vectorization ---
  n <- length(x)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  pex <- rep(pex, length.out = n)
  bex <- rep(bex, length.out = n)

  # --- Parameter Calculation ---
  eps <- .Machine$double.eps # Smallest representable positive number

  # Cutpoints on probability scale
  cutzero <- pex * (1 - bex)
  cutone <- 1 - pex * bex

  # Cutpoints on logit scale
  cutzerolog <- stats::qlogis(cutzero)
  cutonelog <- stats::qlogis(cutone)

  # Beta distribution parameters
  shape1 <- mu * phi * 2
  shape2 <- (1 - mu) * phi * 2

  # Location parameter on logit scale
  mu_ql <- stats::qlogis(mu)

  # --- Probabilities for Outcome Categories ---
  # See rbetagate for derivation
  prob_0 <- stats::plogis(cutzerolog, location = mu_ql, lower.tail = TRUE)
  prob_1 <- stats::plogis(cutonelog, location = mu_ql, lower.tail = FALSE)
  prob_mid <- 1 - prob_0 - prob_1
  prob_mid <- pmax(0, prob_mid) # Handle potential floating point inaccuracies

  # --- Calculate Density ---
  # Initialize density vector
  density <- numeric(n)

  # Indices for different cases
  idx_zero <- which(x == 0)
  idx_one <- which(x == 1)
  idx_mid <- which(x > 0 & x < 1)
  idx_outside <- which(x < 0 | x > 1) # Should have 0 density

  # Density for x = 0
  if (length(idx_zero) > 0) {
    density[idx_zero] <- prob_0[idx_zero]
  }

  # Density for x = 1
  if (length(idx_one) > 0) {
    density[idx_one] <- prob_1[idx_one]
  }

  # Density for 0 < x < 1
  if (length(idx_mid) > 0) {
    # Calculate Beta density component
    # Need to handle case where prob_mid is zero (e.g., pex=1)
    beta_dens <- ifelse(prob_mid[idx_mid] > eps,
                        stats::dbeta(x[idx_mid], shape1 = shape1[idx_mid], shape2 = shape2[idx_mid], log = FALSE),
                        0)
    # Total density is P(middle category) * BetaPDF(x | params)
    density[idx_mid] <- prob_mid[idx_mid] * beta_dens
  }

  # Density for x outside [0, 1] is 0 (already initialized)

  # --- Return Log Density if Requested ---
  if (log) {
    # Use log1p for potentially small probabilities near 0? No, log(prob) is fine.
    # Handle density = 0 case -> log(0) = -Inf
    density <- ifelse(density > 0, log(density), -Inf)
  }

  density
}
