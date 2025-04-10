# Based on the ordbetareg package by Robert Kubinec (2023).

#' @title Beta-Extreme (BEXT) Model
#'
#' @description
#' The BEXT model represents subjective ratings as a mixture of a continuous Beta distribution
#' with additional point masses at the extremes (0 and 1). This structure effectively captures
#' common patterns in subjective rating data where respondents often select extreme values
#' at higher rates than would be expected from a Beta distribution alone.
#'
#' The BeXt model corresponds to a reparametrized ordered beta model ([Kubinec, 2023](https://doi.org/10.1017/pan.2022.20)),
#' which was introduced as an appropriate and parsimonious way of describing data commonly
#' observed in psychological science (such as from slider scales). It is defined
#' with a Beta distribution on the interval 0-1 with additional point masses at 0 and 1.
#'
#' It differs from the Zero-One-Inflated Beta (ZOIB) model in that the ZOIB model has `zoi`
#' and `coi` parameters, directly controlling the likelihood of extreme values. Instead,
#' BeXt uses `pex` and `bex` to define "cutpoints" after which extreme values become likely.
#' In an ordered beta framework, the boundary probabilities arise through a single underlying
#' ordering process (the location of the cutpoints on the latent scale). In a ZOIB framework,
#' the boundaries are more like additional mass points inserted into a beta distribution.
#' In ordered beta models, the parameters represent a coherent underlying latent process.
#' Having cutpoints interact with mu better represents a single unified psychological process.
#' 
#' The BEXT model changes the parameterization of the ordered beta model to be more interpretable:
#' Instead of using fixed cutpoints, this version computes the cutpoints from two intuitive parameters:
#' - `pex` (p-extreme): the overall probability of extreme values (0 or 1).
#' - `bex` (balance-extreme): the balance of extreme probability mass between 0 and 1.
#' Because these parameters directly map onto the expected proportion of responses at the extremes, 
#' this approach provides greater flexibility (e.g., in capturing pronounced endpoints: It can adapt 
#' to various shapes - such as extreme clustering at zero or one) and intuitive interpretability.
#'
#' @param n Number of simulated values. Must be a positive integer.
#' @param mu Mean of the continuous Beta component (between 0 and 1), represents the central
#'   tendency of non-extreme values.
#' @param phi Precision parameter (positive). Note: This corresponds to half of the "typical"
#'   precision parameter used in Beta models (internally, `precision = phi * 2` is used).
#'   This reparametrization means `phi = 1` (with `mu = 0.5`) corresponds to a uniform distribution
#'   between 0 and 1 for the Beta component. It facilitates the usage of priors (e.g., on a log-link,
#'   `log(phi) = 0` corresponds to this uniform case).
#' @param pex Overall probability of extreme values (0 or 1). `0 <= pex <= 1`.
#' @param bex Balances the extreme probability mass between 0 and 1. If `bex = 0.5`, the mass
#'   is distributed equally (`P(0) = P(1) = pex / 2`). If `bex > 0.5`, more mass goes to 0.
#'   If `bex < 0.5`, more mass goes to 1. Specifically, `P(0) = pex * (1 - bex)` and `P(1) = pex * bex`.
#'
#' @details
#' The BeXt model corresponds to a reparametrized ordered beta model.
#' Instead of defining the left and right cutpoints directly, the BeXt parametrization
#' uses the overall probability of extreme values (`pex`) and their balance (`bex`).
#'
#' The probability masses at the extremes are determined by:
#' - Probability mass at 0: `P(0) = pex * (1 - bex)`
#' - Probability mass at 1: `P(1) = pex * bex`
#' The probability of a value falling between 0 and 1 (drawn from the Beta distribution) is `1 - pex`.
#'
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
#'   towards the upper anchor (1), while `bex < 0.5` suggests a bias towards the lower anchor (0).
#'
#' @seealso rchoco
#' 
#' @return A vector of simulated outcomes in the range 0-1.
#'
#' @references
#' - Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with
#'     lower and upper bounds. Political Analysis, 31(4), 519-536.
#'
#' @examples
#' # Simulate data with different parameterizations
#' x <- rbext(10000, mu = 0.5, phi = 2, pex = 0, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#' @export
rbext <- function(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5) {
  # Validate inputs (handle vector inputs)
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1")
  if (any(mu < 0 | mu > 1)) stop("mu must be between 0 and 1")
  if (any(phi <= 0)) stop("phi must be positive")

  # Handle vectorization if n > 1 and parameters are vectors
  # Determine the length needed based on n and parameter lengths
  n_out <- max(n, length(mu), length(phi), length(pex), length(bex))
  if (n_out > 1) {
      mu <- rep(mu, length.out = n_out)
      phi <- rep(phi, length.out = n_out)
      pex <- rep(pex, length.out = n_out)
      bex <- rep(bex, length.out = n_out)
  } else {
      n_out <- n # Ensure n_out is at least n
  }


  # Initialize output vector
  out <- numeric(n_out)
  precision <- phi * 2 # Vectorized calculation

  # --- Handle different cases using logical indexing ---

  # Case 1: pex == 1 (Pure Bernoulli)
  idx_pex1 <- which(abs(pex - 1) < 1e-9) # Use tolerance for float comparison
  if (length(idx_pex1) > 0) {
    out[idx_pex1] <- stats::rbinom(length(idx_pex1), size = 1, prob = bex[idx_pex1])
  }

  # Case 2: pex == 0 (Pure Beta)
  idx_pex0 <- which(abs(pex) < 1e-9 & !(1:n_out %in% idx_pex1)) # Exclude already processed
  if (length(idx_pex0) > 0) {
    out[idx_pex0] <- stats::rbeta(length(idx_pex0),
                                  mu[idx_pex0] * precision[idx_pex0],
                                  (1 - mu[idx_pex0]) * precision[idx_pex0])
  }

  # Case 3: 0 < pex < 1 (Mixture)
  idx_mix <- which(pex > 1e-9 & pex < (1 - 1e-9) & !(1:n_out %in% c(idx_pex1, idx_pex0)))
  if (length(idx_mix) > 0) {
    mu_mix <- mu[idx_mix]
    phi_mix <- phi[idx_mix]
    pex_mix <- pex[idx_mix]
    bex_mix <- bex[idx_mix]
    precision_mix <- precision[idx_mix]
    n_mix <- length(idx_mix)

    # Calculate cutpoints on the logit scale (vectorized)
    left_prob <- pex_mix * (1 - bex_mix)
    right_prob <- pex_mix * bex_mix

    # Avoid issues with qlogis(0) or qlogis(1)
    eps_q <- 1e-9
    left_prob <- pmax(eps_q, pmin(left_prob, 1 - eps_q))
    right_prob <- pmax(eps_q, pmin(right_prob, 1 - eps_q))
    mu_mix_adj <- pmax(eps_q, pmin(mu_mix, 1 - eps_q)) # Adjust mu for qlogis

    mu_logit <- stats::qlogis(mu_mix_adj)
    left_cutpoint <- mu_logit - stats::qlogis(1 - left_prob)
    right_cutpoint <- mu_logit - stats::qlogis(right_prob)

    # Generate latent values (vectorized)
    latent_values <- mu_logit + stats::rlogis(n_mix)

    # Determine category based on cutpoints (vectorized)
    out_mix <- numeric(n_mix)
    idx_mix_0 <- latent_values < left_cutpoint
    idx_mix_1 <- latent_values > right_cutpoint
    idx_mix_mid <- !(idx_mix_0 | idx_mix_1)

    out_mix[idx_mix_0] <- 0
    out_mix[idx_mix_1] <- 1

    # Generate Beta values for middle category (vectorized subset)
    n_mix_mid <- sum(idx_mix_mid)
    if (n_mix_mid > 0) {
      out_mix[idx_mix_mid] <- stats::rbeta(
        n_mix_mid,
        shape1 = mu_mix[idx_mix_mid] * precision_mix[idx_mix_mid],
        shape2 = (1 - mu_mix[idx_mix_mid]) * precision_mix[idx_mix_mid]
      )
    }
    out[idx_mix] <- out_mix
  }

  # Ensure output length matches original n if n was the max dimension
  if (n_out == n && n > 0) {
      return(out[1:n])
  } else if (n_out > n && n == 1) {
      # If n=1 but parameters were vectors, return the vector
      return(out)
  } else {
      # Default case or if parameter vectors determined n_out > n
      return(out)
  }
}

#' @rdname rbext
#' @param x Vector of quantiles (values at which to evaluate the density). Must be between 0 and 1, inclusive.
#' @param log Logical; if TRUE, returns the log-density.
#' @examples
#' x <- seq(0, 1, length.out = 1001)
#' densities <- dbext(x, mu = 0.5, phi = 5, pex = 0.2, bex = 0.5)
#' plot(x, densities, type = "l", main = "Density Function", xlab = "y", ylab = "Density")
#' @export
dbext <- function(x, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5, log = FALSE) {
  # --- Input Validation ---
  if (any(mu < 0 | mu > 1)) {
      warning("mu must be between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(phi <= 0)) {
      warning("phi must be positive. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(pex < 0 | pex > 1)) {
      warning("pex must be between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(bex < 0 | bex > 1)) {
      warning("bex must be between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }

  # Recycle parameters to match length of x
  n <- length(x)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  pex <- rep(pex, length.out = n)
  bex <- rep(bex, length.out = n)

  # Calculate internal precision
  precision <- phi * 2

  # Initialize density vector
  density <- numeric(n)

  # Handle values outside [0, 1] - return 0/-Inf without warning
  outside_idx <- x < 0 | x > 1
  density[outside_idx] <- 0

  # Handle x = 0
  zero_idx <- x == 0 & !outside_idx
  density[zero_idx] <- pex[zero_idx] * (1 - bex[zero_idx])

  # Handle x = 1
  one_idx <- x == 1 & !outside_idx
  density[one_idx] <- pex[one_idx] * bex[one_idx]

  # Handle 0 < x < 1
  middle_idx <- x > 0 & x < 1 & !outside_idx
  if (any(middle_idx)) {
    # Handle edge case pex = 1 for middle values (density should be 0)
    pex_one_middle_idx <- middle_idx & (pex == 1)
    density[pex_one_middle_idx] <- 0

    # Calculate density for valid middle values where pex < 1
    valid_middle_idx <- middle_idx & (pex < 1)
    if(any(valid_middle_idx)) {
        density[valid_middle_idx] <- (1 - pex[valid_middle_idx]) * stats::dbeta(
        x[valid_middle_idx],
        shape1 = mu[valid_middle_idx] * precision[valid_middle_idx],
        shape2 = (1 - mu[valid_middle_idx]) * precision[valid_middle_idx]
      )
    }
  }

  # Return log density if requested
  if (log) {
    # Avoid log(0) warnings; replace 0 density with -Inf
    log_density <- ifelse(density == 0, -Inf, log(density))
    # Ensure any NaN from dbeta (e.g., if mu=0/1 and precision is low) becomes -Inf
    log_density[is.na(log_density)] <- -Inf
    return(log_density)
  } else {
    # Ensure any NaN from dbeta becomes 0
    density[is.na(density)] <- 0
    return(density)
  }
}

