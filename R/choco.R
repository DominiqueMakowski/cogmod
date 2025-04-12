#' @title Choice-Confidence (CHOCO) Model Simulation
#'
#' @description
#' Simulates data from the Choice-Confidence (CHOCO) model. This model is useful for
#' subjective ratings (e.g., Likert-type scales) where responses represent a choice
#' between two underlying categories (e.g., "disagree" vs. "agree") along with a
#' degree of confidence or intensity.
#'
#' The CHOCO model divides the response scale at a `threshold`. Responses below the
#' threshold are modeled by a mirrored [Beta-Extreme (BEXT)][rbext()] distribution,
#' and responses above the threshold are modeled by a standard BEXT distribution.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param p Proportion parameter determining the balance between the left and right sides
#'   *after excluding* the probability mass at the threshold (`pmid`). Specifically,
#'   `P(Right Side | Not Threshold) = p` and `P(Left Side | Not Threshold) = 1 - p`.
#'   Must be in the range `[0, 1]`. In `brms` models, this corresponds to the `mu` parameter.
#' @param conf Mean parameter for the underlying BEXT distribution used for the *right* side
#'   (`threshold` to 1). This represents the central tendency (confidence) of the *raw* BEXT component
#'   (before scaling). Must be strictly between 0 and 1. The mean of the raw BEXT component
#'   for the *left* side (`muleft`) is derived from `conf` before applying `confleft`, such that
#'   `confleft = 0` implies mirrored confidence (`muleft = conf`).
#' @param confleft Difference parameter modifying the mean for the *left* side BEXT distribution.
#'   The mean of the raw BEXT component for the left side is calculated as
#'   `muleft = inv_logit(logit(conf) + confleft)`. A `confleft = 0` implies the raw mean for the
#'   left side is the same as the right side (`muleft = conf`). Positive values increase `muleft`,
#'   resulting in final CHOCO values closer to 0 (higher confidence on the left).
#' @param prec Precision parameter (`phi`) for the underlying BEXT distributions. This sets the
#'   base precision for the right side and is used to derive the left side's precision.
#'   Must be positive. Note: This corresponds to half the typical Beta precision
#'   (`precision = prec * 2`). `prec = 1` corresponds to a uniform distribution
#'   (when the respective raw mean is 0.5).
#' @param precleft Difference parameter modifying the precision for the *left* side BEXT distribution.
#'   The precision (`phi`) of the raw BEXT component for the left side is calculated as
#'   `prec_left = prec * exp(precleft)`. A `precleft = 0` implies the precisions are the same.
#'   Positive values increase the precision on the left side.
#' @param pex Overall probability of extreme values (0 or 1) *within* each underlying BEXT component
#'   (before scaling/mirroring). This applies commonly to both left and right sides. `0 <= pex <= 1`.
#' @param bex Balances the extreme probability mass *within* each underlying BEXT component
#'   between its 0 and 1 anchors. `P(raw=1) = pex * bex`, `P(raw=0) = pex * (1 - bex)`.
#'   Note that for the left CHOCO component, `raw=1` maps to CHOCO value 0, and for the right
#'   CHOCO component, `raw=1` maps to CHOCO value 1. `0 <= bex <= 1`. The interpretation
#'   of how `pex` and `bex` translate to the final CHOCO extremes (0 and 1) depends on the
#'   internal `pex_left` and `pex_right` calculations (see Details).
#' @param pmid Probability mass exactly at the `threshold`. `0 <= pmid <= 1`.
#' @param threshold The point dividing the scale into left and right components. Must be
#'   strictly between 0 and 1.
#'
#' @seealso rbext, dchoco
#'
#' @details
#' The simulation process involves three steps:
#' 1. Decide whether the outcome is exactly at the `threshold` (with probability `pmid`).
#' 2. If not at the threshold, decide whether the outcome falls on the right side (with
#'    conditional probability `p`) or the left side (with conditional probability `1-p`).
#' 3. Simulate the value from the corresponding BEXT distribution:
#'    - **Left Side (0 to `threshold`):** Calculate `mu_left = inv_logit(logit(conf) + confleft)` and
#'      `prec_left = prec * exp(precleft)`. The effective extreme probability for the underlying
#'      `rbext` call (targeting raw value 1, which maps to CHOCO 0) is calculated as
#'      `pex_left = pmin(1, pmax(0, (1 - bex) * (pex * 2)))`. Simulate `y_raw` from
#'      `rbext(mu=mu_left, phi=prec_left, pex=pex_left, bex=1)`. The final value is
#'      `(1 - y_raw) * threshold`.
#'    - **Right Side (`threshold` to 1):** Calculate `mu_right = conf` and `prec_right = prec`.
#'      The effective extreme probability for the underlying `rbext` call (targeting raw value 1,
#'      which maps to CHOCO 1) is calculated as `pex_right = pmin(1, pmax(0, bex * (pex * 2)))`.
#'      Simulate `y_raw` from `rbext(mu=mu_right, phi=prec_right, pex=pex_right, bex=1)`.
#'      The final value is `threshold + y_raw * (1 - threshold)`.
#'
#' The calculation `pex * 2` in `pex_left` and `pex_right` arises from the specific way `rbext`
#' handles its `pex` and `bex` parameters when simulating extremes. It ensures the intended
#' proportions land at CHOCO=0 and CHOCO=1 based on the input `pex` and `bex`.
#'
#' @references
#' - Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with
#'     lower and upper bounds. Political Analysis, 31(4), 519-536. (Describes the underlying ordered beta model)
#'
#' @examples
#' # Simulate data with different parameterizations
#' # 10% at threshold, 50/50 split otherwise, symmetric confidence/precision
#' x1 <- rchoco(n=5000, p = 0.5, conf = 0.5, confleft = 0, prec = 4,
#'   precleft = 0, pex = 0.1, bex = 0.5, pmid = 0, threshold = 0.5)
#' hist(x1, breaks = 50, main = "CHOCO: Symmetric Confidence", xlab = "y")
#'
#' # No threshold mass, 70% probability on right, higher confidence left
#' x2 <- rchoco(n=5000, p = 0.7, conf = 0.5, confleft = 1, prec = 3,
#'   precleft = 1, pex = 0.05, bex = 0.7, pmid = 0, threshold = 0.5)
#' hist(x2, breaks = 50, main = "CHOCO: Asymmetric p, Higher Conf Left", xlab = "y")
#'
#' # Lower confidence overall, high probability in the middle
#' x3 <- rchoco(n=5000, p = 0.5, conf = 0.2, confleft = 0, prec = 3,
#'   precleft = 0, pex = 0, bex = 0.5, pmid = 0.05, threshold = 0.5)
#' hist(x3, breaks = 50, main = "CHOCO: Low confidence overall", xlab = "y")
#' @rdname rchoco
#' @export
rchoco <- function(n,
                   p = 0.5,
                   conf = 0.5,
                   confleft = 0,
                   prec = 4,
                   precleft = 0,
                   pex = 0.1,
                   bex = 0.5,
                   pmid = 0,
                   threshold = 0.5) {

  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) stop("n must be a positive integer.")
  if (any(threshold <= 0 | threshold >= 1)) stop("threshold must be between 0 and 1 (exclusive).")
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1.")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1.")
  if (any(p < 0 | p > 1)) stop("p must be between 0 and 1.")
  if (any(pmid < 0 | pmid > 1)) stop("pmid must be between 0 and 1.")
  if (any(conf <= 0 | conf >= 1)) stop("conf must be between 0 and 1 (exclusive).")
  if (any(prec <= 0)) stop("prec must be positive.")

  # Tolerance for floating point comparisons
  eps <- 1e-9

  # --- Vectorization ---
  # Determine output length based on n and parameter vector lengths
  param_lengths <- c(length(p), length(conf), length(confleft), length(prec),
                     length(precleft), length(pex), length(bex), length(pmid), length(threshold))
  n_params <- max(param_lengths)
  n_out <- max(n, n_params)

  # Recycle parameters to match output length
  p <- rep(p, length.out = n_out)
  conf <- rep(conf, length.out = n_out)
  confleft <- rep(confleft, length.out = n_out)
  prec <- rep(prec, length.out = n_out)
  precleft <- rep(precleft, length.out = n_out)
  pex <- rep(pex, length.out = n_out)
  bex <- rep(bex, length.out = n_out)
  pmid <- rep(pmid, length.out = n_out)
  threshold <- rep(threshold, length.out = n_out)

  # --- Parameter computation ---
  # Right side parameters
  muright <- conf # Direct assignment
  phiright <- prec # Direct assignment

  # Left side parameters
  # Calculate muleft based on conf, adjusted by confleft
  # logit(x) = log(x / (1-x))
  # inv_logit(y) = exp(y) / (1 + exp(y))
  logit_conf <- log(conf / (1 - conf)) # Base logit-mean (mirrored confidence)
  logit_muleft <- logit_conf + confleft # Adjust logit-mean
  muleft <- exp(logit_muleft) / (1 + exp(logit_muleft)) # Transform back to mean
  phileft <- prec * exp(precleft) # Adjust precision

  # Clamp derived parameters to avoid issues at boundaries
  muright <- pmax(eps, pmin(muright, 1 - eps))
  muleft <- pmax(eps, pmin(muleft, 1 - eps))
  phiright <- pmax(eps, phiright) # Ensure phiright is positive
  phileft <- pmax(eps, phileft) # Ensure phileft is positive

  # Effective extreme probabilities for rbext calls (clamped)
  # pex_left corresponds to P(raw=1) for the left BEXT, which maps to CHOCO value 0.
  # pex_right corresponds to P(raw=1) for the right BEXT, which maps to CHOCO value 1.
  # The (1-bex) and bex terms distribute the overall pex probability mass to
  # the respective extremes (CHOCO=0 and CHOCO=1). The pex*2 scaling adjusts
  # for how rbext internally handles pex when bex=1 is specified.
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))

  # --- Simulation ---

  # 1. Decide for each trial: Left (0), Middle (0.5), or Right (1) component
  prob_left <- (1 - pmid) * (1 - p)
  prob_mid <- pmid
  prob_right <- (1 - pmid) * p

  # Sample side choice element-wise using uniform random numbers
  side_choice <- numeric(n_out)
  rand_unif <- stats::runif(n_out)
  side_choice[rand_unif < prob_left] <- 0 # Left component
  side_choice[rand_unif >= prob_left & rand_unif < (prob_left + prob_mid)] <- 0.5 # Middle (threshold)
  side_choice[rand_unif >= (prob_left + prob_mid)] <- 1 # Right component

  # Get indices for each component
  idx_left <- which(side_choice == 0)
  idx_mid <- which(side_choice == 0.5)
  idx_right <- which(side_choice == 1)
  n_left <- length(idx_left)
  n_mid <- length(idx_mid)
  n_right <- length(idx_right)

  # Initialize output vector
  out <- numeric(n_out)

  # 2. Simulate values at threshold
  if (n_mid > 0) {
      out[idx_mid] <- threshold[idx_mid] # Use vectorized threshold
  }

  # 3. Simulate Left Component (Mirrored and Rescaled)
  if (n_left > 0) {
      # Simulate underlying Beta-extreme values (mapped to [0,1])
      # Use bex=1 in rbext to force extremes only to raw=1 (which maps to CHOCO=0).
      # The effective pex for this call is pex_left.
      y_raw_left <- rbext(n = n_left, mu = muleft[idx_left], phi = phileft[idx_left],
                          pex = pex_left[idx_left], bex = 1)
      # Mirror (1 - y_raw) and scale by threshold
      out[idx_left] <- (1 - y_raw_left) * threshold[idx_left] # Use vectorized threshold
  }

  # 4. Simulate Right Component (Rescaled)
  if (n_right > 0) {
      # Simulate underlying Beta-extreme values (mapped to [0,1])
      # Use bex=1 in rbext to force extremes only to raw=1 (which maps to CHOCO=1).
      # The effective pex for this call is pex_right.
      y_raw_right <- rbext(n = n_right, mu = muright[idx_right], phi = phiright[idx_right],
                           pex = pex_right[idx_right], bex = 1)
      # Scale by (1 - threshold) and shift by threshold
      out[idx_right] <- threshold[idx_right] + y_raw_right * (1 - threshold[idx_right]) # Use vectorized threshold
  }

  # Return the combined results
  out
}


#' @rdname rchoco
#' @inheritParams rbext
#' @export
dchoco <- function(x, p = 0.5, conf = 0.5, confleft = 0, prec = 4, precleft = 0,
                   pex = 0.1, bex = 0.5, pmid = 0, threshold = 0.5, log = FALSE) {

  # --- Input Validation ---
  if (any(threshold <= 0 | threshold >= 1)) stop("threshold must be between 0 and 1 (exclusive).")
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1.")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1.")
  if (any(p < 0 | p > 1)) stop("p must be between 0 and 1.")
  if (any(pmid < 0 | pmid > 1)) stop("pmid must be between 0 and 1.")
  if (any(conf <= 0 | conf >= 1)) stop("conf must be between 0 and 1 (exclusive).")
  if (any(prec <= 0)) stop("prec must be positive.")

  # Tolerance for floating point comparisons
  eps <- 1e-9

  # --- Vectorization ---
  n <- length(x)
  p <- rep(p, length.out = n)
  conf <- rep(conf, length.out = n)
  confleft <- rep(confleft, length.out = n)
  prec <- rep(prec, length.out = n)
  precleft <- rep(precleft, length.out = n)
  pex <- rep(pex, length.out = n)
  bex <- rep(bex, length.out = n)
  pmid <- rep(pmid, length.out = n)
  threshold <- rep(threshold, length.out = n)

  # --- Parameter computation ---
  # Right side parameters
  muright <- conf
  phiright <- prec

  # Left side parameters
  logit_conf <- log(conf / (1 - conf))
  logit_muleft <- logit_conf + confleft
  muleft <- exp(logit_muleft) / (1 + exp(logit_muleft))
  phileft <- prec * exp(precleft)

  # Clamp derived parameters for numerical stability
  muright_clamped <- pmax(eps, pmin(muright, 1 - eps))
  muleft_clamped <- pmax(eps, pmin(muleft, 1 - eps))
  phiright_clamped <- pmax(eps, phiright)
  phileft_clamped <- pmax(eps, phileft)

  # Effective extreme probabilities (clamped)
  # These represent the probability mass assigned to the extremes (0 and 1)
  # *within* the respective rbext calls used in the simulation.
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2))) # Corresponds to CHOCO=0 extreme
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))      # Corresponds to CHOCO=1 extreme

  # --- Density Calculation ---
  density <- numeric(n)

  # Handle values outside [0, 1]
  outside_idx <- x < 0 | x > 1
  density[outside_idx] <- 0

  # Handle x = threshold (Point mass)
  thresh_idx <- abs(x - threshold) < eps & !outside_idx
  density[thresh_idx] <- pmid[thresh_idx]

  # Handle x = 0 (Point mass from Left BEXT's raw=1)
  # Probability = P(Left Component) * P(Extreme at raw=1 | Left Component)
  zero_idx <- abs(x - 0) < eps & !outside_idx & !thresh_idx
  density[zero_idx] <- (1 - pmid[zero_idx]) * (1 - p[zero_idx]) * pex_left[zero_idx]

  # Handle x = 1 (Point mass from Right BEXT's raw=1)
  # Probability = P(Right Component) * P(Extreme at raw=1 | Right Component)
  one_idx <- abs(x - 1) < eps & !outside_idx & !thresh_idx
  density[one_idx] <- (1 - pmid[one_idx]) * p[one_idx] * pex_right[one_idx]

  # Handle 0 < x < threshold (Continuous Left)
  # Density(x) = P(Left) * P(Continuous | Left) * dbeta(transformed_x) * |Jacobian|
  # P(Continuous | Left) = 1 - P(Extreme at raw=1 | Left) = 1 - pex_left
  # Transformation: x = (1 - raw_left) * threshold => raw_left = 1 - x / threshold
  # Jacobian: |d(raw_left)/dx| = 1 / threshold
  left_cont_idx <- x > eps & x < threshold - eps & !outside_idx & !thresh_idx & !zero_idx & !one_idx
  if (any(left_cont_idx)) {
    y_raw_left <- 1 - x[left_cont_idx] / threshold[left_cont_idx]
    # Clamp input to dbeta for stability, although less critical here than in simulation
    y_raw_left_clamped <- pmax(eps, pmin(y_raw_left, 1 - eps))
    beta_dens_left <- stats::dbeta(y_raw_left_clamped,
                                   shape1 = muleft_clamped[left_cont_idx] * phileft_clamped[left_cont_idx] * 2,
                                   shape2 = (1 - muleft_clamped[left_cont_idx]) * phileft_clamped[left_cont_idx] * 2)
    # Ensure density is 0 if beta_dens is NA/NaN (e.g., if parameters lead to invalid shapes)
    beta_dens_left[is.na(beta_dens_left)] <- 0
    density[left_cont_idx] <- (1 - pmid[left_cont_idx]) * (1 - p[left_cont_idx]) * (1 - pex_left[left_cont_idx]) *
                              beta_dens_left / threshold[left_cont_idx]
  }

  # Handle threshold < x < 1 (Continuous Right)
  # Density(x) = P(Right) * P(Continuous | Right) * dbeta(transformed_x) * |Jacobian|
  # P(Continuous | Right) = 1 - P(Extreme at raw=1 | Right) = 1 - pex_right
  # Transformation: x = threshold + raw_right * (1 - threshold) => raw_right = (x - threshold) / (1 - threshold)
  # Jacobian: |d(raw_right)/dx| = 1 / (1 - threshold)
  right_cont_idx <- x > threshold + eps & x < 1 - eps & !outside_idx & !thresh_idx & !zero_idx & !one_idx
  if (any(right_cont_idx)) {
    y_raw_right <- (x[right_cont_idx] - threshold[right_cont_idx]) / (1 - threshold[right_cont_idx])
    # Clamp input to dbeta for stability
    y_raw_right_clamped <- pmax(eps, pmin(y_raw_right, 1 - eps))
    beta_dens_right <- stats::dbeta(y_raw_right_clamped,
                                    shape1 = muright_clamped[right_cont_idx] * phiright_clamped[right_cont_idx] * 2,
                                    shape2 = (1 - muright_clamped[right_cont_idx]) * phiright_clamped[right_cont_idx] * 2)
    # Ensure density is 0 if beta_dens is NA/NaN
    beta_dens_right[is.na(beta_dens_right)] <- 0
    density[right_cont_idx] <- (1 - pmid[right_cont_idx]) * p[right_cont_idx] * (1 - pex_right[right_cont_idx]) *
                               beta_dens_right / (1 - threshold[right_cont_idx])
  }

  # Ensure density is non-negative
  density <- pmax(0, density)

  # Return log density if requested
  if (log) {
    # Avoid log(0) issues
    density[density == 0] <- -Inf
    density[density > 0] <- log(density[density > 0])
  } 
  density
}