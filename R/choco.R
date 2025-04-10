#' @title Choice-Confidence (CHOCO) Model
#'
#' @description
#' The Choice-Confidence (CHOCO) model is useful to model data from subjective ratings,
#' such as Likert-type or analog scales, in which the left and the right side correspond
#' to different processes or higher order categorical responses (e.g., "disagree" vs.
#' "agree", "true" vs. "false"). They can be used to jointly model the choice
#' between two underlying categories (e.g., "disagree" vs. "agree") along with a degree
#' of confidence or intensity.
#'
#' The CHOCO model conceptualizes the response scale as being divided at a `threshold` (typically 
#' at 0.5, i.e., the middle point of the scale). Responses below the threshold are modeled by one 
#' [Beta-Extreme (BEXT)][rbext()] distribution, and responses above the threshold are modeled 
#' by a separate, mirrored BEXT distribution. The relative proportions of data on the left and right
#' sides of the threshold are controlled by the `p` parameter, which indicates the probability of 
#' observing a response on the right side of the threshold. There is also a point mass probability 
#' (`pmid`) for responses exactly at the `threshold`.
#' 
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param p Proportion parameter determining the balance between the left and right sides
#'   *after excluding* the probability mass at the threshold (`pmid`). Specifically,
#'   `P(Right Side | Not Threshold) = p` and `P(Left Side | Not Threshold) = 1 - p`.
#'   Must be in the range `[0, 1]`. Due to Stan requirement, this parameter is named `mu` in
#'   brms models and is its primary parameter.
#' @param muleft Mean parameter for the underlying BEXT distribution used for the *left* side
#'   (0 to `threshold`). This represents the central tendency of the *raw* BEXT component
#'   (before mirroring and scaling). Must be strictly between 0 and 1.
#' @param mudelta Difference parameter modifying the mean for the *right* side BEXT distribution.
#'   The mean of the raw BEXT component for the right side is calculated as
#'   `inv_logit(logit(muleft) + mudelta)`. An `mudelta = 0` implies the raw means are the same.
#' @param phileft Precision parameter for the underlying BEXT distribution on the *left* side.
#'   Must be positive. Note: This corresponds to half the typical Beta precision (`precision = phileft * 2`).
#'   `phileft = 1` corresponds to a uniform distribution (when `muleft = 0.5`).
#' @param phidelta Difference parameter modifying the precision for the *right* side BEXT distribution.
#'   The precision of the raw BEXT component for the right side is calculated as
#'   `phiright = phileft * exp(phidelta)`. A `phidelta = 0` implies the precisions are the same.
#' @param pex Overall probability of extreme values (0 or 1) *within* each underlying BEXT component
#'   (before scaling/mirroring). This applies commonly to both left and right sides. `0 <= pex <= 1`.
#' @param bex Balances the extreme probability mass *within* each underlying BEXT component
#'   between its 0 and 1 anchors. `P(raw=1) = pex * bex`, `P(raw=0) = pex * (1 - bex)`.
#'   Note that for the left CHOCO component, `raw=1` maps to CHOCO value 0, and for the right
#'   CHOCO component, `raw=1` maps to CHOCO value 1. `0 <= bex <= 1`.
#' @param pmid Probability mass exactly at the `threshold`. `0 <= pmid <= 1`.
#' @param threshold The point dividing the scale into left and right components. Must be
#'   strictly between 0 and 1.
#' 
#' @seealso rbext
#'
#' @details
#' The simulation process involves three steps:
#' 1. Decide whether the outcome is exactly at the `threshold` (with probability `pmid`).
#' 2. If not at the threshold, decide whether the outcome falls on the right side (with
#'    conditional probability `p`) or the left side (with conditional probability `1-p`).
#' 3. Simulate the value from the corresponding BEXT distribution:
#'    - **Left Side (0 to `threshold`):** Simulate `y_raw` from `rbext(mu=muleft, phi=phileft, pex=pex_left, bex=1)`.
#'      The final value is `(1 - y_raw) * threshold`. Note `pex_left = (1-bex)*pex*2` and `bex=1` is used in `rbext`
#'      to ensure the extreme mass corresponds to CHOCO value 0.
#'    - **Right Side (`threshold` to 1):** Simulate `y_raw` from `rbext(mu=muright, phi=phiright, pex=pex_right, bex=1)`.
#'      The final value is `threshold + y_raw * (1 - threshold)`. Note `pex_right = bex*pex*2` and `bex=1` is used
#'      in `rbext` to ensure the extreme mass corresponds to CHOCO value 1.
#' 
#' @references
#' - Kubinec, R. (2023). Ordered beta regression: a parsimonious, well-fitting model for continuous data with
#'     lower and upper bounds. Political Analysis, 31(4), 519-536. (Describes the underlying ordered beta model)
#'
#' @examples
#' # Simulate data with different parameterizations
#' # 10% at threshold, 50/50 split otherwise, symmetric means/precisions
#' x1 <- rchoco(n=5000, p = 0.5, muleft = 0.7, mudelta = 0, phileft = 5,
#'   phidelta = 0, pex = 0.1, bex = 0.5, pmid = 0.1, threshold = 0.5)
#' hist(x1, breaks = 50, main = "CHOCO: Symmetric, pmid=0.1", xlab = "y")
#'
#' # No threshold mass, 70% probability on right, different means/precisions
#' x2 <- rchoco(n=5000, p = 0.7, muleft = 0.8, mudelta = -1, phileft = 3,
#'   phidelta = 0.5, pex = 0.05, bex = 0.7, pmid = 0, threshold = 0.5)
#' hist(x2, breaks = 50, main = "CHOCO: Asymmetric, pmid=0", xlab = "y")
#' @rdname rchoco
#' @export
rchoco <- function(n,
                   p = 0.5,
                   muleft = 0.3,
                   mudelta = 0,
                   phileft = 4,
                   phidelta = 0,
                   pex = 0.1,
                   bex = 0.5,
                   pmid = 0,
                   threshold = 0.5) {

  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) stop("n must be a positive integer.")
  if (any(threshold <= 0 | threshold >= 1)) stop("threshold must be between 0 and 1 (exclusive).")
  if (any(pex < 0 | pex > 1)) stop("pex must be between 0 and 1.")
  if (any(bex < 0 | bex > 1)) stop("bex must be between 0 and 1.") # Added validation for bex
  if (any(p < 0 | p > 1)) stop("p must be between 0 and 1.") # Added validation for p
  if (any(pmid < 0 | pmid > 1)) stop("pmid must be between 0 and 1.") # Added validation for pmid
  if (any(muleft <= 0 | muleft >= 1)) stop("muleft must be between 0 and 1 (exclusive).")
  if (any(phileft <= 0)) stop("phileft must be positive.")

  # Tolerance for floating point comparisons
  eps <- 1e-9 

  # --- Vectorization ---
  # Determine output length based on n and parameter vector lengths
  param_lengths <- c(length(p), length(muleft), length(mudelta), length(phileft),
                     length(phidelta), length(pex), length(bex), length(pmid), length(threshold))
  n_params <- max(param_lengths)
  n_out <- max(n, n_params)

  # Recycle parameters to match output length
  p <- rep(p, length.out = n_out)
  muleft <- rep(muleft, length.out = n_out)
  mudelta <- rep(mudelta, length.out = n_out)
  phileft <- rep(phileft, length.out = n_out)
  phidelta <- rep(phidelta, length.out = n_out)
  pex <- rep(pex, length.out = n_out)
  bex <- rep(bex, length.out = n_out)
  pmid <- rep(pmid, length.out = n_out)
  threshold <- rep(threshold, length.out = n_out)

  # --- Parameter computation ---
  logit_muleft <- log(muleft / (1 - muleft))
  logit_muright <- logit_muleft + mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))
  # Clamp derived parameters to avoid issues at boundaries
  muright <- pmax(eps, pmin(muright, 1 - eps))
  muleft <- pmax(eps, pmin(muleft, 1 - eps)) # Also clamp original muleft

  phiright <- phileft * exp(phidelta)
  phiright <- pmax(eps, phiright) # Ensure phiright is positive

  # Effective extreme probabilities (clamped)
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2))) # Zeros
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))      # Ones

  # --- Simulation ---

  # 1. Decide for each trial: Left (0), Middle (0.5), or Right (1) component
  # Calculate probabilities for each component for each draw
  prob_left <- (1 - pmid) * (1 - p)
  prob_mid <- pmid
  prob_right <- (1 - pmid) * p

  # Sample side choice element-wise using uniform random numbers
  side_choice <- numeric(n_out)
  rand_unif <- stats::runif(n_out)
  # Assign based on cumulative probability thresholds
  side_choice[rand_unif < prob_left] <- 0
  side_choice[rand_unif >= prob_left & rand_unif < (prob_left + prob_mid)] <- 0.5
  # Anything >= prob_left + prob_mid falls into the right component
  side_choice[rand_unif >= (prob_left + prob_mid)] <- 1

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
      # Note: For the left side, the underlying Beta's "success" (value near 1)
      # corresponds to the final CHOCO value being near 0.
      # We use bex=1 in rbext to force extremes to be only at 1 (which corresponds to 0 after mirroring).
      y_raw_left <- rbext(n = n_left, mu = muleft[idx_left], phi = phileft[idx_left],
                          pex = pex_left[idx_left], bex = 1) # Use subset of params, bex=1
      # Mirror (1 - y_raw) and scale by threshold
      out[idx_left] <- (1 - y_raw_left) * threshold[idx_left] # Use vectorized threshold
  }

  # 4. Simulate Right Component (Rescaled)
  if (n_right > 0) {
      # Simulate underlying Beta-extreme values (mapped to [0,1])
      # Note: For the right side, the underlying Beta's "success" (value near 1)
      # corresponds to the final CHOCO value being near 1.
      # We use bex=1 in rbext to force extremes to be only at 1.
      y_raw_right <- rbext(n = n_right, mu = muright[idx_right], phi = phiright[idx_right],
                           pex = pex_right[idx_right], bex = 1) # Use subset of params, bex=1
      # Scale by (1 - threshold) and shift by threshold
      out[idx_right] <- threshold[idx_right] + y_raw_right * (1 - threshold[idx_right]) # Use vectorized threshold
  }

  # Return the combined results, ensuring correct length if n was the determining factor
  if (n_out == n) {
      return(out)
  } else {
      # This case happens if n=1 but parameters were vectors. Return the full vector.
      return(out)
  }
}



#' @rdname rchoco
#' @inheritParams rbext
#' @export
dchoco <- function(x, p = 0.5, muleft = 0.3, mudelta = 0, phileft = 4, phidelta = 0,
                   pex = 0.1, bex = 0.5, pmid = 0, threshold = 0.5, log = FALSE) {

  eps <- 1e-9 # Tolerance for floating point comparisons

  # --- Input Validation ---
  if (any(p < 0 | p > 1)) {
      warning("p must be between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(muleft <= 0 | muleft >= 1)) {
      warning("muleft must be strictly between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(phileft <= 0)) {
      warning("phileft must be positive. Returning 0 density / -Inf log-density.")
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
  if (any(pmid < 0 | pmid > 1)) {
      warning("pmid must be between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }
  if (any(threshold <= 0 | threshold >= 1) ) {
      warning("threshold must be strictly between 0 and 1. Returning 0 density / -Inf log-density.")
      return(ifelse(log, -Inf, 0))
  }

  # --- Vectorization ---
  n <- length(x)
  p <- rep(p, length.out = n)
  muleft <- rep(muleft, length.out = n)
  mudelta <- rep(mudelta, length.out = n)
  phileft <- rep(phileft, length.out = n)
  phidelta <- rep(phidelta, length.out = n)
  pex <- rep(pex, length.out = n)
  bex <- rep(bex, length.out = n)
  pmid <- rep(pmid, length.out = n)
  threshold <- rep(threshold, length.out = n)

  # --- Parameter computation ---
  logit_muleft <- log(muleft / (1 - muleft))
  logit_muright <- logit_muleft + mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))
  phiright <- phileft * exp(phidelta)
  # Clamp derived parameters to avoid issues at boundaries
  muright <- pmax(eps, pmin(muright, 1 - eps))
  muleft_clamped <- pmax(eps, pmin(muleft, 1 - eps)) # Use clamped version for dbeta
  phiright <- pmax(eps, phiright) # Ensure phiright is positive

  # Effective extreme probabilities (clamped)
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2))) # Zeros
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))      # Ones

  # --- Density Calculation ---
  density <- numeric(n)

  # Handle values outside [0, 1]
  outside_idx <- x < 0 | x > 1
  density[outside_idx] <- 0

  # Handle x = threshold
  thresh_idx <- abs(x - threshold) < eps & !outside_idx
  density[thresh_idx] <- pmid[thresh_idx]

  # Handle x = 0
  zero_idx <- abs(x - 0) < eps & !outside_idx & !thresh_idx
  density[zero_idx] <- (1 - pmid[zero_idx]) * (1 - p[zero_idx]) * pex_left[zero_idx]

  # Handle x = 1
  one_idx <- abs(x - 1) < eps & !outside_idx & !thresh_idx
  density[one_idx] <- (1 - pmid[one_idx]) * p[one_idx] * pex_right[one_idx]

  # Handle 0 < x < threshold
  left_cont_idx <- x > eps & x < threshold - eps & !outside_idx
  if (any(left_cont_idx)) {
    y_raw_left <- 1 - x[left_cont_idx] / threshold[left_cont_idx]
    # Clamp input to dbeta for stability
    y_raw_left <- pmax(eps, pmin(y_raw_left, 1 - eps))
    beta_dens_left <- stats::dbeta(y_raw_left,
                                   shape1 = muleft_clamped[left_cont_idx] * phileft[left_cont_idx] * 2,
                                   shape2 = (1 - muleft_clamped[left_cont_idx]) * phileft[left_cont_idx] * 2)
    # Jacobian for transformation: 1 / threshold
    density[left_cont_idx] <- (1 - pmid[left_cont_idx]) * (1 - p[left_cont_idx]) * (1 - pex_left[left_cont_idx]) *
                              beta_dens_left / threshold[left_cont_idx]
  }

  # Handle threshold < x < 1
  right_cont_idx <- x > threshold + eps & x < 1 - eps & !outside_idx
  if (any(right_cont_idx)) {
    y_raw_right <- (x[right_cont_idx] - threshold[right_cont_idx]) / (1 - threshold[right_cont_idx])
    # Clamp input to dbeta for stability
    y_raw_right <- pmax(eps, pmin(y_raw_right, 1 - eps))
    beta_dens_right <- stats::dbeta(y_raw_right,
                                    shape1 = muright[right_cont_idx] * phiright[right_cont_idx] * 2,
                                    shape2 = (1 - muright[right_cont_idx]) * phiright[right_cont_idx] * 2)
    # Jacobian for transformation: 1 / (1 - threshold)
    density[right_cont_idx] <- (1 - pmid[right_cont_idx]) * p[right_cont_idx] * (1 - pex_right[right_cont_idx]) *
                               beta_dens_right / (1 - threshold[right_cont_idx])
  }

  # Return log density if requested
  if (log) {
    log_density <- ifelse(density == 0, -Inf, log(density))
    # Ensure any NaN from dbeta becomes -Inf
    log_density[is.na(log_density)] <- -Inf
    return(log_density)
  } else {
    # Ensure any NaN from dbeta becomes 0
    density[is.na(density)] <- 0
    return(density)
  }
}