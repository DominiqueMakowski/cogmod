#' @title Discrete Beta Model
#'
#' @description
#' The Discrete Beta (DBT) distribution models ordinal rating data on a fixed
#' integer scale \eqn{R \in \{1, \dots, k\}} by discretizing an underlying
#' continuous Beta distribution at \eqn{k - 1} evenly-spaced thresholds
#' \eqn{\gamma_j = j / k}. Unlike proportional-odds style models, which fix
#' the underlying distribution and estimate the thresholds, the Discrete Beta
#' fixes the thresholds and estimates the two shape parameters of the
#' underlying Beta distribution instead. This keeps the model parsimonious
#' (only 2 parameters) while remaining flexible enough to reproduce "U" and "J"
#' (non-monotonic convex) shapes that are common in rating data and that
#' proportional-odds models cannot capture (Sciandra et al., 2024).
#'
#' @param n Number of simulated values.
#' @param x,q Vector of quantiles (integer ratings between 1 and `k`, or 0 if
#'   `pzero > 0`).
#' @param p Vector of probabilities.
#' @inheritParams rbetagate
#' @param k Number of rating categories (a positive integer, `k >= 1`), i.e.
#'   the response scale runs from 1 to `k`.
#' @param pzero Probability of an additional "hurdle" point mass at 0, on top of
#'   the `1:k` rating scale. Defaults to `0`, in which case the distribution
#'   reduces to the pure Discrete Beta model. Useful for rating scales that
#'   include an extra "zero" category (e.g., "not applicable" or a genuine
#'   zero response) that is not part of the underlying `1:k` continuum.
#' @param log,log.p Logical; if `TRUE`, probabilities/densities are returned
#'   on the log scale.
#' @param lower.tail Logical; if `TRUE` (default), probabilities are
#'   \eqn{P(R \le q)}, otherwise \eqn{P(R > q)}.
#'
#' @details
#' Writing \eqn{\alpha = \mu \phi} and \eqn{\beta = (1 - \mu)\phi} for the
#' shape parameters of the underlying Beta distribution, the probability mass
#' function is (Sciandra et al., 2024, eq. 2)
#' \deqn{P(R = j) = F_B(j/k; \alpha, \beta) - F_B((j-1)/k; \alpha, \beta), \quad j = 1, \dots, k}
#' where \eqn{F_B} is the Beta CDF.
#'
#' `rbetadiscrete()` uses the equivalent, faster generative representation:
#' draw a continuous \eqn{X \sim Beta(\alpha, \beta)} and set
#' \eqn{R = \lceil k X \rceil}, clipped to `[1, k]`.
#'
#' When `pzero > 0`, a hurdle is added at 0: with probability `pzero` the
#' response is 0, and with probability `1 - pzero` it is generated from the
#' Discrete Beta distribution described above, i.e.
#' \deqn{P(R = 0) = \code{pzero}, \quad P(R = j) = (1 - \code{pzero}) \times [F_B(j/k) - F_B((j-1)/k)], \quad j = 1, \dots, k}
#'
#' **Special cases:**
#' - `mu = 0.5`, `phi = 1` (i.e. `alpha = beta = 1`): reduces to the discrete
#'   Uniform distribution on `1:k`.
#' - `alpha, beta < 1`: "U"/"J"-shaped, with mass concentrated in the tails.
#' - `alpha, beta > 1`: concave, with mass concentrated around the middle
#'   category.
#' - `phi -> Inf` (with `mu` fixed): mass concentrates on a single category.
#' - `pzero = 0`: reduces to the pure Discrete Beta model (no hurdle).
#'
#' @return `dbetadiscrete()` returns the probability mass; `pbetadiscrete()`
#'   returns the cumulative probability; `qbetadiscrete()` returns the
#'   quantile (an integer between 0 and `k`); `rbetadiscrete()` returns
#'   simulated ratings. All are vectorized over `x`/`q`/`p`, `mu`, `phi`,
#'   `pzero` and `k`.
#'
#' @references
#' - Sciandra, M., Fasola, S., Albano, A., Di Maria, C., & Plaia, A. (2024).
#'   Discrete Beta and Shifted Beta-Binomial models for rating and ranking
#'   data. Environmental and Ecological Statistics, 31, 317-338.
#'   \doi{10.1007/s10651-023-00592-5}
#'
#' @examples
#' x <- 1:10
#' probs <- dbetadiscrete(x, mu = 0.66, phi = 3.51, k = 10)
#' # barplot(probs, names.arg = x)
#'
#' y <- rbetadiscrete(1000, mu = 0.66, phi = 3.51, k = 10)
#' # hist(y, breaks = 0:10)
#'
#' # discrete Uniform special case
#' dbetadiscrete(1:5, mu = 0.5, phi = 1, k = 5)
#'
#' # hurdle at zero: 20% chance of a 0, otherwise pure Discrete Beta
#' dbetadiscrete(0:5, mu = 0.66, phi = 3.51, k = 5, pzero = 0.2)
#'
#' @export
rbetadiscrete <- function(n, mu = 0.5, phi = 3, k = 5, pzero = 0) {
  # --- Input Validation ---
  if (any(n <= 0 | n != floor(n))) {
    stop("n must be a positive integer.")
  }
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n_out <- max(n, length(mu), length(phi), length(k), length(pzero))
  if (n_out > 1 || n > 1) {
    mu <- rep(mu, length.out = n_out)
    phi <- rep(phi, length.out = n_out)
    k <- rep(k, length.out = n_out)
    pzero <- rep(pzero, length.out = n_out)
  } else {
    n_out <- n
  }

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Simulation ---
  # Discretize a continuous Beta draw at the fixed thresholds gamma_j = j / k
  x <- stats::rbeta(n_out, shape1 = alpha, shape2 = beta)
  out <- ceiling(x * k)
  out <- pmax(1L, pmin(out, k)) # guard against floating-point edge cases at 0 and 1

  # Hurdle: overwrite with 0 wherever the hurdle is crossed
  is_zero <- stats::runif(n_out) < pzero
  out[is_zero] <- 0L

  if (n_out == n && n > 0) {
    return(out[seq_len(n)])
  }
  out
}


#' @rdname rbetadiscrete
#' @export
dbetadiscrete <- function(x, mu = 0.5, phi = 3, k = 5, pzero = 0, log = FALSE) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(x)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Probability Mass ---
  # P(R = j) = F_B(j/k) - F_B((j-1)/k); handled at the boundaries explicitly
  # so that floating-point pbeta(0, ...) / pbeta(1, ...) calls are avoided.
  upper <- ifelse(x >= k, 1, stats::pbeta(x / k, alpha, beta))
  lower <- ifelse(x <= 1, 0, stats::pbeta((x - 1) / k, alpha, beta))

  prob <- (1 - pzero) * (upper - lower)
  # Hurdle point mass at 0
  is_zero <- x == 0
  prob[is_zero] <- pzero[is_zero]

  invalid <- (x < 0 | x > k | x != floor(x))
  prob[invalid] <- 0
  prob <- pmax(0, prob) # guard against tiny negative values from cancellation

  if (log) {
    prob <- ifelse(prob > 0, log(prob), -Inf)
  }
  prob
}


#' @rdname rbetadiscrete
#' @export
pbetadiscrete <- function(
  q,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(q)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Cumulative Probability ---
  qq <- pmin(pmax(floor(q), 0), k)
  p_orig <- ifelse(
    qq <= 0,
    0,
    ifelse(qq >= k, 1, stats::pbeta(qq / k, alpha, beta))
  )
  # P(R <= q) = pzero + (1 - pzero) * F_B(q/k), with P(R <= 0) = pzero
  p <- pzero + (1 - pzero) * p_orig
  p[q < 0] <- 0

  if (!lower.tail) {
    p <- 1 - p
  }
  if (log.p) {
    p <- log(p)
  }
  p
}


#' @rdname rbetadiscrete
#' @export
qbetadiscrete <- function(
  p,
  mu = 0.5,
  phi = 3,
  k = 5,
  pzero = 0,
  lower.tail = TRUE,
  log.p = FALSE
) {
  # --- Input Validation ---
  if (any(mu <= 0 | mu >= 1)) {
    stop("mu must be strictly between 0 and 1.")
  }
  if (any(phi <= 0)) {
    stop("phi must be positive.")
  }
  if (any(k < 1 | k != floor(k))) {
    stop("k must be a positive integer (>= 1).")
  }
  if (any(pzero < 0 | pzero > 1)) {
    stop("pzero must be between 0 and 1.")
  }

  if (log.p) {
    p <- exp(p)
  }
  if (!lower.tail) {
    p <- 1 - p
  }
  if (any(p < 0 | p > 1)) {
    stop("p must be between 0 and 1.")
  }

  # --- Vectorization ---
  n <- length(p)
  mu <- rep(mu, length.out = n)
  phi <- rep(phi, length.out = n)
  k <- rep(k, length.out = n)
  pzero <- rep(pzero, length.out = n)

  # --- Parameter Calculation ---
  alpha <- mu * phi * 2
  beta <- (1 - mu) * phi * 2

  # --- Quantile ---
  # Because the thresholds gamma_j = j/k are fixed, the discretized quantile
  # is simply the continuous Beta quantile mapped through the same thresholds.
  # Below the hurdle mass, the quantile is 0; above it, rescale p onto the
  # [0, 1] scale of the underlying (non-hurdle) Discrete Beta distribution.
  p_resc <- pmax(0, (p - pzero) / (1 - pzero))
  xp <- stats::qbeta(p_resc, alpha, beta)
  out <- pmax(1L, pmin(k, ceiling(k * xp)))
  out[p <= pzero] <- 0L
  out
}
