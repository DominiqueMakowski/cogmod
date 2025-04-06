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
#' as the mirrored value of `muleft` around the threshold. A general parameter `mu` controls the
#' proportion of data on the right-hand side versus the left-hand side of the threshold.
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param mu General proportion of data on the right-hand side of the threshold. Must be in the range `[0, 1]`.
#'   If `mu = 0.3`, 30% of the data will be on the right-hand side.
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
#' x <- rchoco(10000, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5,
#'   phidelta = 0, pex = 0.1, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#'
#' x <- rchoco(10000, mu = 0.6, muleft = 0.4, mudelta = 0.4, phileft = 5,
#'   phidelta = 0, pex = 0.05, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#' @export
rchoco <- function(n, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 4, phidelta = 0, pex = 0.1, bex = 0.5, threshold = 0.5) {
  # Validate inputs
  if (mu < 0 || mu > 1)
    stop("mu must be between 0 and 1")
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

  # Compute muright using the logit link and mirroring
  logit_muleft <- log(muleft / (1 - muleft))  # Logit transformation of muleft
  logit_muright <- -logit_muleft + mudelta    # Mirror muleft and add mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))  # Inverse logit transformation

  # Validate muright
  if (muright <= 0 || muright >= 1)
    stop("Invalid parameterization: muright must be between 0 and 1")

  # Compute phiright
  phiright <- phileft * exp(phidelta)  # Adjust precision for the right-hand side

  # Compute kleft and kright from pex and bex
  kleft <- pex * (1 - bex)  # Probability of zeros
  kright <- pex * bex       # Probability of ones

  # Validate kleft and kright
  if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1 || kleft + kright > 1)
    stop("Invalid parameterization: ensure kleft + kright <= 1")

  # Draw uniform random numbers to determine the outcome category
  u <- stats::runif(n)

  # Initialize outcomes
  outcomes <- numeric(n)

  # Assign zeros and ones
  outcomes[u < kleft] <- 0            # Category 1: outcome = 0
  outcomes[u >= (1 - kright)] <- 1    # Category 3: outcome = 1

  # Identify continuous outcomes
  cont_idx <- which(u >= kleft & u < (1 - kright))

  # Draw from the Beta mixture for continuous outcomes
  if (length(cont_idx) > 0) {
    # Randomly assign to left or right Beta distribution based on `1 - mu`
    left_or_right <- stats::runif(length(cont_idx)) < (1 - mu)

    # Left-hand side Beta distribution
    outcomes[cont_idx[left_or_right]] <- threshold * stats::rbeta(
      sum(left_or_right),
      shape1 = muleft * phileft,
      shape2 = (1 - muleft) * phileft
    )

    # Right-hand side Beta distribution
    outcomes[cont_idx[!left_or_right]] <- threshold + (1 - threshold) * stats::rbeta(
      sum(!left_or_right),
      shape1 = muright * phiright,
      shape2 = (1 - muright) * phiright
    )
  }

  outcomes
}



#' @rdname rchoco
#' @param x A vector of values for which to compute the density.
#' @param log Logical; if TRUE, returns the log-density.
#' @export
dchoco <- function(x, mu = 0.5, muleft = 0.3, mudelta = 0, phileft = 5, phidelta = 0.5, pex = 0.1, bex = 0.5, threshold = 0.5, log = FALSE) {
  # Validate inputs
  if (any(x < 0 | x > 1))
    stop("x must be in the range [0, 1]")
  if (mu < 0 || mu > 1)
    stop("mu must be between 0 and 1")
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

  # Compute muright using the logit link and mirroring
  logit_muleft <- log(muleft / (1 - muleft))  # Logit transformation of muleft
  logit_muright <- -logit_muleft + mudelta    # Mirror muleft and add mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))  # Inverse logit transformation

  # Validate muright
  if (muright <= 0 || muright >= 1)
    stop("Invalid parameterization: muright must be between 0 and 1")

  # Compute phiright
  phiright <- phileft * exp(phidelta)  # Adjust precision for the right-hand side

  # Compute kleft and kright from pex and bex
  kleft <- pex * (1 - bex)  # Probability of zeros
  kright <- pex * bex       # Probability of ones
  p_cont <- 1 - pex         # Probability of continuous outcomes

  # Validate kleft and kright
  if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1 || kleft + kright > 1)
    stop("Invalid parameterization: ensure kleft + kright <= 1")

  # Initialize density vector
  density <- numeric(length(x))

  # Compute density for each region
  for (i in seq_along(x)) {
    if (x[i] == 0) {
      # Point mass at 0
      density[i] <- kleft
    } else if (x[i] == 1) {
      # Point mass at 1
      density[i] <- kright
    } else if (x[i] < threshold) {
      # Left-hand side Beta distribution
      beta_density <- stats::dbeta(x[i] / threshold, shape1 = muleft * phileft, shape2 = (1 - muleft) * phileft)
      density[i] <- p_cont * (1 - mu) * beta_density / threshold
    } else {
      # Right-hand side Beta distribution
      beta_density <- stats::dbeta((x[i] - threshold) / (1 - threshold), shape1 = muright * phiright, shape2 = (1 - muright) * phiright)
      density[i] <- p_cont * mu * beta_density / (1 - threshold)
    }
  }

  # Normalize the density to ensure it integrates to 1
  density <- density / sum(density * (1 / length(x)))

  # Return log-density if requested
  if (log) {
    return(log(density))
  } else {
    return(density)
  }
}



# Stanvars ----------------------------------------------------------------

# Stan Functions and Custom Family for brms

#' @rdname rchoco
#' @export
choco_stanvars <- function() {
  brms::stanvar(scode = "
real choco_lpdf(real y, real mu, real muleft, real mudelta, real phileft, real phidelta, real pex, real bex) {
  real eps = 1e-6;  // Small constant to avoid numerical issues
  real tol = 1e-6;  // Tolerance for floating-point comparisons
  real threshold = 0.5;  // Fixed threshold

  // Compute muright using the logit link and mirroring
  real logit_muleft = log(muleft / (1 - muleft));
  real logit_muright = -logit_muleft + mudelta;
  real muright = exp(logit_muright) / (1 + exp(logit_muright));

  // Validate muright
  if (muright <= eps || muright >= 1 - eps)
    return negative_infinity();

  // Compute phiright
  real phiright = phileft * exp(phidelta);

  // Compute probability masses
  real p0 = pex * (1 - bex);  // Mass at 0
  real p1 = pex * bex;        // Mass at 1
  real p_cont = 1 - pex;      // Mass for continuous outcomes

  // Precompute log values for efficiency
  real log_p0 = log(p0);
  real log_p1 = log(p1);
  real log_p_cont = log(p_cont);
  real log_mu = log(mu);
  real log1m_mu = log1m(mu);

  // Check valid parameterization
  if (p0 < 0 || p0 > 1 || p1 < 0 || p1 > 1 || abs(p0 + p_cont + p1 - 1) > tol)
    return negative_infinity();  // Invalid parameterization

  if (y < eps) {
    // Log-density for zeros
    return log_p0;
  } else if (y > (1 - eps)) {
    // Log-density for ones
    return log_p1;
  } else if (y < threshold) {
    // Log-density for left-hand side Beta distribution
    return log_p_cont + log1m_mu +
           beta_lpdf(y / threshold | muleft * phileft, (1 - muleft) * phileft) -
           log(threshold);
  } else {
    // Log-density for right-hand side Beta distribution
    return log_p_cont + log_mu +
           beta_lpdf((y - threshold) / (1 - threshold) | muright * phiright, (1 - muright) * phiright) -
           log(1 - threshold);
  }
}
", block = "functions")
}


#' @rdname rchoco
#' @param link_mu,link_muleft,link_mudelta,link_phileft,link_phidelta,link_pex,link_bex Link functions for the parameters.
#' @export
choco <- function(link_mu = "logit", link_muleft = "logit", link_mudelta = "identity",
                  link_phileft = "softplus", link_phidelta = "identity",
                  link_pex = "logit", link_bex = "logit") {
  brms::custom_family(
    name = "choco",
    dpars = c("mu", "muleft", "mudelta", "phileft", "phidelta", "pex", "bex"),
    links = c(link_mu, link_muleft, link_mudelta, link_phileft, link_phidelta, link_pex, link_bex),
    type = "real"
  )
}

# brms --------------------------------------------------------------------



#' @rdname rchoco
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_choco <- function(i, prep, ...) {
  # Extract distributional parameters for draw i
  mu <- brms::get_dpar(prep, "mu", i = i)
  muleft <- brms::get_dpar(prep, "muleft", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  phileft <- brms::get_dpar(prep, "phileft", i = i)
  phidelta <- brms::get_dpar(prep, "phidelta", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  # Compute muright using the logit link and mirroring
  logit_muleft <- log(muleft / (1 - muleft))
  logit_muright <- -logit_muleft + mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))

  # Compute phiright
  phiright <- phileft * exp(phidelta)

  # Compute probabilities for the three categories
  kleft <- pex * (1 - bex)  # Probability of zeros
  kright <- pex * bex       # Probability of ones
  p_cont <- 1 - pex         # Probability of continuous outcomes

  # Simulate outcomes
  n <- length(mu)
  u <- stats::runif(n)
  outcomes <- numeric(n)

  # Assign zeros and ones
  outcomes[u < kleft] <- 0
  outcomes[u >= (1 - kright)] <- 1

  # Simulate continuous outcomes
  cont_idx <- which(u >= kleft & u < (1 - kright))
  if (length(cont_idx) > 0) {
    left_or_right <- stats::runif(length(cont_idx)) < (1 - mu[cont_idx])
    outcomes[cont_idx[left_or_right]] <- 0.5 * stats::rbeta(
      sum(left_or_right),
      shape1 = muleft[cont_idx[left_or_right]] * phileft[cont_idx[left_or_right]],
      shape2 = (1 - muleft[cont_idx[left_or_right]]) * phileft[cont_idx[left_or_right]]
    )
    outcomes[cont_idx[!left_or_right]] <- 0.5 + 0.5 * stats::rbeta(
      sum(!left_or_right),
      shape1 = muright[cont_idx[!left_or_right]] * phiright[cont_idx[!left_or_right]],
      shape2 = (1 - muright[cont_idx[!left_or_right]]) * phiright[cont_idx[!left_or_right]]
    )
  }

  as.matrix(outcomes)
}

#' @rdname rchoco
#' @export
log_lik_choco <- function(i, prep) {
  # Extract observed value
  y <- prep$data$Y[i]

  # Extract distributional parameters for draw i
  mu <- brms::get_dpar(prep, "mu", i = i)
  muleft <- brms::get_dpar(prep, "muleft", i = i)
  mudelta <- brms::get_dpar(prep, "mudelta", i = i)
  phileft <- brms::get_dpar(prep, "phileft", i = i)
  phidelta <- brms::get_dpar(prep, "phidelta", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  # Compute muright using the logit link and mirroring
  logit_muleft <- log(muleft / (1 - muleft))
  logit_muright <- -logit_muleft + mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))

  # Compute phiright
  phiright <- phileft * exp(phidelta)

  # Compute probabilities for the three categories
  kleft <- pex * (1 - bex)  # Probability of zeros
  kright <- pex * bex       # Probability of ones
  p_cont <- 1 - pex         # Probability of continuous outcomes

  eps <- 1e-6  # Small constant for numerical stability

  # Compute log-likelihood
  if (y < eps) {
    return(log(kleft))
  } else if (y > (1 - eps)) {
    return(log(kright))
  } else if (y < 0.5) {
    beta_density <- stats::dbeta(y / 0.5, shape1 = muleft * phileft, shape2 = (1 - muleft) * phileft, log = TRUE)
    return(log(p_cont) + log(1 - mu) + beta_density - log(0.5))
  } else {
    beta_density <- stats::dbeta((y - 0.5) / 0.5, shape1 = muright * phiright, shape2 = (1 - muright) * phiright, log = TRUE)
    return(log(p_cont) + log(mu) + beta_density - log(0.5))
  }
}
