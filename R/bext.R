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
#' uses the total probability of extreme values (0 and 1) and their balance (i.e., the
#' relative proportion of zeros and ones).
#'
#'
#' @param n Number of simulated trials. Must be a positive integer.
#' @param mu Mean of the continuous component (between 0 and 1).
#' @param phi Precision parameter (positive).
#' @param pex Total probability of extreme values (zeros and ones). Must be in the range 0-1.
#' @param bex Balance of extreme values (identity link, between 0 and 1). Represents the proportion of ones relative to zeros. Default: 0.5.
#'
#' @return A vector of simulated outcomes in the range 0-1.
#'
#' @examples
#' # Simulate data with different parameterizations
#' x <- rbext(10000, mu = 0.5, phi = 3, pex = 0.5, bex = 0.5)
#' hist(x, breaks = 50, main = "Simulated Outcomes", xlab = "y")
#'
#' @export
rbext <- function(n, mu = 0.5, phi = 3, pex = 0.1, bex = 0.5) {
  # Validate inputs
  if (pex < 0 || pex > 1)
    stop("pex must be between 0 and 1")
  if (bex < 0 || bex > 1)
    stop("bex must be between 0 and 1")

  # Compute kright and kleft from pex and bex
  kright <- pex * bex       # Proportion of `pex` allocated to ones
  kleft <- pex * (1 - bex)  # Proportion of `pex` allocated to zeros

  # Validate kleft and kright
  if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1)
    stop("Invalid parameterization: kleft and kright must be in [0, 1]")

  # Draw uniform random numbers to determine the outcome category
  u <- stats::runif(n)

  # Generate outcomes
  outcomes <- numeric(n)
  outcomes[u < kleft] <- 0  # Category 1: outcome = 0
  outcomes[u >= (1 - kright)] <- 1  # Category 3: outcome = 1

  # Category 2: continuous outcomes
  cont_idx <- which(u >= kleft & u < (1 - kright))
  if (length(cont_idx) > 0) {
    outcomes[cont_idx] <- stats::rbeta(length(cont_idx),
                                       shape1 = mu * phi,
                                       shape2 = (1 - mu) * phi)
  }

  outcomes
}


# --------------------------------------------------------------------
# Stan Functions and Custom Family for brms

#' @rdname rbext
#' @export
bext_stanvars <- function() {
  brms::stanvar(scode = "
  real bext_lpdf(real y, real mu, real phi, real pex, real bex) {
    real eps = 1e-8;  // Small constant to avoid numerical issues
    real kright = pex * bex;       // Probability of ones
    real kleft = pex * (1 - bex);  // Probability of zeros

    if (y < eps) {
      // Log-density for zeros
      return log(kleft);
    } else if (y > (1 - eps)) {
      // Log-density for ones
      return log(kright);
    } else {
      // Log-density for continuous outcomes
      real beta_shape1 = mu * phi;
      real beta_shape2 = (1 - mu) * phi;
      return log(1 - pex) +  // Probability of continuous outcomes
             beta_lpdf(y | beta_shape1, beta_shape2);
    }
  }
  ", block = "functions")
}

#' @rdname rbext
#' @export
bext <- function() {
  # The custom family uses:
  # - mu: logit link (to constrain it to (0,1)),
  # - phi: log link (to ensure positivity),
  # - pex: logit link (to constrain it to (0,1)),
  # - bex: logit link (to constrain it to (0,1)).
  brms::custom_family("bext",
                      dpars = c("mu", "phi", "pex", "bex"),
                      links  = c("logit", "log", "logit", "logit"),
                      lb     = c(NA, 0, NA, NA),
                      type   = "real")
}

# --------------------------------------------------------------------
# Posterior Prediction

#' @rdname rbext
#' @inheritParams choco
#' @export
posterior_predict_bext <- function(i, prep, ...) {
  args <- list(...)
  type <- if (!is.null(args$type) && args$type == "continuous") "continuous" else "combined"

  # Extract draws for each parameter
  mu    <- brms::get_dpar(prep, "mu", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  pex   <- brms::get_dpar(prep, "pex", i = i)
  bex   <- brms::get_dpar(prep, "bex", i = i)

  # Compute kright and kleft from pex and bex
  kright <- pex * bex       # Proportion of `pex` allocated to ones
  kleft <- pex * (1 - bex)  # Proportion of `pex` allocated to zeros

  # Simulate continuous outcomes
  y_beta <- stats::rbeta(n = length(mu),
                         shape1 = mu * phi,
                         shape2 = (1 - mu) * phi)

  if (type == "continuous") {
    # Return only continuous outcomes
    final_out <- y_beta
  } else {
    # Simulate outcomes
    outcomes <- sapply(1:length(mu), function(j) {
      # Sample from the three categories: zeros, continuous, ones
      sample(1:3, size = 1, prob = c(kleft[j], 1 - pex[j], kright[j]))
    })

    # Assign outcomes based on the sampled category
    final_out <- numeric(length(mu))
    for (j in 1:length(mu)) {
      if (outcomes[j] == 1) {
        final_out[j] <- 0  # Zero
      } else if (outcomes[j] == 2) {
        final_out[j] <- y_beta[j]  # Continuous
      } else {
        final_out[j] <- 1  # One
      }
    }
  }

  as.matrix(final_out)
}

# --------------------------------------------------------------------
# Log Likelihood Function

#' @rdname rbext
#' @export
log_lik_bext <- function(i, prep) {
  y <- prep$data$Y[i]

  mu    <- brms::get_dpar(prep, "mu", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  pex   <- brms::get_dpar(prep, "pex", i = i)
  bex   <- brms::get_dpar(prep, "bex", i = i)

  # Compute kright and kleft from pex and bex
  kright <- pex * bex
  kleft <- pex * (1 - bex)

  # Compute log-likelihood
  ll <- if (y == 0) {
    log(kleft)
  } else if (y == 1) {
    log(kright)
  } else {
    log(1 - pex) +  # Probability of continuous outcomes
      stats::dbeta(y,
                   shape1 = mu * phi,
                   shape2 = (1 - mu) * phi,
                   log = TRUE)
  }

  ll
}
