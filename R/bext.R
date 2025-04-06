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
#'
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

  # Special case: pex = 1
  if (pex == 1) {
    if (bex == 0) {
      return(rep(0, n))  # All outcomes are zeros
    } else if (bex == 1) {
      return(rep(1, n))  # All outcomes are ones
    }
  }

  # Compute kleft and kright from pex and bex.
  # kleft: threshold below which outcomes are set to 0 (zeros)
  # kright: threshold above which outcomes are set to 1 (ones)
  kleft <- pex * (1 - bex)
  kright <- 1 - (pex * bex)

  # Validate kleft and kright
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

  # Compute the mass thresholds
  kleft <- pex * (1 - bex)       # Mass at 0
  kright <- 1 - (pex * bex)       # Upper threshold for continuous outcomes

  # If mu or phi are scalar, replicate to match length(x)
  if (length(mu) == 1) mu <- rep(mu, length(x))
  if (length(phi) == 1) phi <- rep(phi, length(x))

  # Compute the beta density for continuous outcomes
  beta_density <- mapply(function(xi, mui, phii) {
    stats::dbeta(xi, shape1 = mui * phii, shape2 = (1 - mui) * phii)
  }, x, mu, phi)

  # Compute the overall density
  density_vals <- sapply(seq_along(x), function(i) {
    if (x[i] == 0) {
      # Point mass at 0
      return(kleft)
    } else if (x[i] == 1) {
      # Point mass at 1; note: mass at 1 = pex * bex
      return(pex * bex)
    } else {
      # Continuous part scaled by the total probability for non-extremes
      # (kright - kleft) equals 1 - pex
      return((kright - kleft) * beta_density[i])
    }
  })

  # Return on log scale if requested
  if (log) {
    return(log(density_vals))
  } else {
    return(density_vals)
  }
}




# Stanvars ----------------------------------------------------------------

# Stan Functions and Custom Family for brms

#' @rdname rbext
#' @export
bext_stanvars <- function() {
  brms::stanvar(scode = "
real bext_lpdf(real y, real mu, real phi, real pex, real bex) {
  real eps = 1e-8;  // Small constant to avoid numerical issues
  // Compute probability masses
  real p0 = pex * (1 - bex);      // mass at 0
  real p1 = pex * bex;            // mass at 1
  real p_cont = 1 - pex;          // mass for continuous outcomes

  // Check valid parameterization: p0 and p1 must be in [0,1] and p0 + p_cont + p1 = 1
  if (p0 < 0 || p0 > 1 || p1 < 0 || p1 > 1)
    return negative_infinity();

  if (y < eps) {
    // Log-density for zeros
    return log(p0);
  } else if (y > (1 - eps)) {
    // Log-density for ones
    return log(p1);
  } else {
    // Log-density for continuous outcomes
    return log(p_cont) + beta_lpdf(y | mu * phi, (1 - mu) * phi);
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
                      links  = c("logit", "softplus", "logit", "logit"),
                      lb     = c(NA, 0, NA, NA),
                      type   = "real")
}

# brms --------------------------------------------------------------------

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

  # Compute probabilities for the three categories
  p0 <- pex * (1 - bex)  # probability of 0
  p1 <- pex * bex        # probability of 1
  p_cont <- 1 - pex      # probability of a continuous outcome

  # Simulate continuous outcomes from the beta distribution
  y_beta <- stats::rbeta(n = length(mu),
                         shape1 = mu * phi,
                         shape2 = (1 - mu) * phi)

  if (type == "continuous") {
    # Return only the continuous outcomes
    final_out <- y_beta
  } else {
    # Simulate outcomes by sampling category for each observation
    outcomes <- sapply(1:length(mu), function(j) {
      sample(1:3, size = 1, prob = c(p0[j], p_cont[j], p1[j]))
    })

    # Assign outcomes based on the sampled category
    final_out <- numeric(length(mu))
    for (j in 1:length(mu)) {
      if (outcomes[j] == 1) {
        final_out[j] <- 0  # zero
      } else if (outcomes[j] == 2) {
        final_out[j] <- y_beta[j]  # continuous outcome
      } else {
        final_out[j] <- 1  # one
      }
    }
  }

  as.matrix(final_out)
}


# Log Likelihood Function

#' @rdname rbext
#' @export
log_lik_bext <- function(i, prep) {
  # Extract observed value
  y <- prep$data$Y[i]

  # Extract model draws
  mu  <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pex <- brms::get_dpar(prep, "pex", i = i)
  bex <- brms::get_dpar(prep, "bex", i = i)

  # Compute probability masses
  p0 <- pex * (1 - bex)  # mass at 0
  p1 <- pex * bex        # mass at 1
  p_cont <- 1 - pex      # mass for continuous outcomes

  eps <- 1e-6  # small constant for numerical stability

  # Compute log-likelihood based on the value of y
  if (y < eps) {
    return(log(p0))
  } else if (y > (1 - eps)) {
    return(log(p1))
  } else {
    return(log(p_cont) + stats::dbeta(y,
                               shape1 = mu * phi,
                               shape2 = (1 - mu) * phi,
                               log = TRUE))
  }
}
