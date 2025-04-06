# Stanvars ----------------------------------------------------------------

# Stan Functions and Custom Family for brms

#' @rdname rbext
#' @export
bext_stanvars <- function() {
  brms::stanvar(scode = "
real bext_lpdf(real y, real mu, real phi, real pex, real bex) {
  real eps = 1e-8;  // Small constant to avoid numerical issues
  
  // Special case: pex = 1 (discrete Bernoulli distribution)
  if (abs(pex - 1) < eps) {
    if (y < eps) {
      // Log-density for zeros
      return log(1 - bex);
    } else if (y > (1 - eps)) {
      // Log-density for ones
      return log(bex);
    } else {
      // All other values have zero density
      return negative_infinity();
    }
  }
  
  // Regular case: mixture of continuous and discrete
  // Compute kleft and kright from pex and bex
  real kleft = pex * (1 - bex);    // probability mass at 0
  real kright = 1 - (pex * bex);   // threshold above which outcomes are set to 1
  
  // Validate parameters
  if (kleft < 0 || kleft > 1 || kright < 0 || kright > 1 || kleft >= kright)
    return negative_infinity();
    
  if (y < eps) {
    // Log-density for zeros
    return log(kleft);
  } else if (y > (1 - eps)) {
    // Log-density for ones
    return log(1 - kright);
  } else {
    // Log-density for continuous outcomes
    return log(kright - kleft) + beta_lpdf(y | mu * phi, (1 - mu) * phi);
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
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_bext <- function(i, prep, ...) {
  args <- list(...)
  type <- if (!is.null(args$type) && args$type == "continuous") "continuous" else "combined"

  # Extract draws for each parameter
  mu    <- brms::get_dpar(prep, "mu", i = i)
  phi   <- brms::get_dpar(prep, "phi", i = i)
  pex   <- brms::get_dpar(prep, "pex", i = i)
  bex   <- brms::get_dpar(prep, "bex", i = i)

  # Special case: pex = 1 (Bernoulli distribution)
  eps <- 1e-8
  pex_is_one <- abs(pex - 1) < eps
  
  if (any(pex_is_one) && type != "continuous") {
    final_out <- numeric(length(mu))
    
    # Handle the special case entries
    special_idx <- which(pex_is_one)
    for (j in special_idx) {
      # Sample from Bernoulli(bex)
      final_out[j] <- sample(c(0, 1), size = 1, prob = c(1-bex[j], bex[j]))
    }
    
    # If all are special case, return
    if (length(special_idx) == length(mu)) {
      return(as.matrix(final_out))
    }
    
    # Otherwise, process the regular cases
    regular_idx <- which(!pex_is_one)
    mu <- mu[regular_idx]
    phi <- phi[regular_idx]
    pex <- pex[regular_idx]
    bex <- bex[regular_idx]
  } else {
    final_out <- numeric(length(mu))
    regular_idx <- seq_along(mu)
  }
  
  # For regular cases (pex < 1)
  # Compute kleft and kright from pex and bex
  kleft <- pex * (1 - bex)   # threshold below which outcomes are set to 0
  kright <- 1 - (pex * bex)  # threshold above which outcomes are set to 1

  # Simulate continuous outcomes from the beta distribution
  y_beta <- stats::rbeta(n = length(mu),
                         shape1 = mu * phi,
                         shape2 = (1 - mu) * phi)

  if (type == "continuous") {
    # Return only the continuous outcomes
    return(as.matrix(y_beta))
  } else {
    # Generate uniform random numbers
    u <- stats::runif(length(mu))
    
    for (j in seq_along(mu)) {
      if (u[j] < kleft[j]) {
        final_out[regular_idx[j]] <- 0  # zero
      } else if (u[j] >= kright[j]) {
        final_out[regular_idx[j]] <- 1  # one
      } else {
        final_out[regular_idx[j]] <- y_beta[j]  # continuous outcome
      }
    }
  }

  as.matrix(final_out)
}

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

  eps <- 1e-8  # small constant for numerical stability
  
  # Initialize log-likelihood vector
  ll <- rep(NA, length(mu))
  
  # Special case: pex = 1 (discrete Bernoulli distribution)
  pex_is_one <- abs(pex - 1) < eps
  
  # Handle the special case (Bernoulli distribution)
  if (y < eps) {
    # Log-density for zeros
    ll[pex_is_one] <- log(1 - bex[pex_is_one])
  } else if (y > (1 - eps)) {
    # Log-density for ones
    ll[pex_is_one] <- log(bex[pex_is_one]) 
  } else {
    # All other values have zero density
    ll[pex_is_one] <- -Inf
  }
  
  # Handle the regular case (mixture of continuous and discrete)
  regular_idx <- which(!pex_is_one)
  
  if (length(regular_idx) > 0) {
    # Compute kleft and kright from pex and bex
    kleft <- pex[regular_idx] * (1 - bex[regular_idx])    # probability mass at 0
    kright <- 1 - (pex[regular_idx] * bex[regular_idx])   # threshold above which outcomes are set to 1
    
    # Validate parameters
    valid_params <- kleft >= 0 & kleft <= 1 & kright >= 0 & kright <= 1 & kleft < kright
    kleft[!valid_params] <- NA  # Mark invalid params
    
    if (y < eps) {
      # Log-density for zeros
      ll[regular_idx] <- log(kleft)
    } else if (y > (1 - eps)) {
      # Log-density for ones
      ll[regular_idx] <- log(1 - kright)
    } else {
      # Log-density for continuous outcomes
      ll[regular_idx] <- log(kright - kleft) + stats::dbeta(y,
                                                shape1 = mu[regular_idx] * phi[regular_idx],
                                                shape2 = (1 - mu[regular_idx]) * phi[regular_idx],
                                                log = TRUE)
    }
    
    # Set invalid parameter combinations to -Inf
    ll[regular_idx][!valid_params] <- -Inf
  }
  
  ll
}