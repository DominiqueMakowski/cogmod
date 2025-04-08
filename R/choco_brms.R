# Stanvars ----------------------------------------------------------------

# Stan Functions and Custom Family for brms

#' @keywords internal
.choco_lpdf <- function() {
"
real choco_lpdf(
  real y,          // Observed value (must be in [0, 1])
  real mu,         // Proportion of data on the right-hand side of the threshold (in [0, 1])
  real muleft,     // Mean of the Beta distribution for the left-hand side (in [0, 1])
  real mudelta,    // Deviation in mean for the right-hand side relative to the left-hand side
  real phileft,    // Precision parameter for the left-hand side Beta distribution (positive)
  real phidelta,   // Deviation in precision for the right-hand side (log scale)
  real pex,        // Probability of extreme values (zeros and ones) (in [0, 1])
  real bex         // Balance of extreme values (proportion of ones relative to zeros) (in [0, 1])
) {
  real eps = 1e-6;  // Small constant to avoid numerical issues
  real threshold = 0.5;  // Fixed threshold separating left and right Beta distributions

  // Special case: pex = 1 (all mass at extremes)
  if (pex >= 1 - eps) {
    real prob0 = (1 - mu) * (1 - bex) / ((1 - mu) * (1 - bex) + mu * bex);
    if (y < eps) {
      return log(prob0);
    } else if (y > 1 - eps) {
      return log(1 - prob0);
    } else {
      return negative_infinity();
    }
  }

  // Transform parameters for right-side distribution
  real logit_muleft = log(muleft / (1 - muleft));
  real logit_muright = -logit_muleft + mudelta;
  real muright = exp(logit_muright) / (1 + exp(logit_muright));
  real phiright = phileft * exp(phidelta);

  if (muright <= 0 || muright >= 1) {
    return negative_infinity();
  }

  // Compute probabilities for zeros and ones
  real zero_prob = pex * (1 - bex) * (1 - mu);
  real one_prob  = pex * bex * mu;
  if (zero_prob + one_prob > 1) {
    return negative_infinity();
  }

  // Handle exact zeros
  if (y < eps) {
    return log(zero_prob);
  }
  // Handle exact ones
  if (y > 1 - eps) {
    return log(one_prob);
  }

  // Continuous mass and scale factors
  real cont_mass = 1 - zero_prob - one_prob;
  real left_scale = cont_mass * (1 - mu);
  real right_scale = cont_mass * mu;

  // Handle threshold boundary
  if (abs(y - threshold) < eps) {
    real overlap = 1e-6;
    if (mu < eps) {
      // only left side
      real left_log_density = beta_lpdf(1.0 - overlap | muleft * phileft, (1 - muleft) * phileft);
      return log(left_scale) + left_log_density - log(threshold);
    } else if (mu > 1 - eps) {
      // only right side
      real right_log_density = beta_lpdf(overlap | muright * phiright, (1 - muright) * phiright);
      return log(right_scale) + right_log_density - log(1 - threshold);
    } else {
      real left_density = left_scale * exp(beta_lpdf(1.0 - overlap | muleft * phileft, (1 - muleft) * phileft)) / threshold;
      real right_density = right_scale * exp(beta_lpdf(overlap | muright * phiright, (1 - muright) * phiright)) / (1 - threshold);
      return log((1 - mu) * left_density + mu * right_density);
    }
  }

  // Left side
  if (y < threshold) {
    real scaled_y = y / threshold;
    return log(left_scale) + beta_lpdf(scaled_y | muleft * phileft, (1 - muleft) * phileft) - log(threshold);
  }
  // Right side
  else {
    real scaled_y = (y - threshold) / (1 - threshold);
    return log(right_scale) + beta_lpdf(scaled_y | muright * phiright, (1 - muright) * phiright) - log(1 - threshold);
  }
}
"
}

#' @rdname rchoco
#' @examples
#' # You can expose the lpdf function as follows:
#' # choco_lpdf <- choco_lpdf_expose()
#' # choco_lpdf(y = 0.5, mu = 0.5, muleft = 0.3, mudelta = 0.2,
#' #' #            phileft = 1, phidelta = 0.5, pex = 0.8, bex = 0.6)
#'
#' @export
choco_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .choco_lpdf(),
    "\n}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$choco_lpdf
}

#' @rdname rchoco
#' @export
choco_stanvars <- function() {
  brms::stanvar(scode = .choco_lpdf(), block = "functions")
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
    lb = c(0, 0, NA, 0, NA, 0, 0),
    type = "real"
  )
}

# brms --------------------------------------------------------------------



##' @rdname rchoco
#' @inheritParams rbext
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
  
  # Fixed threshold in the CHOCO model
  threshold <- 0.5
  
  # Get the number of posterior draws and observation sets
  n_draws <- length(mu)
  outcomes <- numeric(n_draws)
  
  # Small value for numerical stability
  eps <- 1e-6
  
  # Special case: pex ≈ 1 (all mass at extremes)
  pex_is_one <- pex >= (1 - eps)
  
  if (any(pex_is_one)) {
    # Handle observations with pex ≈ 1
    special_idx <- which(pex_is_one)
    for (j in special_idx) {
      # When pex = 1, outcomes are either 0 or 1
      prob0 <- (1 - mu[j]) * (1 - bex[j]) / ((1 - mu[j]) * (1 - bex[j]) + mu[j] * bex[j])
      outcomes[j] <- sample(c(0, 1), size = 1, prob = c(prob0, 1 - prob0))
    }
    
    # If all are special case, return
    if (length(special_idx) == n_draws) {
      return(as.matrix(outcomes))
    }
    
    # Otherwise, process the regular cases
    regular_idx <- which(!pex_is_one)
    mu <- mu[regular_idx]
    muleft <- muleft[regular_idx]
    mudelta <- mudelta[regular_idx]
    phileft <- phileft[regular_idx]
    phidelta <- phidelta[regular_idx]
    pex <- pex[regular_idx]
    bex <- bex[regular_idx]
  } else {
    regular_idx <- seq_along(mu)
  }
  
  # For each regular draw
  for (j in seq_along(mu)) {
    # Transform parameters for the right-side distribution
    logit_muleft <- log(muleft[j] / (1 - muleft[j]))
    logit_muright <- -logit_muleft + mudelta[j]
    muright <- exp(logit_muright) / (1 + exp(logit_muright))
    phiright <- phileft[j] * exp(phidelta[j])
    
    # Check if parameter combination is valid
    if (muright <= 0 || muright >= 1) {
      outcomes[regular_idx[j]] <- NA
      next
    }
    
    # Compute probabilities for zeros and ones
    zero_prob <- pex[j] * (1 - bex[j]) * (1 - mu[j])
    one_prob <- pex[j] * bex[j] * mu[j]
    
    # Check if extreme value probabilities are valid
    if (zero_prob + one_prob > 1) {
      outcomes[regular_idx[j]] <- NA
      next
    }
    
    # Draw uniform random number for classification
    u <- stats::runif(1)
    
    # Assign value based on the drawn uniform number
    if (u < zero_prob) {
      outcomes[regular_idx[j]] <- 0  # zero
    } else if (u >= (1 - one_prob)) {
      outcomes[regular_idx[j]] <- 1  # one
    } else {
      # Generate continuous outcome
      # Rescale uniform to determine left/right assignment
      scaled_u <- (u - zero_prob) / (1 - zero_prob - one_prob)
      
      if (scaled_u < (1 - mu[j])) {
        # Generate Beta variate for left side
        outcomes[regular_idx[j]] <- threshold * stats::rbeta(
          1,
          shape1 = muleft[j] * phileft[j],
          shape2 = (1 - muleft[j]) * phileft[j]
        )
      } else {
        # Generate Beta variate for right side
        outcomes[regular_idx[j]] <- threshold + (1 - threshold) * stats::rbeta(
          1,
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright
        )
      }
    }
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
  
  # Fixed threshold
  threshold <- 0.5
  
  # Initialize log likelihood vector
  ll <- numeric(length(mu))
  
  # Small value for numerical stability
  eps <- 1e-6
  
  # For each posterior draw
  for (j in seq_along(mu)) {
    # Special case: pex ≈ 1 (all mass at extremes)
    if (pex[j] >= (1 - eps)) {
      prob0 <- (1 - mu[j]) * (1 - bex[j]) / ((1 - mu[j]) * (1 - bex[j]) + mu[j] * bex[j])
      
      if (y < eps) {
        ll[j] <- log(prob0)
      } else if (y > (1 - eps)) {
        ll[j] <- log(1 - prob0)
      } else {
        ll[j] <- -Inf  # All other values have zero density
      }
      next
    }
    
    # Transform parameters for the right-side distribution
    logit_muleft <- log(muleft[j] / (1 - muleft[j]))
    logit_muright <- -logit_muleft + mudelta[j]
    muright <- exp(logit_muright) / (1 + exp(logit_muright))
    phiright <- phileft[j] * exp(phidelta[j])
    
    # Check if parameter combination is valid
    if (muright <= 0 || muright >= 1) {
      ll[j] <- -Inf
      next
    }
    
    # Compute probabilities for zeros and ones
    zero_prob <- pex[j] * (1 - bex[j]) * (1 - mu[j])
    one_prob <- pex[j] * bex[j] * mu[j]
    
    # Check if extreme value probabilities are valid
    if (zero_prob + one_prob > 1) {
      ll[j] <- -Inf
      next
    }
    
    # Handle exact zeros (point mass at 0)
    if (y < eps) {
      ll[j] <- log(zero_prob)
      next
    }
    
    # Handle exact ones (point mass at 1)
    if (y > (1 - eps)) {
      ll[j] <- log(one_prob)
      next
    }
    
    # Continuous mass and scale factors
    cont_mass <- 1 - zero_prob - one_prob
    left_scale <- cont_mass * (1 - mu[j])
    right_scale <- cont_mass * mu[j]
    
    # Handle values at exactly the threshold
    if (abs(y - threshold) < eps) {
      overlap <- 1e-6  # Small offset to compute limits
      
      if (mu[j] < eps) {
        # When mu = 0, only use left limit
        left_log_density <- stats::dbeta(
          1.0 - overlap, 
          shape1 = muleft[j] * phileft[j], 
          shape2 = (1 - muleft[j]) * phileft[j],
          log = TRUE
        )
        ll[j] <- log(left_scale) + left_log_density - log(threshold)
      } else if (mu[j] > (1 - eps)) {
        # When mu = 1, only use right limit
        right_log_density <- stats::dbeta(
          0.0 + overlap,
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright,
          log = TRUE
        )
        ll[j] <- log(right_scale) + right_log_density - log(1 - threshold)
      } else {
        # Calculate left and right limits
        left_density <- left_scale * exp(
          stats::dbeta(
            1.0 - overlap,
            shape1 = muleft[j] * phileft[j],
            shape2 = (1 - muleft[j]) * phileft[j],
            log = TRUE
          )
        ) / threshold
        
        right_density <- right_scale * exp(
          stats::dbeta(
            0.0 + overlap,
            shape1 = muright * phiright,
            shape2 = (1 - muright) * phiright,
            log = TRUE
          )
        ) / (1 - threshold)
        
        # Weighted average based on mu
        ll[j] <- log((1 - mu[j]) * left_density + mu[j] * right_density)
      }
      next
    }
    
    # Left side of threshold
    if (y < threshold) {
      scaled_y <- y / threshold
      ll[j] <- log(left_scale) + 
        stats::dbeta(
          scaled_y,
          shape1 = muleft[j] * phileft[j],
          shape2 = (1 - muleft[j]) * phileft[j],
          log = TRUE
        ) - log(threshold)
    }
    # Right side of threshold
    else {
      scaled_y <- (y - threshold) / (1 - threshold)
      ll[j] <- log(right_scale) + 
        stats::dbeta(
          scaled_y,
          shape1 = muright * phiright,
          shape2 = (1 - muright) * phiright,
          log = TRUE
        ) - log(1 - threshold)
    }
  }
  
  ll
}