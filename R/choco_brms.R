# Stan Functions and Custom Family for brms (CHOCO)

#' @keywords internal
.choco_lpdf <- function() {
"
// Log probability density function for the CHOCO distribution (Optimized)
// Accounts for point masses at 0, threshold (fixed at 0.5), and 1, and continuous Beta components.
real choco_lpdf(real y, real mu, real muleft, real mudelta, real phileft,
                 real phidelta, real pex, real bex, real pmid) {

  real eps = 1e-8; // Tolerance for floating point comparisons
  real threshold = 0.5; // Hardcoded threshold
  real log_lik; // Will hold our final log likelihood

  // Validate outcome y (essential)
  if (y < 0.0 || y > 1.0) {
    print(\"Error 1: outcome validation issue\");
    return negative_infinity();
  }

  // --- 1. Parameter Clamping ---
  real p_c = fmax(eps, fmin(mu, 1.0 - eps));
  real muleft_c = fmax(eps, fmin(muleft, 1.0 - eps));
  real phileft_c = fmax(eps, phileft);
  real pex_c = fmax(0.0, fmin(pex, 1.0)); // Allow 0 and 1
  real bex_c = fmax(0.0, fmin(bex, 1.0)); // Allow 0 and 1
  real pmid_c = fmax(0.0, fmin(pmid, 1.0)); // Allow 0 and 1

  // --- 3. Derived Parameter Computation ---
  // Calculate right-side parameters
  real logit_muleft = log(muleft_c / (1.0 - muleft_c));
  real logit_muright = logit_muleft + mudelta;
  real muright;
  // More inverse logit
  if (logit_muright > 15.0) {
    muright = 1.0 - eps;
  } else if (logit_muright < -15.0) {
    muright = eps;
  } else {
    muright = inv_logit(logit_muright);
  }
  real phiright = fmax(eps, phileft_c * exp(fmin(10.0, fmax(-10.0, phidelta))));

  // Calculate extreme probabilities for left and right
  real pex_left = fmin(1.0, fmax(0.0, (1.0 - bex_c) * (pex_c * 2.0)));
  real pex_right = fmin(1.0, fmax(0.0, bex_c * (pex_c * 2.0)));
  
  // Calculate Beta shape parameters
  real precision_left = phileft_c * 2.0;
  real precision_right = phiright * 2.0;
  real shape1_left = muleft_c * precision_left;
  real shape2_left = (1.0 - muleft_c) * precision_left;
  real shape1_right = muright * precision_right;
  real shape2_right = (1.0 - muright) * precision_right;
  
  // --- 3. Calculate component probabilities ---
  // These are the three main probabilities that divide the scale
  real prob_left = (1.0 - pmid_c) * (1.0 - p_c);
  real prob_mid = pmid_c;
  real prob_right = (1.0 - pmid_c) * p_c;
  
  // --- 4. Calculate log-density based on where y falls ---
  
  // Case 1: y is exactly at threshold
  if (abs(y - threshold) < eps) {
    // Handle the case where pmid is effectively zero
    if (pmid_c < eps) {
      return negative_infinity();
    }
    return log(pmid_c);
  }
  
  // Case 2: pmid = 1 but y is not at threshold
  // This is a special case where all probability mass is at the threshold (pmid=1)
  // In this case, any other value should have zero probability
  if (abs(pmid_c - 1.0) < eps) {
    return negative_infinity();
  }
  
  // Case 3: y is at 0
  else if (abs(y) < eps) {
    // Edge case: If pex is very close to 0, return -Inf
    if (pex_c < eps || pex_left < eps) {
      return negative_infinity();
    }

    // P(y=0) = P(left) * P(extreme=0|left)
    real p_zero = prob_left * pex_left;
    if (p_zero < eps) return negative_infinity();
    return log(p_zero);
  }
  
  // Case 4: y is at 1
  else if (abs(y - 1.0) < eps) {
    // Edge case: If pex is very close to 0, return -Inf
    if (pex_c < eps || pex_right < eps) {
      return negative_infinity();
    }

    // P(y=1) = P(right) * P(extreme=1|right)
    real p_one = prob_right * pex_right;
    if (p_one < eps) return negative_infinity();
    return log(p_one);
  }
  
  // Case 5: y is between 0 and threshold
  else if (y < threshold) {
    // Bail out early if the whole left component has zero probability
    if (prob_left < eps || 1.0 - pex_left < eps) {
      return negative_infinity(); 
    }
    
    // Transform y to raw Beta scale and compute density
    real y_raw_left = 1.0 - y / threshold;
    // Clamp for stability
    y_raw_left = fmax(eps, fmin(y_raw_left, 1.0 - eps));
    
    // Log-density = log(P(left) * P(continuous|left) * beta_density * jacobian)
    return log(prob_left) + log1m(pex_left) + 
           beta_lpdf(y_raw_left | shape1_left, shape2_left) + 
           log(1.0 / threshold);
  }
  
  // Case 6: y is between threshold and 1
  else {
    // Bail out early if the whole right component has zero probability
    if (prob_right < eps || 1.0 - pex_right < eps) {
      return negative_infinity(); 
    }
    
    // Transform y to raw Beta scale and compute density
    real y_raw_right = (y - threshold) / (1.0 - threshold);
    // Clamp for stability
    y_raw_right = fmax(eps, fmin(y_raw_right, 1.0 - eps));
    
    // Log-density = log(P(right) * P(continuous|right) * beta_density * jacobian)
    return log(prob_right) + log1m(pex_right) + 
           beta_lpdf(y_raw_right | shape1_right, shape2_right) + 
           log(1.0 / (1.0 - threshold));
  }
}
"
}


#' @rdname rchoco
#' @export
#' @examples
#' \dontrun{
#'   # Requires cmdstanr to be installed and configured
#'   logpdf_func <- choco_lpdf_expose()
#'   # Note: parameter 'p' is now 'mu' in the Stan function
#'   logpdf_func(y = 0.2, mu = 0.6, muleft = 0.3, mudelta = 0.1, phileft = 5,
#'               phidelta = -0.2, pex = 0.1, bex = 0.5, pmid = 0.05)
#'   logpdf_func(y = 0.5, mu = 0.6, muleft = 0.3, mudelta = 0.1, phileft = 5,
#'               phidelta = -0.2, pex = 0.1, bex = 0.5, pmid = 0.05)
#' }
choco_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0("functions {\n", .choco_lpdf(), "\n}")

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
#' @param link_mu,link_muleft,link_mudelta,link_phileft,link_phidelta,link_pex,link_bex,link_pmid Link functions for the parameters.
#' @examples
#' \dontrun{
#' # Example usage in brm formula:
#' # bf(y ~ x1 + (1|group),
#' #    muleft ~ 1,
#' #    mudelta ~ x3,
#' #    phileft ~ 1,
#' #    phidelta ~ 1,
#' #    pex ~ s(age),
#' #    bex ~ 1,
#' #    pmid ~ 1,
#' #    family = choco())
#' }
#' @export
choco <- function(link_mu = "logit", link_muleft = "logit", link_mudelta = "identity",
                  link_phileft = "softplus", link_phidelta = "identity",
                  link_pex = "logit", link_bex = "logit", link_pmid = "logit") {
  brms::custom_family(
    name = "choco",
    dpars = c("mu", "muleft", "mudelta", "phileft", "phidelta", "pex", "bex", "pmid"),
    links = c(link_mu, link_muleft, link_mudelta, link_phileft, link_phidelta, 
              link_pex, link_bex, link_pmid),
    lb = c(0, 0, NA, 0, NA, 0, 0, 0),
    ub = c(1, 1, NA, NA, NA, 1, 1, 1),
    type = "real"  # Continuous outcome variable
  )
}



# brms --------------------------------------------------------------------

#' @rdname rchoco
#' @export
posterior_predict_choco <- function(i, prep, ...) {
  # Extract draws for each parameter for the i-th observation
  mu       <- brms::get_dpar(prep, "mu", i = i) 
  muleft   <- brms::get_dpar(prep, "muleft", i = i)
  mudelta  <- brms::get_dpar(prep, "mudelta", i = i)
  phileft  <- brms::get_dpar(prep, "phileft", i = i)
  phidelta <- brms::get_dpar(prep, "phidelta", i = i)
  pex      <- brms::get_dpar(prep, "pex", i = i)
  bex      <- brms::get_dpar(prep, "bex", i = i)
  pmid     <- brms::get_dpar(prep, "pmid", i = i)

  n_draws <- length(mu) # Changed p to mu

  threshold <- rep(0.5, length.out = n_draws)

  # Simulate using rchoco (vectorized)
  # Note: rchoco uses 'p', so we pass 'mu' draws to the 'p' argument
  final_out <- rchoco(n = n_draws, p = mu, muleft = muleft, mudelta = mudelta,
                       phileft = phileft, phidelta = phidelta, pex = pex,
                       bex = bex, pmid = pmid, threshold = 0.5)

  as.matrix(final_out)
}

#' @rdname rchoco
#' @export
posterior_epred_choco <- function(prep) {
  # Extract draws for parameters (matrices: draws x observations)
  mu       <- brms::get_dpar(prep, "mu") # p
  muleft   <- brms::get_dpar(prep, "muleft")
  mudelta  <- brms::get_dpar(prep, "mudelta")
  pex      <- brms::get_dpar(prep, "pex")
  bex      <- brms::get_dpar(prep, "bex")
  pmid     <- brms::get_dpar(prep, "pmid")

  n_draws <- nrow(mu)
  threshold <- 0.5

  # --- Calculate derived parameters needed for expectation ---
  eps_mu <- 1e-9
  logit_muleft <- log(muleft / (1 - muleft))
  logit_muright <- logit_muleft + mudelta
  muright <- exp(logit_muright) / (1 + exp(logit_muright))
  muleft_c <- pmax(eps_mu, pmin(muleft, 1 - eps_mu))
  muright_c <- pmax(eps_mu, pmin(muright, 1 - eps_mu))
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))

  # --- Calculate components of the expectation ---
  # Probabilities of components
  prob_0 <- (1 - pmid) * (1 - mu) * pex_left 
  prob_1 <- (1 - pmid) * mu * pex_right      
  prob_thresh <- pmid
  prob_left_cont <- (1 - pmid) * (1 - mu) * (1 - pex_left) 
  prob_right_cont <- (1 - pmid) * mu * (1 - pex_right)    

  # Conditional expectations
  mean_left_cont <- threshold * (1 - muleft_c)
  mean_right_cont <- threshold + (1 - threshold) * muright_c

  # Calculate total expectation
  epred <- prob_0 * 0 +
           prob_1 * 1 +
           prob_thresh * threshold +
           prob_left_cont * mean_left_cont +
           prob_right_cont * mean_right_cont

  epred <- pmax(0, pmin(epred, 1))
  epred
}



#' @rdname rchoco
#' @export
log_lik_choco <- function(i, prep) {
  # Extract observed value
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y_scalar <- prep$data$Y[i]

  # Extract model draws
  mu       <- brms::get_dpar(prep, "mu", i = i)
  muleft   <- brms::get_dpar(prep, "muleft", i = i)
  mudelta  <- brms::get_dpar(prep, "mudelta", i = i)
  phileft  <- brms::get_dpar(prep, "phileft", i = i)
  phidelta <- brms::get_dpar(prep, "phidelta", i = i)
  pex      <- brms::get_dpar(prep, "pex", i = i)
  bex      <- brms::get_dpar(prep, "bex", i = i)
  pmid     <- brms::get_dpar(prep, "pmid", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0))

  y_vec <- rep(y_scalar, length.out = n_draws)

  # Calculate log-likelihood using dchoco
  # Note: dchoco uses 'p', so pass 'mu' draws to 'p' argument
  ll <- dchoco(x = y_vec, p = mu, muleft = muleft, mudelta = mudelta, # Changed p = p to p = mu
                phileft = phileft, phidelta = phidelta, pex = pex,
                bex = bex, pmid = pmid, threshold = 0.5, log = TRUE)

  ll[is.nan(ll) | is.na(ll)] <- -Inf
  ll
}