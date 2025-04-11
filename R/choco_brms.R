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

  // --- 1. Parameter Validation ---
  // Check input parameter ranges (essential for preventing downstream errors)
   if (mu < 0.0 || mu > 1.0 || muleft <= 0.0 || muleft >= 1.0 || 
      phileft <= 0.0 || pex < 0.0 || pex > 1.0 || 
      bex < 0.0 || bex > 1.0 || pmid < 0.0 || pmid > 1.0) {
    return negative_infinity();
  }
  // Validate outcome y (essential)
  if (y < 0.0 || y > 1.0) {
    return negative_infinity();
  }

  // --- Case 1: y is exactly at threshold - only need pmid ---
  if (abs(y - threshold) < eps) {
    if (pmid < eps) {
      // print(\"Observation == threshold (0.5), but pmid is 0.0. Relax pmid or nudge mid-values.\");
      return negative_infinity();
    }
    return log(pmid);
  }

  // --- Parameter Clamping ---
  real p_c = fmax(eps, fmin(mu, 1.0 - eps));
  real muleft_c = fmax(eps, fmin(muleft, 1.0 - eps)); 
  real pmid_c = pmid; // Already validated, no need to clamp
  
  // --- Case 2: y is at 0 - only need left component and extreme probability ---
  if (abs(y) < eps) {
    real prob_left = (1.0 - pmid_c) * (1.0 - p_c);
    real pex_left = (1.0 - bex) * (pex * 2.0);
    
    // Clamp extreme probability to [0,1] range
    pex_left = fmin(1.0, fmax(0.0, pex_left));
    
    // If either probability is effectively zero, return -Inf
    if (prob_left < eps || pex_left < eps) return negative_infinity();
    
    // Only need log of extreme probability at 0
    return log(prob_left) + log(pex_left);
  }
  
  // --- Case 3: y is at 1 - only need right component and extreme probability ---
  if (abs(y - 1.0) < eps) {
    real prob_right = (1.0 - pmid_c) * p_c;
    real pex_right = bex * (pex * 2.0);
    
    // Clamp extreme probability to [0,1] range
    pex_right = fmin(1.0, fmax(0.0, pex_right));
    
    // If either probability is effectively zero, return -Inf
    if (prob_right < eps || pex_right < eps) return negative_infinity();
    
    // Only need log of extreme probability at 1
    return log(prob_right) + log(pex_right);
  }
  
  // --- Case 4: y is between 0 and threshold ---
  if (y < threshold) {
    // Only compute left side parameters
    real prob_left = (1.0 - pmid_c) * (1.0 - p_c);
    real pex_left = (1.0 - bex) * (pex * 2.0);

    // Clamp extreme probability to [0,1] range
    pex_left = fmin(1.0, fmax(0.0, pex_left));
    
    // Bail out early if component has zero probability
    if (prob_left < eps || pex_left >= 1.0 - eps) return negative_infinity();
    
    // Calculate Beta parameters for left side only
    real phileft_c = fmax(eps, phileft);
    real precision_left = phileft_c * 2.0;
    real shape1_left = muleft_c * precision_left;
    real shape2_left = (1.0 - muleft_c) * precision_left;
    
    // Transform y to raw Beta scale
    real y_raw_left = 1.0 - y / threshold;
    y_raw_left = fmax(eps, fmin(y_raw_left, 1.0 - eps));
    
    // Compute log density for left side
    return log(prob_left) + log1m(pex_left) +
           beta_lpdf(y_raw_left | shape1_left, shape2_left) +
           log(1.0 / threshold);
  }
  
  // --- Case 5: y is between threshold and 1 ---
  else {
    // Only compute right side parameters
    real prob_right = (1.0 - pmid_c) * p_c;
    real pex_right = bex * (pex * 2.0);

    // Clamp extreme probability to [0,1] range
    pex_right = fmin(1.0, fmax(0.0, pex_right));
    
    // Bail out early if component has zero probability
    if (prob_right < eps || pex_right >= 1.0 - eps) return negative_infinity();
    
    // Calculate right-side parameters only when needed
    real phileft_c = fmax(eps, phileft);
    real phiright = fmax(eps, phileft_c * exp(fmin(10.0, fmax(-10.0, phidelta))));
    
    // Calculate derived parameters for right side
    real logit_muleft = log(muleft_c / (1.0 - muleft_c));
    real logit_muright = logit_muleft + mudelta;
    real muright;
    
    // Stable inverse logit
    if (logit_muright > 15.0) muright = 1.0 - eps;
    else if (logit_muright < -15.0) muright = eps;
    else muright = inv_logit(logit_muright);
    
    // Calculate Beta parameters for right side
    real precision_right = phiright * 2.0;
    real shape1_right = muright * precision_right;
    real shape2_right = (1.0 - muright) * precision_right;
    
    // Transform y to raw Beta scale
    real y_raw_right = (y - threshold) / (1.0 - threshold);
    y_raw_right = fmax(eps, fmin(y_raw_right, 1.0 - eps));
    
    // Compute log density for right side
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



