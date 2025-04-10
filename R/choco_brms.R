# Stan Functions and Custom Family for brms (CHOCO)

#' @keywords internal
.choco_lpdf <- function() {
"
// Log probability density function for the CHOCO distribution (Optimized)
// Accounts for point masses at 0, threshold (fixed at 0.5), and 1, and continuous Beta components.
real choco_lpdf(real y, real mu, real muleft, real mudelta, real phileft,
                 real phidelta, real pex, real bex, real pmid) {

  real eps = 1e-9; // Tolerance for floating point comparisons
  real threshold = 0.5; // Hardcoded threshold

  // --- 1. Parameter Validation ---
  // Check input parameter ranges (essential for preventing downstream errors)
  if (mu < 0.0 || mu > 1.0 ||
      muleft <= 0.0 || muleft >= 1.0 || // muleft must be strictly > 0 and < 1
      phileft <= 0.0 || // phileft must be strictly > 0
      pex < 0.0 || pex > 1.0 ||
      bex < 0.0 || bex > 1.0 ||
      pmid < 0.0 || pmid > 1.0) {
    return negative_infinity();
  }
  // Validate outcome y
  if (y < 0.0 || y > 1.0) {
     return negative_infinity();
  }

  // --- 2. Derived Parameter Computation ---
  // Calculate right-side parameters and clamp means for beta shapes
  real logit_muleft = log(muleft / (1.0 - muleft));
  real logit_muright = logit_muleft + mudelta;
  real muright = inv_logit(logit_muright);
  real muleft_c = fmax(eps, fmin(muleft, 1.0 - eps)); // Clamped version for shape calculation
  real muright_c = fmax(eps, fmin(muright, 1.0 - eps)); // Clamped version for shape calculation

  real phiright = phileft * exp(phidelta);
  // Check derived precision parameter
  if (phiright <= 0.0) {
      return negative_infinity();
  }

  // Calculate effective extreme probabilities for left (0) and right (1) components
  real pex_left = fmin(1.0, fmax(0.0, (1.0 - bex) * (pex * 2.0)));
  real pex_right = fmin(1.0, fmax(0.0, bex * (pex * 2.0)));

  // Calculate Beta distribution shapes
  real precision_left = phileft * 2.0;
  real precision_right = phiright * 2.0;
  real shape1_left = muleft_c * precision_left;
  real shape2_left = (1.0 - muleft_c) * precision_left;
  real shape1_right = muright_c * precision_right;
  real shape2_right = (1.0 - muright_c) * precision_right;

  // --- 3. Pre-calculate Log Probabilities of Components ---
  // Calculate logs needed for mixture components, handling edge cases (log(0) = -Inf)
  real log_mu = (mu <= eps) ? negative_infinity() : log(mu);
  real log1m_mu = (mu >= 1.0 - eps) ? negative_infinity() : log1m(mu);
  real log_pmid = (pmid <= eps) ? negative_infinity() : log(pmid);
  real log1m_pmid = (pmid >= 1.0 - eps) ? negative_infinity() : log1m(pmid);
  real log_pex_left = (pex_left <= eps) ? negative_infinity() : log(pex_left);
  real log1m_pex_left = (pex_left >= 1.0 - eps) ? negative_infinity() : log1m(pex_left);
  real log_pex_right = (pex_right <= eps) ? negative_infinity() : log(pex_right);
  real log1m_pex_right = (pex_right >= 1.0 - eps) ? negative_infinity() : log1m(pex_right);

  // Pre-calculate log jacobians for transformations
  real log_jac_left = -log(threshold); // log(1/threshold)
  real log_jac_right = -log1m(threshold); // log(1/(1-threshold))

  // --- 4. Calculate Log Probability/Density based on y ---
  real log_lik;

  if (abs(y - 0.0) < eps) { // Case 1: y = 0 (Point Mass)
    // Need (1-pmid) > 0, (1-mu) > 0, pex_left > 0
    if (log1m_pmid == negative_infinity() || log1m_mu == negative_infinity() || log_pex_left == negative_infinity()) {
      log_lik = negative_infinity();
    } else {
      // log( (1-pmid)*(1-mu)*pex_left )
      log_lik = log1m_pmid + log1m_mu + log_pex_left;
    }
  } else if (abs(y - 1.0) < eps) { // Case 2: y = 1 (Point Mass)
    // Need (1-pmid) > 0, mu > 0, pex_right > 0
    if (log1m_pmid == negative_infinity() || log_mu == negative_infinity() || log_pex_right == negative_infinity()) {
      log_lik = negative_infinity();
    } else {
      // log( (1-pmid)*mu*pex_right )
      log_lik = log1m_pmid + log_mu + log_pex_right;
    }
  } else if (abs(y - threshold) < eps) { // Case 3: y = threshold (Point Mass)
    // Need pmid > 0. log_pmid will be -Inf if pmid was <= eps.
    log_lik = log_pmid;
  } else if (y > eps && y < threshold - eps) { // Case 4: 0 < y < threshold (Left Continuous Part)
    // Need (1-pmid) > 0, (1-mu) > 0, (1-pex_left) > 0
    if (log1m_pmid == negative_infinity() || log1m_mu == negative_infinity() || log1m_pex_left == negative_infinity()) {
      log_lik = negative_infinity();
    } else {
      // Transform y back to the [0, 1] scale of the underlying Beta
      real y_raw_left = 1.0 - y / threshold;
      // Clamp input to beta_lpdf for stability, although y range check helps
      y_raw_left = fmax(eps, fmin(y_raw_left, 1.0 - eps));
      real log_beta_dens = beta_lpdf(y_raw_left | shape1_left, shape2_left);

      // Check if beta density calculation failed (e.g., invalid shapes, though unlikely with checks)
      if (log_beta_dens == negative_infinity()) {
         log_lik = negative_infinity();
      } else {
         // log( (1-pmid)*(1-mu)*(1-pex_left) * beta_dens(y_raw|...) * jacobian )
         log_lik = log1m_pmid + log1m_mu + log1m_pex_left + log_beta_dens + log_jac_left;
      }
    }
  } else if (y > threshold + eps && y < 1.0 - eps) { // Case 5: threshold < y < 1 (Right Continuous Part)
    // Need (1-pmid) > 0, mu > 0, (1-pex_right) > 0
    if (log1m_pmid == negative_infinity() || log_mu == negative_infinity() || log1m_pex_right == negative_infinity()) {
      log_lik = negative_infinity();
    } else {
      // Transform y back to the [0, 1] scale of the underlying Beta
      real y_raw_right = (y - threshold) / (1.0 - threshold);
      // Clamp input to beta_lpdf for stability
      y_raw_right = fmax(eps, fmin(y_raw_right, 1.0 - eps));
      real log_beta_dens = beta_lpdf(y_raw_right | shape1_right, shape2_right);

      // Check if beta density calculation failed
      if (log_beta_dens == negative_infinity()) {
         log_lik = negative_infinity();
      } else {
        // log( (1-pmid)*mu*(1-pex_right) * beta_dens(y_raw|...) * jacobian )
        log_lik = log1m_pmid + log_mu + log1m_pex_right + log_beta_dens + log_jac_right;
      }
    }
  } else { // Should not be reached due to initial y validation, but acts as a safeguard
    log_lik = negative_infinity();
  }

  return log_lik;
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
    # lb = c(0, 0, NA, 0, NA, 0, 0, 0),
    # ub = c(1, 1, NA, NA, NA, 1, 1, 1),
    lb = c(NA, NA, NA, 0, NA, NA, NA, NA),
    ub = c(NA, NA, NA, NA, NA, NA, NA, NA),
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