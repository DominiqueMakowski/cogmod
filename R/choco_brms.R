# Stan Functions and Custom Family for brms (CHOCO)

#' @keywords internal
.choco_lpdf <- function() {
"
// Log probability density function for the CHOCO distribution (Optimized)
// Accounts for point masses at 0, threshold (fixed at 0.5), and 1, and continuous Beta components.
// Uses conf, confleft, prec, precleft parameterization.
real choco_lpdf(real y, real mu, real conf, real confleft, real prec,
                 real precleft, real pex, real bex, real pmid) {
  real eps = 1e-10;   // Tolerance for floating point comparisons
  real threshold = 0.5; // Hardcoded threshold for the confidence split

  // --- Parameter Validation ---
  // mu in [0,1]; conf in (0,1); prec > 0; pex, bex, pmid in [0,1]
  if (mu < 0.0 || mu > 1.0 ||
      conf <= 0.0 || conf >= 1.0 ||
      prec <= 0.0 ||
      pex < 0.0 || pex > 1.0 ||
      bex < 0.0 || bex > 1.0 ||
      pmid < 0.0 || pmid > 1.0) {
    return negative_infinity();
  }
  // Validate outcome y in [0,1]
  if (y < 0.0 || y > 1.0) return negative_infinity();

  // --- Case 1: y is exactly at threshold ---
  if (abs(y - threshold) < eps) {
    if (pmid < eps)
      return negative_infinity();
    return log(pmid);
  }

  // --- Parameter Clamping ---
  real p_c = fmax(eps, fmin(mu, 1.0 - eps));
  real conf_c = fmax(eps, fmin(conf, 1.0 - eps));   // Clamp conf
  real prec_c = fmax(eps, prec);                     // Clamp prec; lower bound eps
  real pmid_c = pmid; // pmid already verified to be in [0,1]

  // --- Calculate Derived Parameters ---
  // Right side parameters are simpler.
  real muright = conf_c;
  real phiright = prec_c;

  // Left side parameters: derive muleft via a logit transformation.
  real logit_conf = log(conf_c / (1.0 - conf_c));
  real logit_muleft = logit_conf + confleft;
  real muleft;
  if (logit_muleft > 15.0) 
    muleft = 1.0 - eps;
  else if (logit_muleft < -15.0) 
    muleft = eps;
  else 
    muleft = inv_logit(logit_muleft);
  muleft = fmax(eps, fmin(muleft, 1.0 - eps));

  real phileft = prec_c * exp(fmin(10.0, fmax(-10.0, precleft)));
  phileft = fmax(eps, phileft);

  // Effective extreme probabilities (for left/right endpoints)
  real pex_left = fmin(1.0, fmax(0.0, (1.0 - bex) * (pex * 2.0)));
  real pex_right = fmin(1.0, fmax(0.0, bex * (pex * 2.0)));

  // --- Case 2: y is 0 (extreme on left) ---
  if (abs(y) < eps) {
    real prob_left = (1.0 - pmid_c) * (1.0 - p_c);
    if (prob_left < eps || pex_left < eps)
      return negative_infinity();
    return log(prob_left) + log(pex_left);
  }

  // --- Case 3: y is 1 (extreme on right) ---
  if (abs(y - 1.0) < eps) {
    real prob_right = (1.0 - pmid_c) * p_c;
    if (prob_right < eps || pex_right < eps)
      return negative_infinity();
    return log(prob_right) + log(pex_right);
  }

  // --- Case 4: y is between 0 and threshold (left side) ---
  if (y < threshold) {
    real prob_left = (1.0 - pmid_c) * (1.0 - p_c);
    if (prob_left < eps || pex_left >= 1.0 - eps)
      return negative_infinity();

    real precision_left = phileft * 2.0;
    real shape1_left = muleft * precision_left;
    real shape2_left = (1.0 - muleft) * precision_left;

    // Transform y to a Beta–scale variable: y_raw = 1 - (y / threshold)
    real y_raw_left = 1.0 - y / threshold;
    y_raw_left = fmax(eps, fmin(y_raw_left, 1.0 - eps));

    return log(prob_left) + log1m(pex_left) +
           beta_lpdf(y_raw_left | shape1_left, shape2_left) +
           log(1.0 / threshold); // Jacobian adjustment
  }
  // --- Case 5: y is between threshold and 1 (right side) ---
  else {
    real prob_right = (1.0 - pmid_c) * p_c;
    if (prob_right < eps || pex_right >= 1.0 - eps)
      return negative_infinity();

    real precision_right = phiright * 2.0;
    real shape1_right = muright * precision_right;
    real shape2_right = (1.0 - muright) * precision_right;

    // Transform y to Beta scale: y_raw = (y - threshold) / (1.0 - threshold)
    real y_raw_right = (y - threshold) / (1.0 - threshold);
    y_raw_right = fmax(eps, fmin(y_raw_right, 1.0 - eps));

    return log(prob_right) + log1m(pex_right) +
           beta_lpdf(y_raw_right | shape1_right, shape2_right) +
           log(1.0 / (1.0 - threshold)); // Jacobian adjustment
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
#'   logpdf_func(y = 0.2, mu = 0.6, conf = 0.3, confleft = 0.1, prec = 5,
#'               precleft = -0.2, pex = 0.1, bex = 0.5, pmid = 0.05)
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
#' @param link_mu,link_conf,link_confleft,link_prec,link_precleft,link_pex,link_bex,link_pmid Link functions for the parameters.
#' @examples
#' \dontrun{
#' # Example usage in brm formula:
#' # bf(y ~ x1 + (1|group),
#' #    conf ~ 1,
#' #    confleft ~ x3,
#' #    prec ~ 1,
#' #    precleft ~ 1,
#' #    pex ~ s(age),
#' #    bex ~ 1,
#' #    pmid ~ 1,
#' #    family = choco())
#' }
#' @export
choco <- function(link_mu = "logit", link_conf = "logit", link_confleft = "identity",
                  link_prec = "softplus", link_precleft = "identity",
                  link_pex = "logit", link_bex = "logit", link_pmid = "logit") {
  brms::custom_family(
    name = "choco",
    dpars = c("mu", "conf", "confleft", "prec", "precleft", "pex", "bex", "pmid"), # Updated dpars
    links = c(link_mu, link_conf, link_confleft, link_prec, link_precleft, # Updated links
              link_pex, link_bex, link_pmid),
    lb = c(0, 0, NA, 0, NA, 0, 0, 0), # Updated lb (conf > 0, prec > 0)
    ub = c(1, 1, NA, NA, NA, 1, 1, 1), # Updated ub (conf < 1)
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

  # Extract model draws using new parameter names
  mu       <- brms::get_dpar(prep, "mu", i = i)
  conf     <- brms::get_dpar(prep, "conf", i = i)
  confleft <- brms::get_dpar(prep, "confleft", i = i)
  prec     <- brms::get_dpar(prep, "prec", i = i)
  precleft <- brms::get_dpar(prep, "precleft", i = i)
  pex      <- brms::get_dpar(prep, "pex", i = i)
  bex      <- brms::get_dpar(prep, "bex", i = i)
  pmid     <- brms::get_dpar(prep, "pmid", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0))

  y_vec <- rep(y_scalar, length.out = n_draws)

  # Calculate log-likelihood using dchoco with new parameter names
  # Note: dchoco uses 'p', so pass 'mu' draws to 'p' argument
  ll <- dchoco(x = y_vec, p = mu, conf = conf, confleft = confleft, # Updated parameters
                prec = prec, precleft = precleft, pex = pex,
                bex = bex, pmid = pmid, threshold = 0.5, log = TRUE)

  ll[is.nan(ll) | is.na(ll)] <- -Inf
  ll
}



#' @rdname rchoco
#' @export
posterior_predict_choco <- function(i, prep, ...) {
  # Extract draws for each parameter for the i-th observation using new names
  mu       <- brms::get_dpar(prep, "mu", i = i)
  conf     <- brms::get_dpar(prep, "conf", i = i)
  confleft <- brms::get_dpar(prep, "confleft", i = i)
  prec     <- brms::get_dpar(prep, "prec", i = i)
  precleft <- brms::get_dpar(prep, "precleft", i = i)
  pex      <- brms::get_dpar(prep, "pex", i = i)
  bex      <- brms::get_dpar(prep, "bex", i = i)
  pmid     <- brms::get_dpar(prep, "pmid", i = i)

  n_draws <- length(mu)

  # Simulate using rchoco (vectorized) with new parameter names
  # Note: rchoco uses 'p', so we pass 'mu' draws to the 'p' argument
  final_out <- rchoco(n = n_draws, p = mu, conf = conf, confleft = confleft, # Updated parameters
                       prec = prec, precleft = precleft, pex = pex,
                       bex = bex, pmid = pmid, threshold = 0.5)

  as.matrix(final_out)
}






#' @rdname rchoco
#' @export
posterior_epred_choco <- function(prep) {
  # Extract draws for parameters (matrices: draws x observations) using new names
  mu       <- brms::get_dpar(prep, "mu") # p
  conf     <- brms::get_dpar(prep, "conf")
  confleft <- brms::get_dpar(prep, "confleft")
  prec     <- brms::get_dpar(prep, "prec")
  precleft <- brms::get_dpar(prep, "precleft")
  pex      <- brms::get_dpar(prep, "pex")
  bex      <- brms::get_dpar(prep, "bex")
  pmid     <- brms::get_dpar(prep, "pmid")

  threshold <- 0.5
  eps_mu <- 1e-9 # Tolerance for clamping

  # --- Calculate derived parameters needed for expectation using new logic ---
  # Clamp input parameters first
  conf_c <- pmax(eps_mu, pmin(conf, 1 - eps_mu))
  prec_c <- pmax(eps_mu, prec)

  # Right side parameters
  muright <- conf_c
  # phiright <- prec_c # Not directly needed for mean

  # Left side parameters
  logit_conf <- log(conf_c / (1 - conf_c))
  logit_muleft <- logit_conf + confleft
  muleft <- exp(logit_muleft) / (1 + exp(logit_muleft))
  muleft_c <- pmax(eps_mu, pmin(muleft, 1 - eps_mu)) # Clamp derived muleft
  # phileft <- prec_c * exp(precleft) # Not directly needed for mean

  # Effective extreme probabilities
  pex_left <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
  pex_right <- pmin(1, pmax(0, bex * (pex * 2)))

  # --- Calculate components of the expectation ---
  # Probabilities of components
  prob_0 <- (1 - pmid) * (1 - mu) * pex_left
  prob_1 <- (1 - pmid) * mu * pex_right
  prob_thresh <- pmid
  prob_left_cont <- (1 - pmid) * (1 - mu) * (1 - pex_left)
  prob_right_cont <- (1 - pmid) * mu * (1 - pex_right)

  # Conditional expectations of the continuous parts (based on underlying Beta means)
  # E[X | Left, Continuous] = E[(1 - Y_raw_left) * threshold] = threshold * (1 - E[Y_raw_left]) = threshold * (1 - muleft_c)
  mean_left_cont <- threshold * (1 - muleft_c)
  # E[X | Right, Continuous] = E[threshold + Y_raw_right * (1 - threshold)] = threshold + (1 - threshold) * E[Y_raw_right] = threshold + (1 - threshold) * muright
  mean_right_cont <- threshold + (1 - threshold) * muright # Use muright (which is conf_c)

  # Calculate total expectation: Sum of (Probability * Conditional Expectation) for each component
  epred <- prob_0 * 0 +
           prob_1 * 1 +
           prob_thresh * threshold +
           prob_left_cont * mean_left_cont +
           prob_right_cont * mean_right_cont

  # Ensure the final expectation is within [0, 1]
  epred <- pmax(0, pmin(epred, 1))
  epred
}