# Stan Functions and Custom Family for brms (CHOCO Model)

#' @keywords internal
.choco2_lpdf <- function() {
"
// mu is the probability p of the right side (0 < mu < 1) (named 'mu' as per Stan convention for main parameters)
real choco2_lpdf(
    real y,
    real mu,          // P(right | not-mid)
    real confright,   // mu for right Beta-Gate
    real precright,   // phi for right Beta-Gate
    real confleft,    // mu_left = 1 - confleft
    real precleft,    // phi for left Beta-Gate
    real pex,         // proportion-extreme total
    real bex,         // balance of extremes (→ right vs left)
    real pmid         // P(mid)
) {
    // Hardcoded threshold 
    real threshold = 0.5;
    real eps = 1e-10;

    // 1) Domain checks
    if (
        y < 0 || y > 1 ||
        mu < 0 || mu > 1 ||
        pmid < 0 || pmid > 1 ||
        confright <= 0 || confright >= 1 ||
        precright <= 0 ||
        confleft <= 0 || confleft >= 1 ||
        precleft <= 0 ||
        pex < 0 || pex > 1 ||
        bex < 0 || bex > 1
    ) {
        return negative_infinity();
    }

    // 2) Mixture weights
    real p_not_mid = 1 - pmid;
    real p_left = p_not_mid * (1 - mu);
    real p_right = p_not_mid * mu;

    // 3) Point mass at threshold
    if (abs(y - threshold) < eps) {
        return (pmid > 0) ? log(pmid) : negative_infinity();
    }

    // 4) Left segment: y in [0, threshold)
    if (y < threshold) {
        if (p_left <= 0) return negative_infinity();

        // Rescale to [0,1]
        real y_rescaled = y / threshold;

        // Beta-Gate parameters for the left side
        real mu_left = 1 - confleft;
        real phi_left = precleft;
        real pex_left = pex * (1 - bex);
        real bex_left = 0;

        // --- Inlined betagate_lpdf for left side (bex=0) ---
        real cutzerolog = logit(pex_left * (1.0 - bex_left)); // = logit(pex_left)
        real cutonelog = logit(1.0 - pex_left * bex_left);    // = logit(1.0)
        real mu_ql = logit(mu_left);
        real log_density_scaled;
        if (abs(y_rescaled - 0.0) < eps) {
            log_density_scaled = log1m_inv_logit(mu_ql - cutzerolog);
        } else if (abs(y_rescaled - 1.0) < eps) {
            log_density_scaled = log_inv_logit(mu_ql - cutonelog);
        } else {
            real log_prob_middle = log_inv_logit(mu_ql - cutzerolog);
            real shape1 = mu_left * phi_left * 2.0;
            real shape2 = (1.0 - mu_left) * phi_left * 2.0;
            real log_beta_dens = beta_lpdf(y_rescaled | shape1, shape2);
            log_density_scaled = log_prob_middle + log_beta_dens;
        }

        // Jacobian for y → y_rescaled (skip at exact 0)
        real log_jacobian = (y_rescaled > eps) ? -log(threshold) : 0;

        return log(p_left) + log_density_scaled + log_jacobian;
    }

    // 5) Right segment: y in (threshold, 1]
    if (y > threshold) {
        if (p_right <= 0) return negative_infinity();

        // Rescale to [0,1]
        real y_rescaled = (y - threshold) / (1 - threshold);

        // Beta-Gate parameters for the right side
        real mu_right = confright;
        real phi_right = precright;
        real pex_right = pex * bex;
        real bex_right = 1;

        // --- Inlined betagate_lpdf for right side (bex=1) ---
        real cutzerolog = logit(pex_right * (1.0 - bex_right)); // = logit(0)
        real cutonelog = logit(1.0 - pex_right * bex_right);    // = logit(1.0 - pex_right)
        real mu_ql = logit(mu_right);
        real log_density_scaled;
        if (abs(y_rescaled - 0.0) < eps) {
            log_density_scaled = log1m_inv_logit(mu_ql - cutzerolog);
        } else if (abs(y_rescaled - 1.0) < eps) {
            log_density_scaled = log_inv_logit(mu_ql - cutonelog);
        } else {
            real log_prob_middle = log_diff_exp(
                log_inv_logit(mu_ql - cutzerolog),
                log_inv_logit(mu_ql - cutonelog)
            );
            real shape1 = mu_right * phi_right * 2.0;
            real shape2 = (1.0 - mu_right) * phi_right * 2.0;
            real log_beta_dens = beta_lpdf(y_rescaled | shape1, shape2);
            log_density_scaled = log_prob_middle + log_beta_dens;
        }

        // Jacobian for y → y_rescaled (skip at exact 1)
        real log_jacobian = (y_rescaled < 1 - eps) ? -log1m(threshold) : 0;

        return log(p_right) + log_density_scaled + log_jacobian;
    }

    // If no branch matched, return -Inf
    return negative_infinity();
}
"
}

#' @rdname rchoco2
#' @export
choco2_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0(
    "functions {\n",
    .choco2_lpdf(),
    "\n}" )
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$choco2_lpdf
}



#' @rdname rchoco2
#' @export
choco2_stanvars <- function() {
  brms::stanvar(scode = .choco2_lpdf(), block = "functions")
}



#' @rdname rchoco2
#' @param link_mu,link_confright,link_precright,link_confleft,link_precleft,link_pex,link_bex,link_pmid Link functions for the parameters.
#' @examples
#' \dontrun{
#' # Example usage in brm formula:
#' # bf(y ~ x1 + (1|group),
#' #    confright ~ x3,
#' #    confleft ~ x3,
#' #    precright ~ 1,
#' #    precleft ~ 1,
#' #    pex ~ s(age),
#' #    bex ~ 1,
#' #    pmid ~ 1,
#' #    family = choco2())
#' }
#' @export
choco2 <- function(
  link_mu = "logit", link_confright = "logit", link_precright = "softplus",
  link_confleft = "logit", link_precleft = "softplus",
  link_pex = "logit", link_bex = "logit",
  link_pmid = "logit"
) {
  brms::custom_family(
    name  = "choco2",
    dpars = c("mu","confright","precright","confleft",
              "precleft","pex","bex","pmid"),
    links = c(link_mu, link_confright, link_precright,
              link_confleft, link_precleft, link_pex, link_bex,
              link_pmid),
    lb = c(0,0,0,0,0,0,0,0),
    ub = c(1,1,NA,1,NA,1,1,1),
    type = "real"
  )
}


# brms Post-processing Functions -------------------------------------------


#' @rdname rchoco2
#' @export
log_lik_choco2 <- function(i, prep) {
  # Extract observed value for the i-th observation
  if (!"Y" %in% names(prep$data)) stop("Outcome variable 'Y' not found in prep$data.")
  y_scalar <- prep$data$Y[i]

  # Extract model draws (vectors) for the i-th observation
  mu         <- brms::get_dpar(prep, "mu", i = i)
  confright  <- brms::get_dpar(prep, "confright", i = i)
  precright  <- brms::get_dpar(prep, "precright", i = i)
  confleft   <- brms::get_dpar(prep, "confleft", i = i)
  precleft   <- brms::get_dpar(prep, "precleft", i = i)
  pex        <- brms::get_dpar(prep, "pex", i = i)
  bex        <- brms::get_dpar(prep, "bex", i = i)
  pmid       <- brms::get_dpar(prep, "pmid", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) return(numeric(0))

  y_vec <- rep(y_scalar, length.out = n_draws)

  # Calculate log-likelihood using the vectorized dchoco2 function
  ll <- dchoco2(
    x = y_vec,
    p = mu,
    confright = confright,
    precright = precright,
    confleft = confleft,
    precleft = precleft,
    pex = pex,
    bex = bex,
    pmid = pmid,
    threshold = 0.5,
    log = TRUE
  )

  ll[is.nan(ll) | is.na(ll)] <- -Inf
  ll
}


#' @rdname rchoco2
#' @export
posterior_predict_choco2 <- function(i, prep, ...) {
  # Extract draws for each parameter for the i-th observation
  # brms::get_dpar should return vectors respecting prep$ndraws
  mu         <- brms::get_dpar(prep, "mu", i = i)
  confright  <- brms::get_dpar(prep, "confright", i = i)
  precright  <- brms::get_dpar(prep, "precright", i = i)
  confleft   <- brms::get_dpar(prep, "confleft", i = i)
  precleft   <- brms::get_dpar(prep, "precleft", i = i)
  pex        <- brms::get_dpar(prep, "pex", i = i)
  bex        <- brms::get_dpar(prep, "bex", i = i)
  pmid       <- brms::get_dpar(prep, "pmid", i = i)

  # Determine number of draws from the length of the parameter vectors
  # This length should correspond to prep$ndraws if ndraws was specified
  n_draws <- length(mu)
  if (n_draws == 0) {
      # If get_dpar returns zero length, return an empty matrix.
      # The number of rows should ideally match prep$ndraws, but 0 is safer if prep$ndraws is unreliable.
      return(matrix(numeric(0), nrow = 0, ncol = 1))
  }

  # Generate exactly n_draws predictions using the extracted parameters
  final_out <- rchoco2(
    n = n_draws, # Use n_draws determined from get_dpar length
    p = mu,
    confright = confright,
    precright = precright,
    confleft = confleft,
    precleft = precleft,
    pex = pex,
    bex = bex,
    pmid = pmid,
    threshold = 0.5
  )

  # Return as a single-column matrix, as required by brms
  as.matrix(final_out)
}


#' @rdname rchoco2
#' @export
posterior_epred_choco2 <- function(prep) {
  # Fetch parameters - these might be matrices (draws x obs) or vectors (draws)
  mu         <- brms::get_dpar(prep, "mu")
  confright  <- brms::get_dpar(prep, "confright")
  confleft   <- brms::get_dpar(prep, "confleft")
  pex        <- brms::get_dpar(prep, "pex")
  bex        <- brms::get_dpar(prep, "bex")
  pmid       <- brms::get_dpar(prep, "pmid")
  # Note: precright and precleft are not needed for the expected value calculation

  # Determine dimensions
  # Use mu as reference; assume others conform or are vectors
  if (is.matrix(mu)) {
    n_draws <- nrow(mu)
    n_obs <- ncol(mu)
  } else { # Assume mu is a vector if not a matrix
    n_draws <- length(mu)
    n_obs <- prep$data$nobs # Get number of observations from prep data
    # If mu is vector, reshape it to matrix for consistent indexing later
    mu <- matrix(mu, nrow = n_draws, ncol = n_obs)
  }

  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  # Helper function to get parameter values for observation j
  # Handles both matrix and vector cases
  get_param_j <- function(param, j) {
    if (is.matrix(param)) {
      return(param[, j, drop = FALSE]) # Use drop=FALSE to keep matrix structure if n_draws=1
    } else {
      return(param) # If vector, return the whole vector
    }
  }

  for (j in seq_len(n_obs)) {
    # Get parameters for observation j using the helper
    mu_j <- get_param_j(mu, j)
    confright_j <- get_param_j(confright, j)
    confleft_j <- get_param_j(confleft, j)
    pex_j <- get_param_j(pex, j)
    bex_j <- get_param_j(bex, j)
    pmid_j <- get_param_j(pmid, j)

    # Side probabilities
    p_not_mid <- 1 - pmid_j
    p_left  <- p_not_mid * (1 - mu_j)
    p_right <- p_not_mid * mu_j
    prob_mid   <- pmid_j

    # Underlying Beta-Gate params means and extreme probs
    # Note: mu_left is 1 - confleft
    mu_right_latent  <- confright_j
    mu_left_latent   <- 1 - confleft_j
    pex_right <- pex_j * bex_j
    pex_left  <- pex_j * (1 - bex_j)

    # Expected value for latent left and right Beta-Gate components
    # E[BetaGate(mu, phi, pex, bex=0)] = (1 - pex) * mu
    # E[BetaGate(mu, phi, pex, bex=1)] = (1 - pex) * mu + pex * 1
    e_latent_left <- (1 - pex_left) * mu_left_latent
    e_latent_right <- (1 - pex_right) * mu_right_latent + pex_right

    # Rescale expected values from latent [0,1] to observed scale
    # Left side: [0, threshold] -> threshold * E_latent
    # Right side: [threshold, 1] -> threshold + (1 - threshold) * E_latent
    threshold <- 0.5 # Assuming fixed threshold
    e_left_rescaled <- threshold * e_latent_left
    e_right_rescaled <- threshold + (1 - threshold) * e_latent_right

    # Combine weighted expected values
    epred[, j] <- p_left * e_left_rescaled +
                  p_right * e_right_rescaled +
                  prob_mid * threshold # E[mid] is just the threshold value
  }

  epred
}