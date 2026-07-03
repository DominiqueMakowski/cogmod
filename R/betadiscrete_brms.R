# Stan Functions and Custom Family for brms (Discrete Beta)

# Stanvars ----------------------------------------------------------------

#' @keywords internal
.betadiscrete_lpmf <- function() {
  "
// Log probability mass function for the Discrete Beta distribution
// (Sciandra et al., 2024, Sect. 3.1)
//   y   : observed rating, integer in {1, ..., k}
//   mu  : mean of the underlying Beta distribution (0 < mu < 1); the
//         'liking' indicator on the logit scale
//   phi : precision of the underlying Beta distribution (alpha + beta > 0);
//         the 'agreement' indicator on the log scale
//   k   : number of rating categories (fixed, passed in as data)
real betadiscrete_lpmf(int y, real mu, real phi, int k) {
  real alpha;
  real beta_par;
  real upper_lcdf;
  real lower_lcdf;

  if (y < 1 || y > k) {
    reject(\"betadiscrete_lpmf: y must be an integer between 1 and k; found y = \", y);
  }

  alpha = mu * phi * 2;
  beta_par = (1 - mu) * phi * 2;

  // P(R = y) = F_B(y/k) - F_B((y-1)/k), computed on the log scale for
  // numerical stability via log_diff_exp(log(upper), log(lower)).
  upper_lcdf = (y == k) ? 0.0 : beta_lcdf(y * 1.0 / k | alpha, beta_par);
  lower_lcdf = (y == 1) ? negative_infinity() : beta_lcdf((y - 1) * 1.0 / k | alpha, beta_par);

  return log_diff_exp(upper_lcdf, lower_lcdf);
}
"
}


#' @rdname rbetadiscrete
#' @examples
#' # You can expose the lpmf function as follows:
#' # betadiscrete_lpmf <- betadiscrete_lpmf_expose()
#' # betadiscrete_lpmf(y = 7, mu = 0.66, phi = 3.51, k = 10)
#'
#' @export
betadiscrete_lpmf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Build the final Stan code string
  stancode <- paste0(
    "functions {\n",
    .betadiscrete_lpmf(),
    "\n}"
  )

  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$betadiscrete_lpmf
}

#' @rdname rbetadiscrete
#' @export
betadiscrete_stanvars <- function() {
  brms::stanvar(scode = .betadiscrete_lpmf(), block = "functions")
}


#' @rdname rbetadiscrete
#' @param link_mu,link_phi Link functions for the parameters.
#'
#' @examples
#' # Fitting with brms
#' # Note: Because `k` is fixed data rather than a distributional parameter, it
#' # must be passed to the model via the brms::vint() addition term
#' # fit <- brms::brm(
#' #   brms::bf(rating | vint(k) ~ predictor),
#' #   data = data,
#' #   family = betadiscrete(),
#' #   stanvars = betadiscrete_stanvars()
#' # )
#'
#' @export
betadiscrete <- function(link_mu = "logit", link_phi = "log") {
  brms::custom_family(
    name = "betadiscrete",
    dpars = c("mu", "phi"),
    links = c(link_mu, link_phi),
    lb = c(NA, 0), # phi > 0
    ub = c(NA, NA),
    type = "int", # discrete outcome
    vars = "vint1[n]" # k, passed via vint()
  )
}

# brms Post-processing Functions -------------------------------------------

#' @rdname rbetadiscrete
#' @export
log_lik_betadiscrete <- function(i, prep) {
  y <- prep$data$Y[i]
  k <- prep$data$vint1[i]

  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(numeric(0))
  }

  y_vec <- rep(y, length.out = n_draws)

  ll <- dbetadiscrete(x = y_vec, mu = mu, phi = phi, k = k, log = TRUE)
  ll[is.nan(ll) | is.na(ll)] <- -Inf
  ll
}


#' @rdname rbetadiscrete
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_betadiscrete <- function(i, prep, ...) {
  k <- prep$data$vint1[i]

  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(matrix(numeric(0), ncol = 1))
  }

  final_out <- rbetadiscrete(n = n_draws, mu = mu, phi = phi, k = k)
  as.matrix(final_out)
}


#' @rdname rbetadiscrete
#' @export
posterior_epred_betadiscrete <- function(prep) {
  mu <- brms::get_dpar(prep, "mu") # draws x obs
  phi <- brms::get_dpar(prep, "phi") # draws x obs
  k <- prep$data$vint1 # length = obs

  n_draws <- nrow(mu)
  n_obs <- ncol(mu)
  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  # E[R | alpha, beta] = k - sum_{j=1}^{k-1} F_B(j/k; alpha, beta)  (eq. 4)
  # Grouping by unique k avoids redundant computation when (as is typical)
  # all observations share the same number of rating categories.
  for (k_val in unique(k)) {
    idx <- which(k == k_val)

    if (k_val <= 1) {
      epred[, idx] <- k_val
      next
    }

    alpha_sub <- mu[, idx, drop = FALSE] * phi[, idx, drop = FALSE] * 2
    beta_sub <- (1 - mu[, idx, drop = FALSE]) * phi[, idx, drop = FALSE] * 2

    cum_sum <- matrix(0, nrow = n_draws, ncol = length(idx))
    for (j in seq_len(k_val - 1)) {
      cum_sum <- cum_sum + stats::pbeta(j / k_val, alpha_sub, beta_sub)
    }

    epred[, idx] <- k_val - cum_sum
  }

  epred
}
