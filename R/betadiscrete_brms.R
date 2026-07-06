# Stan Functions and Custom Family for brms (Discrete Beta)

# Stanvars ----------------------------------------------------------------

#' @keywords internal
.betadiscrete_lpmf <- function() {
  "
// Log probability mass function for the (hurdle) Discrete Beta distribution
// (Sciandra et al., 2024, Sect. 3.1)
//   y     : observed rating, integer in {0, 1, ..., k}. y = 0 is only valid
//           when pzero > 0 (hurdle point mass below the 1..k rating scale)
//   mu    : mean of the underlying Beta distribution (0 < mu < 1); the
//           'liking' indicator on the logit scale
//   phi   : precision of the underlying Beta distribution (alpha + beta > 0);
//           the 'agreement' indicator on the log scale
//   pzero : probability of the hurdle point mass at 0 (0 <= pzero < 1)
//   k     : number of rating categories (fixed, passed in as data)
real betadiscrete_lpmf(int y, real mu, real phi, real pzero, int k) {
  real alpha;
  real beta_par;
  real upper_lcdf;
  real lower_lcdf;

  if (y < 0 || y > k) {
    reject(\"betadiscrete_lpmf: y must be an integer between 0 and k; found y = \", y);
  }

  if (y == 0) {
    return log(pzero);
  }

  alpha = mu * phi * 2;
  beta_par = (1 - mu) * phi * 2;

  // P(R = y) = F_B(y/k) - F_B((y-1)/k), computed on the log scale for
  // numerical stability via log_diff_exp(log(upper), log(lower)).
  upper_lcdf = (y == k) ? 0.0 : beta_lcdf(y * 1.0 / k | alpha, beta_par);
  lower_lcdf = (y == 1) ? negative_infinity() : beta_lcdf((y - 1) * 1.0 / k | alpha, beta_par);

  return log1m(pzero) + log_diff_exp(upper_lcdf, lower_lcdf);
}
"
}


#' @rdname rbetadiscrete
#' @examples
#' # You can expose the lpmf function as follows:
#' # betadiscrete_lpmf <- betadiscrete_lpmf_expose()
#' # betadiscrete_lpmf(y = 7, mu = 0.66, phi = 3.51, pzero = 0, k = 10)
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
#' @param link_mu,link_phi,link_pzero Link functions for the parameters. `pzero`
#'   defaults to a `"logit"` link. By default (i.e., if `pzero` is not included
#'   in the `brms::bf()` formula), it is estimated as a single, intercept-only
#'   value shared across all observations (as is done for `pmid` in
#'   [choco()]); it can instead be given predictors to let it vary (`pzero ~ x`),
#'   or fixed to a constant -- e.g., `pzero = 0`, recovering the pure Discrete
#'   Beta model -- directly in `brms::bf()` (as is done for `pmid` in [choco()]).
#'
#' @details
#' Note that `y = 0` is always handled by `pzero` alone, and `k` always refers to the
#' number of categories of the *non-zero* `1:k` part of the scale. What
#' *does* require some care is deciding what `k` should be and whether to
#' estimate or fix `pzero`, depending on how the zero in your data arose:
#' - Scale is `0:N` and 0 is *not* a hurdle (just the lowest of `N + 1`
#'   ordinary ordinal categories, e.g., a 0-10 rating scale with no excess of
#'   zeros): recode the data to `1:(N + 1)` (add 1 to every response), use
#'   `vint(N + 1)`, and fix `pzero = 0` as shown below.
#' - Scale is `0:N` and 0 *is* a hurdle (e.g., a mix of a genuine/excess "zero"
#'   response with an ordinal `1:N` scale): keep the data as-is, use `vint(N)`
#'   (i.e., the total number of categories *minus* the hurdle category), and
#'   let `pzero` be estimated (optionally with predictors, `pzero ~ x`).
#' - Scale is `1:N` with excess responses piling up at the low end (e.g., a
#'   floor effect at the lowest category): recode by subtracting 1
#'   (`1:N` becomes `0:(N - 1)`), then proceed as in the previous bullet,
#'   i.e., `vint(N - 1)` and estimate `pzero`.
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
#' # To also model the hurdle probability (e.g., proportion of zero ratings):
#' # fit <- brms::brm(
#' #   brms::bf(rating | vint(k) ~ predictor, pzero ~ predictor),
#' #   data = data,
#' #   family = betadiscrete(),
#' #   stanvars = betadiscrete_stanvars()
#' # )
#'
#' # To fix pzero at exactly 0, e.g. because your scale has no hurdle:
#' # fit <- brms::brm(
#' #   brms::bf(rating | vint(k) ~ predictor, pzero = 0),
#' #   data = data,
#' #   family = betadiscrete(),
#' #   stanvars = betadiscrete_stanvars()
#' # )
#'
#' @export
betadiscrete <- function(
  link_mu = "logit",
  link_phi = "log",
  link_pzero = "logit"
) {
  brms::custom_family(
    name = "betadiscrete",
    dpars = c("mu", "phi", "pzero"),
    links = c(link_mu, link_phi, link_pzero),
    lb = c(NA, 0, 0), # phi > 0; 0 <= pzero <= 1
    ub = c(NA, NA, 1),
    type = "int", # discrete outcome
    vars = "vint1[n]" # k, passed via vint()
  )
}

# brms Post-processing Functions -------------------------------------------

# Subsets a dpar that may be a draws x obs matrix (estimated) or a
# draws-length vector (fixed, e.g. via `pzero = 0` in `brms::bf()`, or any
# other intercept-only dpar), into a draws x length(idx) matrix.
#' @keywords internal
.subset_dpar <- function(x, idx, n_draws) {
  if (is.matrix(x)) {
    return(x[, idx, drop = FALSE])
  }
  matrix(x, nrow = n_draws, ncol = length(idx))
}

#' @rdname rbetadiscrete
#' @export
log_lik_betadiscrete <- function(i, prep) {
  y <- prep$data$Y[i]
  k <- prep$data$vint1[i]

  mu <- brms::get_dpar(prep, "mu", i = i)
  phi <- brms::get_dpar(prep, "phi", i = i)
  pzero <- brms::get_dpar(prep, "pzero", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(numeric(0))
  }

  y_vec <- rep(y, length.out = n_draws)

  ll <- dbetadiscrete(
    x = y_vec,
    mu = mu,
    phi = phi,
    k = k,
    pzero = pzero,
    log = TRUE
  )
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
  pzero <- brms::get_dpar(prep, "pzero", i = i)

  n_draws <- length(mu)
  if (n_draws == 0) {
    return(matrix(numeric(0), ncol = 1))
  }

  final_out <- rbetadiscrete(
    n = n_draws,
    mu = mu,
    phi = phi,
    k = k,
    pzero = pzero
  )
  as.matrix(final_out)
}


#' @rdname rbetadiscrete
#' @export
posterior_epred_betadiscrete <- function(prep) {
  mu <- brms::get_dpar(prep, "mu") # draws x obs, or draws-length vector if no predictors
  phi <- brms::get_dpar(prep, "phi") # draws x obs, or draws-length vector if no predictors
  pzero <- brms::get_dpar(prep, "pzero") # draws x obs, or draws-length vector if fixed/intercept-only
  k <- prep$data$vint1 # length = obs

  n_draws <- prep$ndraws
  n_obs <- length(k)
  epred <- matrix(NA_real_, nrow = n_draws, ncol = n_obs)

  # E[R | alpha, beta] = k - sum_{j=1}^{k-1} F_B(j/k; alpha, beta)  (eq. 4)
  # Grouping by unique k avoids redundant computation when (as is typical)
  # all observations share the same number of rating categories. The hurdle
  # simply rescales the non-zero part of the expectation by (1 - pzero).
  for (k_val in unique(k)) {
    idx <- which(k == k_val)
    mu_sub <- .subset_dpar(mu, idx, n_draws)
    phi_sub <- .subset_dpar(phi, idx, n_draws)
    pzero_sub <- .subset_dpar(pzero, idx, n_draws)

    if (k_val <= 1) {
      epred[, idx] <- k_val * (1 - pzero_sub)
      next
    }

    alpha_sub <- mu_sub * phi_sub * 2
    beta_sub <- (1 - mu_sub) * phi_sub * 2

    cum_sum <- matrix(0, nrow = n_draws, ncol = length(idx))
    for (j in seq_len(k_val - 1)) {
      cum_sum <- cum_sum + stats::pbeta(j / k_val, alpha_sub, beta_sub)
    }

    epred[, idx] <- (1 - pzero_sub) * (k_val - cum_sum)
  }

  epred
}
