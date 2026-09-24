#' @title Gaussbit (Gaussian-Probit) Model of Choice and Reaction Time
#'
#' @description
#' Density, random generation, and `brms` custom family for a descriptive joint
#' model of a reaction time and a binary choice: the RT is **Gaussian**, the
#' choice is a **probit**, and the two are linked by a trial-level correlation
#' `rho`. It is the choice-and-RT counterpart of fitting `brms::gaussian()` to
#' the RTs and `brms::bernoulli("probit")` to the choices, and it reduces to
#' exactly that at `rho = 0`. The name is the two halves run together:
#' **Gauss**ian + pro**bit**.
#'
#' It is not a process model. It is the baseline the process models
#' ([cogmod_lnr()], [cogmod_rdm()], [cogmod_ddm()], [cogmod_lba2()]) can be
#' measured against: it takes the same `rt | dec(response)` data, so
#' `brms::loo_compare()` puts it side by side with them, and the difference
#' says how much a sequential sampling model gains over the analysis most
#' people run by default.
#'
#' Functions:
#' - `rcogmod_gaussbit()`: Simulates random draws.
#' - `dcogmod_gaussbit()`: Computes the density (likelihood).
#' - `cogmod_gaussbit()`: Creates a `brms::custom_family()`.
#' - `cogmod_gaussbit_stanvars()`: Generates the `stanvars` to pass to `brm()`.
#'
#' @details
#' # Parameterization
#'
#' On every trial the reaction time is `t = mu + sigma * u`, with `u` a standard
#' Normal, and the choice is `dec = 1` whenever a latent decision variable
#' `w ~ Normal(mudec, 1)` is positive. `u` and `w` are bivariate Normal with
#' correlation `tanh(rho)`. The joint density of a response is then the
#' Gaussian density of `t` times the probit probability of the choice given `t`:
#'
#' \deqn{f(t, 1) = \frac{1}{\sigma} \phi(u) \, \Phi\left(\frac{mudec + r u}{\sqrt{1 - r^2}}\right), \quad r = \tanh(rho)}
#'
#' and `f(t, 0)` the same with the argument of `Phi` negated. Two things follow
#' from the construction, whatever `rho` is:
#'
#' - The **marginal RT distribution is exactly** `Normal(mu, sigma)`. `mu` is
#'   the mean RT and `sigma` its SD, as in `brms::gaussian()`.
#' - The **marginal choice probability is exactly** `pnorm(mudec)`. `mudec` is
#'   the probit of `P(dec = 1)`, as in `brms::bernoulli("probit")`, so its
#'   coefficients are the same quantities a probit regression reports.
#'
#' `rho` changes neither marginal. What it moves is how the choice depends on
#' the time, which is the pattern the process models are built to explain:
#' the probit of `P(dec = 1 | t)` is linear in `t`, with a slope set by `rho`.
#' With `dec` coding errors, `rho > 0` gives **slow errors** - the slower the
#' trial, the more likely an error - and `rho < 0` **fast errors**. Recoding
#' `dec` the other way round flips the sign of both `mudec` and `rho` and leaves
#' the likelihood unchanged.
#'
#' # `rho` is on the Fisher-z scale
#'
#' `brms` has no link onto `(-1, 1)`, so the correlation is estimated as its
#' Fisher-z transform: the correlation itself is `tanh(rho)`, and `rho` is free
#' on the whole real line behind an `identity` link. The two are close for
#' small values (`tanh(0.3) = 0.29`) and part company above that
#' (`tanh(0.5) = 0.46`, `tanh(1) = 0.76`); report `tanh()` of the estimate.
#'
#' There is only one scale to think about: `rho` means the same thing whether
#' it is modelled in `bf()` or left out and estimated as a constant, which is
#' not true of any dpar behind a non-identity link. `rho = 0` is independence
#' either way.
#'
#' The Stan code never forms `tanh(rho)` at all. Since `sqrt(1 - tanh(z)^2) =
#' 1 / cosh(z)`, the argument of `Phi` above is `mudec * cosh(rho) + u *
#' sinh(rho)`: no division, no square root, and no loss of precision as the
#' correlation approaches one.
#'
#' # The same model as `gaussian()` + `bernoulli("probit")`
#'
#' With `rho` fixed at zero this family *is* the multivariate model below, with
#' the same likelihood, observation by observation:
#'
#' ```r
#' # these two are the same model
#' brms::bf(RT | dec(Error) ~ Condition, mudec ~ Condition, rho = 0,
#'          family = cogmod_gaussbit())
#' brms::bf(RT ~ Condition) +
#'   brms::bf(Error ~ Condition, family = brms::bernoulli("probit")) +
#'   brms::set_rescor(FALSE)
#' ```
#'
#' The priors match too. `sigma` has the `log` link `brms::gaussian()` gives it,
#' and [cogmod_priors()] leaves `mu` and `sigma` to the `brms` defaults, gives
#' `mudec` the `student_t(3, 0, 2.5)` intercept `brms::bernoulli()` would, and
#' leaves its slopes flat as `brms` does. The one prior it adds is on `rho`.
#'
#' Freeing `rho` (`rho ~ 1`, or leaving it out of `bf()`) therefore adds exactly
#' one thing to the default analysis - the trial-level coupling of speed and
#' accuracy - and the difference in `elpd` between the two fits measures what
#' that coupling alone is worth, with the RT distribution held Gaussian. A
#' process model's gain over the `rho = 0` fit is what it adds over the
#' default analysis; its gain over the free-`rho` fit is what it adds beyond a
#' linear speed-accuracy dependence.
#'
#' The multivariate form is worth knowing about, but it is not a substitute
#' once `rho` is free. `brms::bf(Error ~ Condition + RT, family =
#' brms::bernoulli("probit"))` reaches the same likelihood by entering the RT
#' as a predictor, but `posterior_predict()` then simulates the choices from
#' the *observed* RTs, so its posterior predictive checks are not joint, and its
#' choice coefficients become effects conditional on the RT rather than the
#' marginal ones above.
#'
#' # What it does not do
#'
#' - **No `ndt`.** The Gaussian's location is free, so a shift of the whole
#'   distribution is `mu` itself: a non-decision time could not be told apart
#'   from it at all, not even weakly.
#' - **No outlier component.** The Gaussian has support on the whole real line
#'   and needs no floor, and the baseline is meant to be the plain default
#'   analysis.
#' - **Negative reaction times have density**, as they do under
#'   `brms::gaussian()`. The mass involved is small for ordinary RT data -
#'   `pnorm(0, 0.6, 0.15)` is 3e-5 - and truncating at zero would cost a
#'   normal tail per observation and break the exact correspondence above.
#' - **No `cens()`.** The choice is already modelled, through `dec()`.
#'
#' # Fitting
#'
#' ```r
#' f <- brms::bf(RT | dec(Error) ~ Condition, mudec ~ Condition, rho ~ 1,
#'               family = cogmod_gaussbit())
#' brms::brm(f, data = df,
#'           prior    = cogmod_priors(f, df),
#'           stanvars = cogmod_stanvars(f))
#' ```
#'
#' `sigma` left out of `bf()` is a constant SD; `sigma ~ Condition` models it,
#' on the log scale.
#'
#' `posterior_predict()` returns simulated reaction times and choices jointly, a
#' draws x 2 matrix as for the other choice families, and `posterior_epred()`
#' the expected reaction time, which is `mu`.
#'
#' @param n Number of simulated trials. If `length(n) > 1`, the length is taken
#'   to be the number required.
#' @param mu Mean reaction time, in seconds. Unbounded, being a location.
#' @param sigma SD of the reaction time, in seconds. Must be positive.
#' @param mudec The probit of the probability of choosing option 1: that
#'   probability is `pnorm(mudec)`. Unbounded.
#' @param rho The trial-level correlation between the reaction time and the
#'   latent decision variable, on the **Fisher-z** scale: the correlation is
#'   `tanh(rho)`. Unbounded; `0` is independence.
#'
#' @return `rcogmod_gaussbit()` returns a data frame with `n` rows and two
#'   columns:
#'   \item{rt}{The simulated reaction time.}
#'   \item{response}{The choice, coded `0` or `1`, matching the `dec()` coding
#'     used by the `brms` family.}
#'
#'   `dcogmod_gaussbit()` returns the joint density of each reaction time and
#'   choice - the log density if `log = TRUE` - recycled to the length of the
#'   longest argument. `cogmod_gaussbit()` returns a `brms::custom_family`
#'   object, to put on a `brms::bf()` formula. `cogmod_gaussbit_stanvars()`
#'   returns a `brms::stanvars` object holding the family's Stan `functions`
#'   block, to pass to `brms::brm()`, and `cogmod_gaussbit_lpdf_expose()`
#'   compiles that Stan code and returns it as an R function, for checking the
#'   density outside of a model. The remaining functions are `brms`
#'   post-processing methods, called by `brms` rather than directly:
#'   `log_lik_cogmod_gaussbit()` returns a numeric vector holding one
#'   log-likelihood value per posterior draw for observation `i`,
#'   `posterior_predict_cogmod_gaussbit()` a draws x 2 matrix of reaction
#'   times and choices simulated for observation `i`, and
#'   `posterior_epred_cogmod_gaussbit()` a draws x observations matrix of
#'   expected reaction times.
#'
#' @seealso [rcogmod_lnr()], [rcogmod_rdm()], [rcogmod_ddm()], [rcogmod_lba2()]
#'
#' @examples
#' # 10% errors (dec = 1), slower trials more error-prone
#' d <- rcogmod_gaussbit(5000, mu = 0.6, sigma = 0.15, mudec = -1.3,
#'                       rho = 0.4)
#' mean(d$response)
#' tapply(d$rt, d$response, mean)
#'
#' # The correlation itself is tanh(rho)
#' tanh(0.4)
#'
#' # At rho = 0 the density is a Gaussian times a probit
#' dcogmod_gaussbit(0.7, mu = 0.6, sigma = 0.15, mudec = -1.3, rho = 0,
#'                  response = 1)
#' dnorm(0.7, 0.6, 0.15) * pnorm(-1.3)
#'
#' @export
rcogmod_gaussbit <- function(n, mu = 0.6, sigma = 0.15, mudec = -1.3,
                             rho = 0) {
  p <- .prepare_gaussbit(n = n, mu = mu, sigma = sigma, mudec = mudec,
                         rho = rho)
  m <- p$ndraws
  u <- stats::rnorm(m)
  # w = mudec + tanh(rho) * u + sqrt(1 - tanh(rho)^2) * e, the last factor
  # written as 1 / cosh(rho), which is the same number without the cancellation
  # as the correlation approaches one.
  w <- p$mudec + tanh(p$rho) * u + stats::rnorm(m) / cosh(p$rho)
  data.frame(rt = p$mu + p$sigma * u, response = as.numeric(w > 0))
}


#' @rdname rcogmod_gaussbit
#' @param x The observed reaction time.
#' @param response The observed choice, 0 or 1.
#' @param log Logical; if TRUE, returns the log-density. Default: FALSE.
#' @export
dcogmod_gaussbit <- function(x, mu = 0.6, sigma = 0.15, mudec = -1.3,
                             rho = 0, response, log = FALSE) {
  p <- .prepare_gaussbit(x = x, response = response, mu = mu, sigma = sigma,
                         mudec = mudec, rho = rho)
  ld <- .gaussbit_ldens(p$x, p$response, p$mu, p$sigma, p$mudec, p$rho)
  if (log) ld else exp(ld)
}


# Internals ---------------------------------------------------------------

# Log-density of (t, k), elementwise. The mirror of cogmod_gaussbit_lpdf()
# below, line for line: the Gaussian log-density of t, plus the log probit of
# the choice given t, whose argument mudec * cosh(rho) + u * sinh(rho) is
# (mudec + r * u) / sqrt(1 - r^2) at r = tanh(rho).
#
# A time or a parameter that is not a finite number gives -Inf rather than NA,
# the convention of the mixture densities (see .dens_mask() in core_shifted.R):
# an infinite time has zero density, and NA * 0 would otherwise leak a NaN
# out of the u * sinh(rho) term at rho = 0.
#' @keywords internal
.gaussbit_ldens <- function(t, k, mu, sigma, mudec, rho) {
  u <- (t - mu) / sigma
  a <- mudec * cosh(rho) + u * sinh(rho)
  # (2k - 1) rather than ifelse(k == 1, a, -a), whose result takes the length
  # of `k` rather than of `a`.
  ld <- -0.5 * u^2 - log(sigma) - 0.5 * log(2 * pi) +
    stats::pnorm((2 * k - 1) * a, log.p = TRUE)
  ld[!is.finite(ld)] <- -Inf
  ld
}


# Validate and recycle. `mu`, `mudec` and `rho` are locations on the whole real
# line and are not checked; only `sigma` has a bound, and the response must be
# one of the two options.
#' @keywords internal
.prepare_gaussbit <- function(x = NULL, n = NULL, response = NULL, mu, sigma,
                              mudec, rho) {
  if (any(sigma <= 0, na.rm = TRUE)) stop("`sigma` must be positive.", call. = FALSE)
  if (!is.null(response) && any(!response %in% c(0, 1), na.rm = TRUE)) {
    stop("`response` must be 0 or 1.", call. = FALSE)
  }
  lens <- c(length(mu), length(sigma), length(mudec), length(rho))
  if (!is.null(x)) {
    if (is.null(response)) {
      stop("`response` must be provided alongside `x`.", call. = FALSE)
    }
    # A zero-length quantile gives a zero-length answer, as in base R.
    m <- if (length(x) == 0L) 0L else max(length(x), length(response), lens)
  } else {
    if (length(n) > 1) n <- length(n)
    if (length(n) != 1 || n < 0 || n != floor(n)) {
      stop("n must be a single non-negative integer.", call. = FALSE)
    }
    m <- if (n == 0) 0L else max(n, lens)
  }
  p <- list(mu = rep_len(mu, m), sigma = rep_len(sigma, m),
            mudec = rep_len(mudec, m), rho = rep_len(rho, m), ndraws = m)
  if (!is.null(x)) {
    p$x <- rep_len(x, m)
    p$response <- rep_len(response, m)
  }
  p
}


# Family ------------------------------------------------------------------

#' @rdname rcogmod_gaussbit
#' @param link_mu,link_sigma,link_mudec,link_rho Link functions for the
#'   parameters. `sigma` defaults to `"log"`, the link `brms::gaussian()` uses,
#'   so that the family matches the default analysis exactly. `rho` is on the
#'   Fisher-z scale behind its `"identity"` link; see Details.
#' @export
cogmod_gaussbit <- function(link_mu = "identity", link_sigma = "log",
                            link_mudec = "identity", link_rho = "identity") {
  brms::custom_family(
    name = "cogmod_gaussbit",
    dpars = c("mu", "sigma", "mudec", "rho"),
    links = c(link_mu, link_sigma, link_mudec, link_rho),
    lb = c(NA, 0, NA, NA),
    ub = c(NA, NA, NA, NA),
    type = "real",
    vars = "dec[n]"
  )
}


# Stanvars ----------------------------------------------------------------

#' @keywords internal
.cogmod_gaussbit_lpdf <- function() {
  paste0(.LOG_PHI_STAN_PRELUDE, "
// Log-likelihood for one observation from the Gaussian-probit model of choice
// and RT: a Gaussian reaction time and a probit choice, correlated tanh(rho)
// at the trial level. dcogmod_gaussbit() in R, line for line.
// Y: observed reaction time.
// mu, sigma: mean and SD of the reaction time (sigma > 0).
// mudec: probit of P(dec = 1), marginally over the reaction time.
// rho: Fisher-z of the correlation between the reaction time and the latent
//      decision variable.
// dec: the observed choice, 0 or 1.
//
// The choice term is Phi((mudec + r u) / sqrt(1 - r^2)) at r = tanh(rho),
// written as Phi(mudec * cosh(rho) + u * sinh(rho)) - the same number, since
// 1 / sqrt(1 - tanh(z)^2) = cosh(z), with no division, no square root and no
// cancellation as the correlation approaches one. The tail goes through
// cogmod_log_Phi(), not std_normal_lcdf(), whose partials are approximate:
// a slow response at a large rho puts the argument far out in either tail.
//
// The Gaussian term is written out rather than left to normal_lpdf() so that
// the standardized time is formed once and shared with the choice term.
real cogmod_gaussbit_lpdf(real Y, real mu, real sigma, real mudec, real rho,
                          int dec) {
    if (sigma <= 0) return negative_infinity();
    if (dec < 0 || dec > 1) return negative_infinity();
    real u = (Y - mu) / sigma;
    real a = mudec * cosh(rho) + u * sinh(rho);
    return -0.5 * square(u) - log(sigma) - 0.91893853320467274
           + cogmod_log_Phi(dec == 1 ? a : -a);
}
")
}


#' @rdname rcogmod_gaussbit
#' @examples
#' \dontrun{
#' # Needs cmdstanr and a CmdStan toolchain, which live outside CRAN - see the
#' # package website to install them. Not run under R CMD check, which executes
#' # every example in one R session: once brms has fitted a model there (the
#' # cogmod_inits() and p_outlier() examples do), rstan is live in the process
#' # and loading an exposed Stan function next to it segfaults on Linux.
#' lpdf <- cogmod_gaussbit_lpdf_expose()
#' lpdf(Y = 0.7, mu = 0.6, sigma = 0.15, mudec = -1.3, rho = 0.4, dec = 1)
#' }
#'
#' @export
cogmod_gaussbit_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")
  stancode <- paste0("functions {\n", .cogmod_gaussbit_lpdf(), "\n}")
  mod <- cmdstanr::cmdstan_model(cmdstanr::write_stan_file(stancode))
  mod$expose_functions()
  mod$functions$cogmod_gaussbit_lpdf
}


#' @rdname rcogmod_gaussbit
#' @export
cogmod_gaussbit_stanvars <- function() {
  brms::stanvar(scode = .cogmod_gaussbit_lpdf(), block = "functions")
}


# brms methods ------------------------------------------------------------

#' @rdname rcogmod_gaussbit
#' @inheritParams rcogmod_betagate
#' @export
log_lik_cogmod_gaussbit <- function(i, prep) {
  y <- prep$data$Y[i]
  if (is.na(y)) return(NA_real_)
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  mudec <- brms::get_dpar(prep, "mudec", i = i)
  rho <- brms::get_dpar(prep, "rho", i = i)
  n_draws <- max(length(mu), length(sigma), length(mudec), length(rho))
  if (n_draws == 0) return(numeric(0))
  k <- .dec_from_prep(prep, i)
  .gaussbit_ldens(rep_len(y, n_draws), rep_len(k, n_draws),
                  rep_len(mu, n_draws), rep_len(sigma, n_draws),
                  rep_len(mudec, n_draws), rep_len(rho, n_draws))
}


#' @rdname rcogmod_gaussbit
#' @inheritParams rcogmod_betagate
#' @export
posterior_predict_cogmod_gaussbit <- function(i, prep, ...) {
  mu <- brms::get_dpar(prep, "mu", i = i)
  sigma <- brms::get_dpar(prep, "sigma", i = i)
  mudec <- brms::get_dpar(prep, "mudec", i = i)
  rho <- brms::get_dpar(prep, "rho", i = i)
  n_draws <- max(length(mu), length(sigma), length(mudec), length(rho))
  # A plain matrix rather than rcogmod_gaussbit()'s data frame: brms calls
  # this once per observation, and building a data frame only to convert it
  # back is most of the cost of the call (see .rchoice() in core_choice.R).
  u <- stats::rnorm(n_draws)
  w <- mudec + tanh(rho) * u + stats::rnorm(n_draws) / cosh(rho)
  cbind(rt = mu + sigma * u, response = as.numeric(w > 0))
}


#' @rdname rcogmod_gaussbit
#' @export
posterior_epred_cogmod_gaussbit <- function(prep) {
  # The marginal RT distribution is exactly Normal(mu, sigma), whatever rho is.
  brms::get_dpar(prep, "mu")
}
