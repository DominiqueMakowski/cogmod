# Simulation ---------------------------------------------------------

#' @title Choice-Confidence (CHOCO) Model
#'
#' @description
#' The Choice-Confidence (CHOCO) model is useful to model data from subjective scales, such as Likert-type or analog scales, in which the left and the right side correspond to different processes or higher order categorical responses (e.g., "disagree" vs. "agree", "true" vs. "false"). They can be used to jointly model choice (left or right) and confidence (the degree of left or right).
#'
#' To better represent potentially bimodal bounded data with zeros and ones, CHOCO models consist of two separate ordered beta distributions (see Kubinec, 2023) modeling the left and the right hand side of the scale.
#'
#' Regarding Ordered Beta distributions, see:
#' - Stan code: https://github.com/saudiwin/ordbetareg/blob/master/beta_logit.stan
#' - Python code: https://github.com/saudiwin/ordbetareg_py/blob/main/ordbetareg/model.py
#'
#' The CHOCO model is available following two different parametrizations:
#' - **CHOCO7**: A 7-parameter model where the right-side parameters (`muright`, `phiright`, `kright`) are independent of the left-side parameters.
#' - **CHOCO "7d"** (default CHOCO): A deviation-based 7-parameter model where the right-side parameters are derived as deviations from the left-side parameters (`mud`, `phid`, `kd`).
#'
#' Functions:
#' - `rchoco()`: simulates random draws from the CHOCO distribution, which models bimodal analog/Likert scale data.
#' - `choco_stanvars()`: Generates a `stanvars` object to pass custom functions to `brms`.
#' - `choco()`: Creates a custom family to be used with `brms`.
#'
#' @param n Number of random draws. Must be a positive integer.
#' @param mu Probability of choosing the right side (relative to the left side).
#'   Determines the overall bias toward the right side. Must be in the range 0-1.
#' @param muleft The center of the left-side beta distribution.
#'   Determines the central tendency of the left-side responses. Must be in the range 0-1.
#' @param phileft The shape parameter of the left-side beta distribution.
#'   Controls the precision (inverse variance) of the left-side responses. Must be positive. Plausible range: (0, Inf).
#' @param kleft The cutoff parameter for the left-side beta distribution.
#'   Determines the probability of extreme values (close to 0) on the left side. Must be in the range 0-1.
#' @param muright The center of the right-side beta distribution.
#'   Determines the central tendency of the right-side responses. Overridden by `mud` if NULL. Must be in the range 0-1.
#' @param phiright The shape parameter of the right-side beta distribution.
#'   Controls the precision (inverse variance) of the right-side responses. Overridden by `phid` if NULL. Must be positive. Plausible range: (0, Inf).
#' @param kright The cutoff parameter for the right-side beta distribution.
#'   Determines the probability of extreme values (close to 1) on the right side. Overridden by `kd` if NULL. Must be in the range 0-1.
#' @param mud Deviation for `muright` on the logit scale relative to `muleft`.
#'   Positive values shift the right-side center further from the left-side center. Plausible range: (-Inf, Inf).
#' @param phid Deviation for `phiright` as a log-multiplier relative to `phileft`.
#'   Positive values increase the precision of the right-side responses relative to the left side. Plausible range: (-Inf, Inf).
#' @param kd Deviation for `kright` on the logit scale relative to `kleft`.
#'   Positive values increase the probability of extreme values on the right side relative to the left side. Plausible range: (-Inf, Inf).
#'
#' @examples
#' hist(rchoco(3000, mu=0.5, muleft=0.5, phileft=3, kleft=0.95), breaks=100)
#' hist(rchoco(3000, mu=0.6, muleft=0.5, phileft=5, kleft=0.99), breaks=100)
#' hist(rchoco(3000, mu=0.4, muleft=0.5, phileft=3, kleft=0.85), breaks=100)
#'
#' @export
rchoco <- function(n, mu = 0.5, muleft = 0.5, phileft = 3, kleft = 0.95,
                   muright = NULL, phiright = NULL, kright = NULL,
                   mud = 0, phid = 0, kd = 0) {
  # Overall probabilities for choosing left or right side
  p_left <- 1 - mu    # probability for the left side
  p_right <- mu       # probability for the right side

  # Use provided right-side parameters if not NULL, else compute using deviation parameters:
  if (is.null(muright)) muright <- stats::plogis(stats::qlogis(muleft) + mud)
  if (is.null(phiright)) phiright <- phileft * exp(phid)
  if (is.null(kright)) kright <- stats::plogis(stats::qlogis(kleft) + kd)

  # Compute the discrete endpoint probabilities for 0 (left side) and 1 (right side)
  p0 <- 1 - stats::plogis(stats::qlogis(muleft) - stats::qlogis(1 - kleft))      # left-side endpoint probability
  p1 <- stats::plogis(stats::qlogis(muright) - stats::qlogis(kright))           # right-side endpoint probability

  # Pre-allocate output vector for efficiency
  y <- numeric(n)

  # Vectorized random assignment of sides for all n draws
  sides <- sample(c("left", "right"), size = n, replace = TRUE, prob = c(p_left, p_right))

  ## Process draws assigned to the left side
  left_idx <- which(sides == "left")
  if (length(left_idx) > 0) {
    left_u <- stats::runif(length(left_idx))
    y[left_idx[left_u < p0]] <- 0
    left_cont_idx <- left_idx[left_u >= p0]
    if (length(left_cont_idx) > 0) {
      x0 <- stats::rbeta(length(left_cont_idx),
                         shape1 = muleft * phileft,
                         shape2 = (1 - muleft) * phileft)
      y[left_cont_idx] <- 0.5 - x0 / 2
    }
  }

  ## Process draws assigned to the right side
  right_idx <- which(sides == "right")
  if (length(right_idx) > 0) {
    right_u <- stats::runif(length(right_idx))
    y[right_idx[right_u < p1]] <- 1
    right_cont_idx <- right_idx[right_u >= p1]
    if (length(right_cont_idx) > 0) {
      x1 <- stats::rbeta(length(right_cont_idx),
                         shape1 = muright * phiright,
                         shape2 = (1 - muright) * phiright)
      y[right_cont_idx] <- 0.5 + x1 / 2
    }
  }

  y
}

# Stan --------------------------------------------------------------------
#' @rdname rchoco
#' @param type The type of CHOCO model. "choco" is the deviation based model and "choco7" the independent parametrization.
#' @examples
#' choco_stanvars("choco7")
#' @export
choco_stanvars <- function(type = "choco") {
  if (type == "choco7") { # 7-parameter model with independent right-side parameters
    stancode <- brms::stanvar(scode = "
real choco7_lpdf(real y, real mu, real muleft, real phileft,
                real muright, real phiright, real kleft, real kright) {
  real eps = 1e-8;
  real p_left = 1 - mu;
  real p_right = mu;
  if (y < eps) {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    return log(p_left) + log(fmax(p0, eps));
  } else if (y > (1 - eps)) {
    real p1 = inv_logit(logit(muright) - logit(kright));
    return log(p_right) + log(fmax(p1, eps));
  } else {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    real p1 = inv_logit(logit(muright) - logit(kright));
    real log1m_p0 = log1m(fmin(p0, 1 - eps));
    real log1m_p1 = log1m(fmin(p1, 1 - eps));
    if (abs(y - 0.5) < eps) {
      real x0 = 2 * (0.5 - (0.5 - eps));
      real x1 = 2 * ((0.5 + eps) - 0.5);
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    } else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    } else {
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")
  } else if (type == "choco") { # 7-parameter deviation-based model
    stancode <- brms::stanvar(scode = "
real choco_lpdf(real y, real mu, real muleft, real mud, real phileft, real phid, real kleft, real kd) {
  real eps = 1e-8;
  real p_left = 1 - mu;
  real p_right = mu;
  // Derive right-side parameters as deviations from the left-side:
  real muright = inv_logit(logit(muleft) + mud);
  real phiright = phileft * exp(phid);
  real kright = inv_logit(logit(kleft) + kd);
  if (y < eps) {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    return log(p_left) + log(fmax(p0, eps));
  } else if (y > (1 - eps)) {
    real p1 = inv_logit(logit(muright) - logit(kright));
    return log(p_right) + log(fmax(p1, eps));
  } else {
    real p0 = 1 - inv_logit(logit(muleft) - logit(1 - kleft));
    real p1 = inv_logit(logit(muright) - logit(kright));
    real log1m_p0 = log1m(fmin(p0, 1 - eps));
    real log1m_p1 = log1m(fmin(p1, 1 - eps));
    if (abs(y - 0.5) < eps) {
      real x0 = 2 * (0.5 - (0.5 - eps));
      real x1 = 2 * ((0.5 + eps) - 0.5);
      real dens_left = p_left * exp(log1m_p0 + log(2) +
                                    beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft));
      real dens_right = p_right * exp(log1m_p1 + log(2) +
                                      beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright));
      real avg_dens = (dens_left + dens_right) / 2;
      return log(fmax(avg_dens, eps));
    } else if (y < 0.5) {
      real x0 = 2 * (0.5 - y);
      return log(p_left) + log1m_p0 + log(2) +
             beta_proportion_lpdf(x0 | fmin(muleft, 1 - eps), phileft);
    } else {
      real x1 = fmax(2 * (y - 0.5), eps);
      return log(p_right) + log1m_p1 + log(2) +
             beta_proportion_lpdf(x1 | fmin(muright, 1 - eps), phiright);
    }
  }
}
", block = "functions")
  } else {
    stop("Invalid type. Choose either 'choco7' or 'choco'.")
  }
  stancode
}



#' @rdname rchoco
#' @param link_mu,link_muleft,link_mud,link_phileft,link_phid,link_kleft,link_kd,link_muright,link_phiright,link_kright Character strings specifying the link functions for the parameters.
#' @examples
#' choco()
#' @export
choco <- function(link_mu = "logit", link_muleft = "logit", link_mud = "identity",
                    link_phileft = "softplus", link_phid = "identity",
                    link_kleft = "logit", link_kd = "identity") {
  brms::custom_family(
    "choco",
    dpars = c("mu", "muleft", "mud", "phileft", "phid", "kleft", "kd"),
    links = c(link_mu, link_muleft, link_mud,
              link_phileft, link_phid,
              link_kleft, link_kd)
  )
}


#' @rdname rchoco
#' @examples
#' choco7()
#' @export
choco7 <- function(link_mu = "logit", link_muleft = "logit", link_muright = "logit",
                   link_phileft = "softplus", link_phiright = "softplus",
                   link_kleft = "logit", link_kright = "logit") {
  brms::custom_family(
    "choco7",
    dpars = c("mu", "muleft", "muright", "phileft", "phiright", "kleft", "kright"),
    links = c(link_mu, link_muleft, link_muright,
              link_phileft, link_phiright,
              link_kleft, link_kright)
  )
}


# Predict -----------------------------------------------------------------

#' @rdname rchoco
#' @param i,prep For brms' functions to run: index of the observation and a `brms` preparation object.
#' @param ... Additional arguments.
#' @export
posterior_predict_choco7 <- function(i, prep, ...) {
  mu       <- brms::get_dpar(prep, "mu", i = i)
  muleft   <- brms::get_dpar(prep, "muleft", i = i)
  muright  <- brms::get_dpar(prep, "muright", i = i)
  phileft  <- brms::get_dpar(prep, "phileft", i = i)
  phiright <- brms::get_dpar(prep, "phiright", i = i)
  kleft    <- brms::get_dpar(prep, "kleft", i = i)
  kright   <- brms::get_dpar(prep, "kright", i = i)

  yrep <- mapply(function(mu_i, muleft_i, muright_i, phileft_i, phiright_i, kleft_i, kright_i) {
    rchoco(1,
           mu       = mu_i,
           muleft   = muleft_i,
           muright  = muright_i,
           phileft  = phileft_i,
           phiright = phiright_i,
           kleft    = kleft_i,
           kright   = kright_i)
  }, mu, muleft, muright, phileft, phiright, kleft, kright)

  yrep
}

#' @rdname rchoco
#' @export
posterior_predict_choco <- function(i, prep, ...) {
  mu      <- brms::get_dpar(prep, "mu", i = i)
  muleft  <- brms::get_dpar(prep, "muleft", i = i)
  mud     <- brms::get_dpar(prep, "mud", i = i)
  phileft <- brms::get_dpar(prep, "phileft", i = i)
  phid    <- brms::get_dpar(prep, "phid", i = i)
  kleft   <- brms::get_dpar(prep, "kleft", i = i)
  kd      <- brms::get_dpar(prep, "kd", i = i)

  yrep <- mapply(function(mu_i, muleft_i, mud_i, phileft_i, phid_i, kleft_i, kd_i) {
    rchoco(1,
           mu      = mu_i,
           muleft  = muleft_i,
           mud     = mud_i,
           phileft = phileft_i,
           phid    = phid_i,
           kleft   = kleft_i,
           kd      = kd_i)
  }, mu, muleft, mud, phileft, phid, kleft, kd)

  yrep
}


#' @rdname rchoco
#' @param y The observed response value. Must be in the range 0-1.
#' @param log Logical. If TRUE, returns the log-density. Default: TRUE.
#' @export
dchoco <- function(y, mu, muleft, phileft, kleft, muright = NULL, phiright = NULL, kright = NULL,
                   mud = 0, phid = 0, kd = 0, log = TRUE) {
  eps <- 1e-6

  # Compute right-side parameters if not provided
  if (is.null(muright)) muright <- stats::plogis(stats::qlogis(muleft) + mud)
  if (is.null(phiright)) phiright <- phileft * exp(phid)
  if (is.null(kright)) kright <- stats::plogis(stats::qlogis(kleft) + kd)

  # Compute probabilities for left and right sides
  p_left <- 1 - mu
  p_right <- mu

  # Compute densities
  if (y < eps) {
    log_density <- log(p_left) + log(1 - stats::plogis(stats::qlogis(muleft) - stats::qlogis(1 - kleft)))
  } else if (y > (1 - eps)) {
    log_density <- log(p_right) + log(stats::plogis(stats::qlogis(muright) - stats::qlogis(kright)))
  } else if (y < 0.5) {
    x0 <- 2 * (0.5 - y)
    log_density <- log(p_left) + log(1 - stats::plogis(stats::qlogis(muleft) - stats::qlogis(1 - kleft))) +
      log(2) + stats::dbeta(x0, shape1 = muleft * phileft, shape2 = (1 - muleft) * phileft, log = TRUE)
  } else {
    x1 <- 2 * (y - 0.5)
    log_density <- log(p_right) + log(stats::plogis(stats::qlogis(muright) - stats::qlogis(kright))) +
      log(2) + stats::dbeta(x1, shape1 = muright * phiright, shape2 = (1 - muright) * phiright, log = TRUE)
  }

  if (log) log_density else exp(log_density)
}

#' @rdname rchoco
#' @export
log_lik_choco <- function(i, prep) {
  # Extract the i-th observed response
  y <- prep$data$Y[i]

  # Extract posterior draws for each parameter for observation i
  mu      <- brms::get_dpar(prep, "mu", i = i)
  muleft  <- brms::get_dpar(prep, "muleft", i = i)
  phileft <- brms::get_dpar(prep, "phileft", i = i)
  kleft   <- brms::get_dpar(prep, "kleft", i = i)
  mud     <- brms::get_dpar(prep, "mud", i = i)
  phid    <- brms::get_dpar(prep, "phid", i = i)
  kd      <- brms::get_dpar(prep, "kd", i = i)

  # Compute log-likelihood using dchoco
  dchoco(y = y, mu = mu, muleft = muleft, phileft = phileft, kleft = kleft,
         mud = mud, phid = phid, kd = kd, log = TRUE)
}
