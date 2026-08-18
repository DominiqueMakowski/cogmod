# Reading the family off what the user passed
# ===========================================
#
# cogmod_priors(), cogmod_inits() and cogmod_stanvars() all take the model
# rather than the family, so the family is named once - in bf() - and cannot
# drift out of step with the rest of the call. These pull it back out.

# A brmsformula, a brmsfit and a family object are all lists carrying the
# family under `$family`; a family object is its own answer. A plain formula
# carries nothing, and gives NULL so the caller can say so.
#' @keywords internal
.cogmod_family <- function(x) {
  if (inherits(x, c("customfamily", "brmsfamily", "family"))) return(x)
  if (is.list(x)) return(x[["family"]])
  NULL
}

# brms families carry $name; stats families only $family.
#' @keywords internal
.family_name <- function(family) {
  if (is.null(family)) return(NULL)
  if (is.character(family)) return(family)
  nm <- family[["name"]]
  if (is.null(nm)) nm <- family[["family"]]
  nm
}

# The link of every dpar, named. brms stores the first dpar's link as `$link`
# and the rest as `$link_<dpar>`, so the two forms are folded back together
# here. Reading them off the family rather than the registry is what makes
# cogmod_gamma(link_mu = "log") behave.
#' @keywords internal
.family_links <- function(family) {
  dpars <- family[["dpars"]]
  if (is.null(dpars)) return(stats::setNames(character(0), character(0)))
  out <- vapply(seq_along(dpars), function(i) {
    v <- family[[paste0("link_", dpars[i])]]
    if (is.null(v) && i == 1L) v <- family[["link"]]
    if (is.null(v)) NA_character_ else as.character(v)[1]
  }, character(1))
  stats::setNames(out, dpars)
}


# Stable computation of log(1 - exp(x))
# A stable version of log1m_exp (assumes x is scalar or vector)
#' @keywords internal
.log1m_exp <- function(x) {
  # Computes log(1 - exp(x)) in a numerically stable way.
  # x should be negative; for vectorized input, use ifelse.
  ifelse(x < log(0.5), log1p(-exp(x)), log(-expm1(x)))
}






# Stable log of a two-component mixture, mirroring Stan's log_mix():
#   log(theta * exp(lp1) + (1 - theta) * exp(lp2))
# Handles lp2 = -Inf (a component with zero density), which is the normal case
# for shifted distributions evaluated below their shift.
#' @keywords internal
.log_mix <- function(theta, lp1, lp2) {
  a <- log(theta) + lp1
  b <- log1p(-theta) + lp2
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  # Both components vanish; pmax() is -Inf and the shift would give NaN.
  out[!is.finite(m)] <- -Inf
  out
}


# Rejection sampling algorithm by Robert (Stat. Comp (1995), 5, 121-5)
# for simulating from the truncated normal distribution.
# Copied from the msm package
#' @keywords internal
.rnorm_truncated <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  if (length(n) > 1) {
    n <- length(n)
  }
  mean <- rep(mean, length = n)
  sd <- rep(sd, length = n)
  lower <- rep(lower, length = n)
  upper <- rep(upper, length = n)
  ret <- numeric(n)
  ind <- seq(length.out = n)

  sdzero <- sd < .Machine$double.eps
  ## return the mean, unless mean is outside the range, then return nan
  sdna <- sdzero & ((mean < lower) | (mean > upper))

  lower <- (lower - mean) / sd ## Algorithm works on mean 0, sd 1 scale
  upper <- (upper - mean) / sd
  nas <- is.na(mean) | is.na(sd) | is.na(lower) | is.na(upper) | sdna
  if (any(nas)) warning("NAs produced")
  ## Different algorithms depending on where upper/lower limits lie.
  alg <- ifelse(
    ((lower > upper) | nas),
    -1, # return NaN
    ifelse(
      sdzero,
      4, # SD zero, so set the sampled value to the mean.
      ifelse(
        ((lower < 0 & upper == Inf) |
          (lower == -Inf & upper > 0) |
          (is.finite(lower) & is.finite(upper) & (lower < 0) & (upper > 0) & (upper - lower > sqrt(2 * pi)))
        ),
        0, # standard "simulate from normal and reject if outside limits" method. Use if bounds are wide.
        ifelse(
          (lower >= 0 & (upper > lower + 2 * sqrt(exp(1)) /
            (lower + sqrt(lower^2 + 4)) * exp((lower * 2 - lower * sqrt(lower^2 + 4)) / 4))),
          1, # rejection sampling with exponential proposal. Use if lower >> mean
          ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1)) /
            (-upper + sqrt(upper^2 + 4)) * exp((upper * 2 - -upper * sqrt(upper^2 + 4)) / 4)),
          2, # rejection sampling with exponential proposal. Use if upper << mean.
          3
          )
        )
      )
    )
  ) # rejection sampling with uniform proposal. Use if bounds are narrow and central.

  ind.nan <- ind[alg == -1]
  ind.no <- ind[alg == 0]
  ind.expl <- ind[alg == 1]
  ind.expu <- ind[alg == 2]
  ind.u <- ind[alg == 3]
  ind.sd0 <- ind[alg == 4]
  ret[ind.nan] <- NaN
  ret[ind.sd0] <- 0 # SD zero, so set the sampled value to the mean.
  while (length(ind.no) > 0) {
    y <- stats::rnorm(length(ind.no))
    done <- which(y >= lower[ind.no] & y <= upper[ind.no])
    ret[ind.no[done]] <- y[done]
    ind.no <- setdiff(ind.no, ind.no[done])
  }
  stopifnot(length(ind.no) == 0)
  while (length(ind.expl) > 0) {
    a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expl), a) + lower[ind.expl]
    u <- stats::runif(length(ind.expl))
    done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= upper[ind.expl]))
    ret[ind.expl[done]] <- z[done]
    ind.expl <- setdiff(ind.expl, ind.expl[done])
  }
  stopifnot(length(ind.expl) == 0)
  while (length(ind.expu) > 0) {
    a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expu), a) - upper[ind.expu]
    u <- stats::runif(length(ind.expu))
    done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= -lower[ind.expu]))
    ret[ind.expu[done]] <- -z[done]
    ind.expu <- setdiff(ind.expu, ind.expu[done])
  }
  stopifnot(length(ind.expu) == 0)
  while (length(ind.u) > 0) {
    z <- stats::runif(length(ind.u), lower[ind.u], upper[ind.u])
    rho <- ifelse(lower[ind.u] > 0,
      exp((lower[ind.u]^2 - z^2) / 2), ifelse(upper[ind.u] < 0,
        exp((upper[ind.u]^2 - z^2) / 2),
        exp(-z^2 / 2)
      )
    )
    u <- stats::runif(length(ind.u))
    done <- which(u <= rho)
    ret[ind.u[done]] <- z[done]
    ind.u <- setdiff(ind.u, ind.u[done])
  }
  stopifnot(length(ind.u) == 0)
  ret * sd + mean
}
