#' @title Simulate from the Two-Accumulator Racing Diffusion Model (RDM)
#'
#' @description
#' Simulates choice reaction times from a two-accumulator Racing Diffusion Model (RDM).
#' This is a specialized version where exactly two accumulators race towards a
#' common threshold. The model assumes variability in the starting point of the
#' diffusion process, drawn from a uniform distribution. This version is optimized
#' for performance using vectorized operations, and allows drift rates to be zero
#' (a zero-drift accumulator is slow but still finishes; see Details).
#'
#' @details
#' The RDM implemented here follows the formulation where the two accumulators
#' have drift rates `vzero` and `vone`. The diffusion process starts at a point `z` drawn
#' from `Uniform(0, bias)`. The process terminates when either accumulator reaches
#' a threshold `b`. The parameter `bs` is defined as `bs = b - bias`, representing
#' the distance from the maximum starting point `bias` to the threshold `b`.
#' Therefore, the effective distance to threshold for a given trial is
#' `bs = b - z = bs + bias - z`.
#'
#' The finishing time for a single accumulator, given its drift rate `v`, `bs`, `bias`, and `ndt`,
#' is simulated by:
#' 1. Sampling a starting point `z ~ Uniform(0, bias)`.
#' 2. Calculating the distance `bs = bs + bias - z`.
#' 3. If `v > 0`, simulating the time to reach `bs` from an Inverse Gaussian distribution
#'    with mean `bs / v` and shape `bs^2`. This simulation uses an
#'    internal implementation based on Michael et al. (1976).
#' 4. If `v = 0`, drawing the driftless first passage time, which is Levy
#'    distributed: `bs^2 / Z^2` with `Z ~ Normal(0, 1)`. A zero-drift
#'    accumulator is slow, but it still finishes with probability one, so it can
#'    win the race.
#' 5. Adding the non-decision time `ndt` to the finishing times.
#'
#' The function simulates this process for both accumulators using vectorized operations.
#' The accumulator that finishes first determines the response (0 for `vzero`, 1 for `vone`)
#' and the reaction time (RT) for that trial.
#'
#' This implementation is based on the description and parameters used in:
#' Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#' without random between-trial variability: The racing diffusion model of
#' speeded decision making. *Psychonomic Bulletin & Review*, *27*, 911-936.
#' \doi{10.3758/s13423-020-01719-6} (specifically matching the `WaldA` component
#' used within their RDM simulation).
#'
#' @param n Number of trials to simulate. Must be a positive integer.
#' @param vzero Drift rate for the first accumulator. Must be a single non-negative number.
#' @param vone Drift rate for the second accumulator. Must be a single non-negative number.
#' @param bs Threshold parameter, defined as `bs = b - bias`, where `b` is the
#'   decision threshold and `bias` is the maximum starting point. Must be a single
#'   positive number.
#' @param bias Maximum starting point parameter. The starting point for each
#'   accumulator on each trial is drawn from `Uniform(0, bias)`. Must be a single
#'   positive number.
#' @param ndt Non-decision time (encoding and motor time offset). Must be a
#'   single non-negative number.
#'
#' @return A data frame with `n` rows and two columns:
#'   \item{rt}{The simulated reaction time (minimum finishing time across the two accumulators).}
#'   \item{response}{The winning accumulator, coded `0` for `vzero` and `1` for
#'     `vone`, matching the `dec()` coding used by the `brms` families.}
#'
#' @references
#' - Michael, J. R., Schucany, W. R., & Haas, R. W. (1976). Generating Random Variates Using
#'     Transformations with Multiple Roots. *The American Statistician*, *30*(2), 88–90. \doi{10.2307/2683801}
#' - Tillman, G., Van Zandt, T., & Logan, G. D. (2020). Sequential sampling models
#'     without random between-trial variability: The racing diffusion model of
#'     speeded decision making. *Psychonomic Bulletin & Review*, *27*, 911-936.
#'     \doi{10.3758/s13423-020-01719-6}
#' - Folks, J. L., & Chhikara, R. S. (1978). The inverse Gaussian distribution and its
#'     statistical application—a review. *Journal of the Royal Statistical Society Series B:
#'     Statistical Methodology*, *40*(3), 263-275.
#'
#' @seealso `rrt_invgaussian`
#' @examples
#' rdm_pos <- rrdm(n = 1000, vzero = 0.8, vone = 0.6, bs = 0.5, bias = 0.2, ndt = 0.15)
#'
#' @export
rrdm <- function(n, vzero, vone, bs, bias, ndt) {
  # Prepare and validate parameters using the helper function
  # This ensures all parameters are recycled to length 'n' (stored in params$n)
  params <- .prepare_rrdm(n = n, vzero = vzero, vone = vone, bs = bs, bias = bias, ndt = ndt)
  nobs <- params$ndraws

  num_accumulators_total <- nobs * 2

  # Simulate starting points (z) for all trials and accumulators
  # params$bias is already recycled to length nobs
  # We need 2*nobs starting points, one for each accumulator in each trial
  z_vec <- stats::runif(num_accumulators_total, min = 0, max = rep(params$bias, each = 2))

  # Calculate distance to threshold (alpha = b - z = bs + bias - z)
  # Interleave bs and bias vectors to match the structure (trial1_acc1, trial1_acc2, ...)
  k_interleaved <- rep(params$bs, each = 2)
  A_interleaved <- rep(params$bias, each = 2)
  alpha_vec <- k_interleaved + A_interleaved - z_vec # Length = nobs * 2

  # Prepare drift rate vector (interleaving vzero_i, vone_i for each trial i)
  v_matrix <- cbind(params$vzero, params$vone) # nobs x 2 matrix
  v_vec <- as.vector(t(v_matrix)) # Length = nobs * 2, ordered (v0_1, v1_1, v0_2, v1_2, ...)

  # Prepare non-decision time vector (interleaving ndt_i, ndt_i for each trial i)
  ndt_vec_interleaved <- rep(params$ndt, each = 2) # Length = nobs * 2

  # Simulate finishing times. Positive drifts give an inverse Gaussian; a zero
  # drift is a driftless diffusion, whose first passage time to `b` is Levy
  # distributed and drawn as b^2 / Z^2 with Z ~ N(0, 1). (Zero drift means the
  # accumulator is slow, not that it never finishes: driftless Brownian motion
  # still reaches any positive level with probability one.)
  finish_time_vec <- numeric(num_accumulators_total)
  pos <- v_vec > 0
  if (any(pos)) {
    finish_time_vec[pos] <- rrt_invgaussian(
      n = sum(pos),
      drift = v_vec[pos],
      bs = alpha_vec[pos],
      ndt = ndt_vec_interleaved[pos]
    )
  }
  if (any(!pos)) {
    nz <- sum(!pos)
    finish_time_vec[!pos] <- (alpha_vec[!pos] / stats::rnorm(nz))^2 +
      ndt_vec_interleaved[!pos]
  }

  # Reshape into a matrix: nobs rows, 2 columns (acc1, acc2)
  finish_time_matrix <- matrix(finish_time_vec, nrow = nobs, ncol = 2, byrow = TRUE)

  # Extract times for each accumulator
  ft_acc1 <- finish_time_matrix[, 1]
  ft_acc2 <- finish_time_matrix[, 2]

  # Find minimum RT using pmin
  rt <- pmin(ft_acc1, ft_acc2)

  # Determine the winner: 0 for the `vzero` accumulator, 1 for `vone`.
  response <- ifelse(ft_acc1 < ft_acc2, 0L, 1L)

  # --- Return Results ---
  data.frame(rt = rt, response = response)
}



#' @rdname rrdm
#' @inheritParams rrt_invgaussian
#' @param response Accumulator whose finishing time is being scored: `0` for the
#'   `vzero` accumulator, `1` for the `vone` accumulator. This gives the
#'   *defective* density `f_response(x) * S_other(x)`, which is what a race
#'   likelihood needs. The default `NULL` instead returns the *marginal* density
#'   of the winning RT, ignoring which accumulator won.
#' @export
drdm <- function(x, response = NULL, vzero, vone, bs, bias, ndt, log = FALSE) {
  n <- length(x)
  params <- .prepare_rrdm(
    n = n,
    vzero = vzero,
    vone = vone,
    bs = bs,
    bias = bias,
    ndt = ndt
  )
  x_vec <- rep(x, length.out = params$ndraws)

  if (!is.null(response)) {
    response <- rep_len(response, params$ndraws)
    if (!all(response %in% c(0, 1))) {
      stop("'response' must contain only 0 (vzero) or 1 (vone).")
    }
    v_win <- ifelse(response == 0, params$vzero, params$vone)
    v_loss <- ifelse(response == 0, params$vone, params$vzero)
    out <- .dwald(
      x_vec,
      v_win,
      params$bs,
      params$bias,
      params$ndt,
      log = TRUE
    ) +
      .swald(
        x_vec,
        v_loss,
        params$bs,
        params$bias,
        params$ndt,
        log.p = TRUE
      )
    return(if (log) out else exp(out))
  }

  # Marginal density of the winning RT: f1 * S2 + f2 * S1.
  logf1 <- .dwald(
    x_vec,
    params$vzero,
    params$bs,
    params$bias,
    params$ndt,
    log = TRUE
  )
  logf2 <- .dwald(
    x_vec,
    params$vone,
    params$bs,
    params$bias,
    params$ndt,
    log = TRUE
  )
  logS1 <- .swald(
    x_vec,
    params$vzero,
    params$bs,
    params$bias,
    params$ndt,
    log.p = TRUE
  )
  logS2 <- .swald(
    x_vec,
    params$vone,
    params$bs,
    params$bias,
    params$ndt,
    log.p = TRUE
  )

  out <- .log_add_exp(logf1 + logS2, logf2 + logS1)
  if (log) out else exp(out)
}


#' @rdname rrdm
#' @param q Vector of quantiles (reaction times).
#' @param lower.tail If `TRUE` (default) return `P(RT <= q)`, otherwise the
#'   survival `P(RT > q)`.
#' @param log.p If `TRUE`, probabilities are returned on the log scale.
#' @details
#' `prdm()` describes the RT of the race as a whole (whichever accumulator
#' wins), since `P(min(T0, T1) > q) = S0(q) * S1(q)`. There is no comparably
#' simple closed form for the per-response defective CDF, so `prdm()` takes no
#' `response` argument.
#' @export
prdm <- function(
  q,
  vzero,
  vone,
  bs,
  bias,
  ndt,
  lower.tail = TRUE,
  log.p = FALSE
) {
  n <- length(q)
  params <- .prepare_rrdm(
    n = n,
    vzero = vzero,
    vone = vone,
    bs = bs,
    bias = bias,
    ndt = ndt
  )
  q_vec <- rep(q, length.out = params$ndraws)

  logS <- .swald(
    q_vec,
    params$vzero,
    params$bs,
    params$bias,
    params$ndt,
    log.p = TRUE
  ) +
    .swald(
      q_vec,
      params$vone,
      params$bs,
      params$bias,
      params$ndt,
      log.p = TRUE
    )

  out <- if (lower.tail) .log1m_exp(logS) else logS
  if (log.p) out else exp(out)
}


# Internals ---------------------------------------------------------------

# Log-space building blocks -----------------------------------------------
#
# The Wald-with-uniform-start-point density and survival are both *difference
# quotients in the threshold* `b`: writing `b = bs + bias`, the start point
# `z ~ U(0, bias)` means
#
#   f_A(t) = (1/bias) * Integral_{bs}^{bs+bias} f_Wald(t; b) db
#   S_A(t) = (1/bias) * Integral_{bs}^{bs+bias} S_Wald(t; b) db
#
# so both take the form `[G(bs + bias) - G(bs)] / bias`. That structure drives
# every numerical choice below:
#
#  * the whole computation is carried in log space, because for moderate drift
#    the individual terms underflow to zero long before the answer does (at
#    drift 3 the survival hits exactly 0 by t = 12, while its true value is
#    4e-26);
#  * `exp(2 * drift * b) * Phi(.)` is evaluated as `exp(2 * drift * b + lcdf)`,
#    since the product is representable when neither factor is;
#  * as `bias / sqrt(t)` shrinks the difference quotient cancels, so we fall
#    back to the plain Wald evaluated at the *midpoint* threshold
#    `b = bs + bias / 2`, which is second-order accurate and therefore exact to
#    double precision well before cancellation sets in.

# Smallest `bias / sqrt(t)` handled by the exact difference; below this the
# midpoint (plain Wald) limit is both cheaper and more accurate.
.RDM_EPS_A <- 1e-4
# Below this |drift| the survival's 1 / (2 * drift) factor is replaced by its
# driftless limit. Measured crossover between the two errors is near 3e-8.
.RDM_EPS_V <- 1e-7


#' log(exp(x) + exp(y)), stable
#' @noRd
.log_add_exp <- function(x, y) {
  m <- pmax(x, y)
  out <- m + log(exp(x - m) + exp(y - m))
  out[is.infinite(m) & m < 0] <- -Inf
  out
}


#' log(exp(x) - exp(y)) for x >= y, stable (reuses `.log1m_exp`)
#'
#' `y >= x` collapses to `-Inf` rather than `NaN`; the guard is applied *before*
#' evaluating `.log1m_exp()` so a masked-out element cannot raise a warning.
#'
#' @noRd
.log_sub_exp <- function(x, y) {
  n <- max(length(x), length(y))
  x <- rep_len(x, n)
  y <- rep_len(y, n)
  out <- rep(-Inf, n)
  bad <- is.na(x) | is.na(y)
  ok <- !bad & (y < x)
  if (any(ok)) {
    out[ok] <- x[ok] + .log1m_exp(y[ok] - x[ok])
  }
  out[bad] <- NA_real_
  out
}


#' Signed log-domain sum
#'
#' Given terms as (log-magnitude, sign) pairs, returns the log-magnitude and
#' sign of their sum. Used to assemble the density (2 terms) and survival
#' (6 terms), each of which mixes signs but is known to be positive overall.
#'
#' @noRd
.log_sum_signed <- function(logs, signs) {
  n <- length(logs[[1]])
  acc_l <- rep(-Inf, n)
  acc_s <- rep(0, n)
  for (i in seq_along(logs)) {
    li <- logs[[i]]
    si <- signs[[i]]
    li[si == 0 | is.na(li)] <- -Inf
    si[is.na(si)] <- 0
    new_l <- rep(-Inf, n)
    new_s <- numeric(n)
    # Same sign (or either term vanishing): magnitudes add.
    same <- (acc_s == si) | (acc_s == 0) | (si == 0)
    if (any(same)) {
      new_l[same] <- .log_add_exp(acc_l[same], li[same])
      new_s[same] <- ifelse(acc_s[same] != 0, acc_s[same], si[same])
    }
    # Opposite signs: subtract the smaller magnitude, keep the larger's sign.
    op <- !same
    if (any(op)) {
      keep_acc <- acc_l[op] >= li[op]
      new_l[op] <- .log_sub_exp(
        pmax(acc_l[op], li[op]),
        pmin(acc_l[op], li[op])
      )
      new_s[op] <- ifelse(keep_acc, acc_s[op], si[op])
    }
    new_s[is.infinite(new_l) & new_l < 0] <- 0
    acc_l <- new_l
    acc_s <- new_s
  }
  list(log = acc_l, sign = acc_s)
}


#' Stable `log(u * Phi(u) + phi(u))`
#'
#' This is `log Integral_{-Inf}^{u} Phi(s) ds`, the antiderivative that appears
#' when integrating the Wald survival over the threshold. For very negative `u`
#' both terms underflow *and* they cancel to order `u^-2`, so we factor out
#' `phi(u)` and expand `1 - |u| * MillsRatio(|u|)` asymptotically.
#'
#' @noRd
.log_gfun <- function(u) {
  out <- numeric(length(u))
  # Seam at -10: the direct form still cancels only ~100-fold there, while the
  # series has converged to ~1e-12, so the two agree to well past what the
  # likelihood needs.
  hi <- u > -10
  if (any(hi)) {
    uh <- u[hi]
    out[hi] <- log(uh * stats::pnorm(uh) + stats::dnorm(uh))
  }
  lo <- !hi
  if (any(lo)) {
    x2 <- u[lo]^2
    s <- 0
    term <- rep(1, length(x2))
    for (j in seq_len(12)) {
      term <- term * (2 * j - 1) / x2
      s <- if (j %% 2 == 1) s + term else s - term
    }
    out[lo] <- stats::dnorm(u[lo], log = TRUE) + log(s)
  }
  out
}


#' Stable `log(Phi(b) - Phi(a))` for `b >= a`
#'
#' Picks whichever tail keeps both arguments away from a saturating normal CDF.
#'
#' @noRd
.log_diff_pnorm <- function(a, b) {
  out <- numeric(length(a))
  up <- a >= 0
  if (any(up)) {
    out[up] <- .log_sub_exp(
      stats::pnorm(a[up], lower.tail = FALSE, log.p = TRUE),
      stats::pnorm(b[up], lower.tail = FALSE, log.p = TRUE)
    )
  }
  dn <- b <= 0
  if (any(dn)) {
    out[dn] <- .log_sub_exp(
      stats::pnorm(b[dn], log.p = TRUE),
      stats::pnorm(a[dn], log.p = TRUE)
    )
  }
  mid <- !up & !dn # straddles zero: both probabilities are O(1)
  if (any(mid)) {
    out[mid] <- log(stats::pnorm(b[mid]) - stats::pnorm(a[mid]))
  }
  out
}


# Wald with uniform start point: log density and log survival --------------

#' Log density of the first passage time `t` (already net of non-decision time)
#'
#' @param t Decision time, `> 0`.
#' @param nu Drift rate.
#' @param k Threshold offset `bs`.
#' @param A Start point range `bias` (`>= 0`).
#' @noRd
.wald_lpdf_core <- function(t, nu, k, A) {
  st <- sqrt(t)

  # Small start-point range (or large t): plain Wald at the midpoint threshold.
  small <- (A / st) < .RDM_EPS_A
  out <- numeric(length(t))
  if (any(small)) {
    b <- k[small] + A[small] / 2
    ts <- t[small]
    out[small] <- log(b) -
      0.5 * (log(2 * pi) + 3 * log(ts)) -
      (b - nu[small] * ts)^2 / (2 * ts)
  }

  gen <- !small
  if (any(gen)) {
    tg <- t[gen]
    stg <- st[gen]
    nug <- nu[gen]
    kg <- k[gen]
    Ag <- A[gen]
    alpha <- (kg - nug * tg) / stg
    beta <- (kg + Ag - nug * tg) / stg

    # T1 = drift * (Phi(beta) - Phi(alpha))
    l1 <- log(abs(nug)) + .log_diff_pnorm(alpha, beta)
    s1 <- sign(nug)

    # T2 = (phi(alpha) - phi(beta)) / sqrt(t); phi is even, so the sign is
    # decided by which argument is closer to zero.
    lpa <- stats::dnorm(alpha, log = TRUE)
    lpb <- stats::dnorm(beta, log = TRUE)
    s2 <- ifelse(abs(alpha) <= abs(beta), 1, -1)
    l2 <- ifelse(
      s2 > 0,
      .log_sub_exp(lpa, lpb),
      .log_sub_exp(lpb, lpa)
    ) -
      log(stg)

    res <- .log_sum_signed(list(l1, l2), list(s1, s2))
    out[gen] <- res$log - log(Ag)
  }
  out
}


#' Log survival of the first passage time `t` (already net of non-decision time)
#'
#' Uses the closed-form antiderivative of the Wald survival in the threshold,
#' rather than `log(1 - CDF)`, which underflows to `-Inf` at ordinary parameter
#' values and would stall the sampler with a zero gradient.
#'
#' @inheritParams .wald_lpdf_core
#' @noRd
.wald_lccdf_core <- function(t, nu, k, A) {
  st <- sqrt(t)
  out <- numeric(length(t))

  # --- Small start-point range: plain Wald survival at the midpoint. ---
  small <- (A / st) < .RDM_EPS_A
  if (any(small)) {
    b <- k[small] + A[small] / 2
    ts <- t[small]
    sts <- st[small]
    nus <- nu[small]
    l1 <- stats::pnorm((b - nus * ts) / sts, log.p = TRUE)
    l2 <- 2 * nus * b + stats::pnorm(-(b + nus * ts) / sts, log.p = TRUE)
    out[small] <- .log_sub_exp(l1, l2)
  }

  gen <- !small
  if (any(gen)) {
    tg <- t[gen]
    stg <- st[gen]
    nug <- nu[gen]
    kg <- k[gen]
    Ag <- A[gen]
    bg <- kg + Ag
    alpha <- (kg - nug * tg) / stg
    beta <- (bg - nug * tg) / stg

    # `res` accumulates log(S * A); the shared 1 / A is divided out at the end.
    res <- rep(NA_real_, length(tg))

    # --- Driftless limit: S = 1 + (2 sqrt(t) / A) * (g(beta0) - g(alpha0)). ---
    # This form already carries its own 1 / A, so log(A) is added back to keep
    # it on the same `log(S * A)` footing as the general branch below.
    z <- abs(nug) < .RDM_EPS_V
    if (any(z)) {
      a0 <- -kg[z] / stg[z]
      b0 <- -bg[z] / stg[z]
      lgap <- .log_sub_exp(.log_gfun(a0), .log_gfun(b0))
      res[z] <- .log_sum_signed(
        list(rep(0, sum(z)), log(2) + log(stg[z]) - log(Ag[z]) + lgap),
        list(rep(1, sum(z)), rep(-1, sum(z)))
      )$log +
        log(Ag[z])
    }

    # --- General case: six signed terms, positive overall. ---
    g <- !z
    if (any(g)) {
      tt <- tg[g]
      stt <- stg[g]
      nn <- nug[g]
      kk <- kg[g]
      AA <- Ag[g]
      bb <- bg[g]
      al <- alpha[g]
      be <- beta[g]
      sv <- sign(nn)
      linv <- -log(2 * abs(nn)) # log(1 / (2 |drift|))
      wa <- -(kk + nn * tt) / stt
      wb <- -(bb + nn * tt) / stt

      logs <- list(
        0.5 * log(tt) + .log_gfun(be), # + sqrt(t) g(beta)
        0.5 * log(tt) + .log_gfun(al), # - sqrt(t) g(alpha)
        linv + 2 * nn * bb + stats::pnorm(wb, log.p = TRUE),
        linv + 2 * nn * kk + stats::pnorm(wa, log.p = TRUE),
        linv + stats::pnorm(be, log.p = TRUE),
        linv + stats::pnorm(al, log.p = TRUE)
      )
      m <- length(tt)
      signs <- list(
        rep(1, m),
        rep(-1, m),
        -sv,
        sv,
        -sv,
        sv
      )
      res[g] <- .log_sum_signed(logs, signs)$log
    }

    out[gen] <- res - log(Ag)
  }
  pmin(out, 0)
}


#' Analytical PDF for Wald with Uniform Start Point Variability
#'
#' Calculates the probability density function (PDF) for the time it takes a
#' diffusion process with drift `drift`, starting point `z ~ U(0, bias)`, to reach
#' threshold `b = bs + bias`, shifted by `ndt`. Based on Tillman et al. (2020).
#' This function calculates the density for the *unadjusted* time `x`.
#'
#' @noRd
.dwald <- function(x, drift, bs, bias, ndt, log = FALSE) {
  n <- length(x)
  if (!all(lengths(list(drift, bs, bias, ndt)) == n)) {
    stop("All inputs must be numeric vectors of the same length.")
  }

  out <- rep(NA_real_, n)
  invalid <- is.na(x) |
    is.na(drift) |
    is.na(bs) |
    is.na(bias) |
    is.na(ndt) |
    (bs <= 0) |
    (bias <= 0) |
    (ndt < 0)

  ok <- !invalid
  if (any(ok)) {
    t_adj <- x[ok] - ndt[ok]
    res <- rep(-Inf, sum(ok)) # zero density at or before the non-decision time
    live <- t_adj > 0
    if (any(live)) {
      res[live] <- .wald_lpdf_core(
        t_adj[live],
        drift[ok][live],
        bs[ok][live],
        bias[ok][live]
      )
    }
    out[ok] <- res
  }
  if (log) out else exp(out)
}


#' Analytical survival function for Wald with Uniform Start Point Variability
#'
#' Complement of `.pwald()`, computed directly rather than as `1 - CDF`. See the
#' note above `.wald_lccdf_core()` for why that distinction matters.
#'
#' @noRd
.swald <- function(x, drift, bs, bias, ndt, log.p = FALSE) {
  n <- length(x)
  if (!all(lengths(list(drift, bs, bias, ndt)) == n)) {
    stop("All inputs must be numeric vectors of the same length.")
  }

  out <- rep(NA_real_, n)
  invalid <- is.na(x) |
    is.na(drift) |
    is.na(bs) |
    is.na(bias) |
    is.na(ndt) |
    (bs <= 0) |
    (bias <= 0) |
    (ndt < 0)

  ok <- !invalid
  if (any(ok)) {
    t_adj <- x[ok] - ndt[ok]
    res <- numeric(sum(ok)) # log(1): nothing can finish before the ndt
    live <- t_adj > 0
    if (any(live)) {
      res[live] <- .wald_lccdf_core(
        t_adj[live],
        drift[ok][live],
        bs[ok][live],
        bias[ok][live]
      )
    }
    out[ok] <- res
  }
  if (log.p) out else exp(out)
}


#' Analytical CDF for Wald with Uniform Start Point Variability
#'
#' Calculates the cumulative distribution function (CDF) for the time it takes a
#' diffusion process with drift `drift`, starting point `z ~ U(0, bias)`, to reach
#' threshold `b = bs + bias`, shifted by `ndt`. Based on Tillman et al. (2020).
#' This function calculates the CDF for the *unadjusted* time `x`.
#' 
#' @noRd
.pwald <- function(
  x,
  drift,
  bs,
  bias,
  ndt,
  lower.tail = TRUE,
  log.p = FALSE
) {
  n <- length(x)
  if (!all(lengths(list(drift, bs, bias, ndt)) == n)) {
    stop("All inputs must be numeric vectors of the same length.")
  }

  # The survival is the primitive: it is the quantity the race likelihood
  # needs, and it is the one that stays accurate far into the upper tail.
  logS <- .swald(x, drift, bs, bias, ndt, log.p = TRUE)

  if (!lower.tail) {
    return(if (log.p) logS else exp(logS))
  }

  # Lower tail. Near t = 0 the CDF is astronomically small and `1 - S` cannot
  # resolve it, so use the direct series there; elsewhere `1 - S` is fine.
  out <- .log1m_exp(logS)
  tiny <- !is.na(logS) & logS > -1e-3
  if (any(tiny)) {
    out[tiny] <- .wald_lcdf_direct(
      x[tiny],
      drift[tiny],
      bs[tiny],
      bias[tiny],
      ndt[tiny]
    )
  }
  out[is.na(logS)] <- NA_real_
  if (log.p) out else exp(out)
}


#' Direct log-CDF, used only in the lower tail where `log(1 - S)` cancels
#'
#' Same closed form as Tillman et al. (2020), with `exp(2 * drift * b) * Phi(.)`
#' folded into a single exponent so it stays finite when the factors separately
#' overflow and underflow.
#'
#' @noRd
.wald_lcdf_direct <- function(x, drift, bs, bias, ndt) {
  t <- x - ndt
  out <- rep(-Inf, length(t))
  live <- !is.na(t) & t > 0
  if (!any(live)) {
    return(out)
  }

  tt <- t[live]
  stt <- sqrt(tt)
  nn <- drift[live]
  kk <- bs[live]
  AA <- bias[live]
  bb <- kk + AA

  # Driftless limit: C = (2 sqrt(t) / A) * (g(-k/sqrt(t)) - g(-b/sqrt(t))).
  # NOTE: the factor of 2 here was missing before; the branch was unreachable
  # from `drdm()` (which short-circuits when a drift is zero) so it went
  # unnoticed, but the Stan version reaches it whenever drift approaches 0.
  res <- rep(NA_real_, length(tt))
  z <- abs(nn) < .RDM_EPS_V
  if (any(z)) {
    lgap <- .log_sub_exp(
      .log_gfun(-kk[z] / stt[z]),
      .log_gfun(-bb[z] / stt[z])
    )
    res[z] <- log(2) + log(stt[z]) - log(AA[z]) + lgap
  }

  g <- !z
  if (any(g)) {
    tg <- tt[g]
    sg <- stt[g]
    ng <- nn[g]
    kg <- kk[g]
    Ag <- AA[g]
    bg <- bb[g]
    sv <- sign(ng)
    linv <- -log(2 * abs(ng))
    a1 <- (ng * tg - bg) / sg
    a2 <- (ng * tg - kg) / sg
    w1 <- (-ng * tg - bg) / sg
    w2 <- (-ng * tg - kg) / sg
    m <- length(tg)

    logs <- list(
      linv + .log_diff_pnorm(a1, a2), #  (1/2vA)(Phi(a2) - Phi(a1))
      log(sg) + .log_gfun(a2), # +sqrt(t)/A * g(a2)
      log(sg) + .log_gfun(a1), # -sqrt(t)/A * g(a1)
      linv + 2 * ng * kg + stats::pnorm(w2, log.p = TRUE),
      linv + 2 * ng * bg + stats::pnorm(w1, log.p = TRUE)
    )
    signs <- list(sv, rep(1, m), rep(-1, m), -sv, sv)
    res[g] <- .log_sum_signed(logs, signs)$log - log(Ag)
  }
  out[live] <- res
  pmin(out, 0)
}



#' @noRd
.prepare_rrdm <- function(n = NULL, vzero, vone, bs, bias, ndt) {
  # Validate parameters
  if (any(vzero < 0)) stop("Drift rate 'vzero' must be non-negative.")
  if (any(vone < 0)) stop("Drift rate 'vone' must be non-negative.")
  if (any(bs <= 0)) stop("Threshold parameter 'bs' must be positive.")
  if (any(bias <= 0)) stop("Maximum starting point 'bias' must be positive.")
  if (any(ndt < 0)) stop("Non-decision time 'ndt' must be non-negative.")

  # Determine target length
  if (!is.null(n)) {
    if (length(n) != 1 || n <= 0 || n != floor(n)) {
      stop("n must be a single positive integer.")
    }
    m <- n
  } else {
    stop("Parameter 'n' must be provided for simulation.")
  }

  # Recycle parameters to length m
  params <- list(
    vzero = rep_len(vzero, m),
    vone  = rep_len(vone,  m),
    bs     = rep_len(bs,     m),
    bias     = rep_len(bias,     m),
    ndt   = rep_len(ndt,   m)
  )

  params$ndraws <- m
  params
}