#' @keywords internal
.ddm_lpdf <- function() {
  "
// Reparameterized Wiener diffusion log-PDF for a single response, using
// Stan's `wiener_lpdf()` (7-parameter form, aka `wiener_full_lpdf()`).
// mu: drift rate parameter
// bs: boundary separation parameter > 0
// bias: initial bias parameter in [0, 1]
// tau: non-decision time proportion in [0, 1]
// minrt: minimum reaction time > 0
// sigmadrift: inter-trial variability in drift (sv), >= 0
// sigmabias: fraction (in [0, 1)) of the maximum allowed starting-point range;
//            sw = sigmabias * fmin(2*bias, 2*(1-bias)) keeps bias +/- sw/2 inside (0, 1)
// sigmatau: scales minrt to yield the non-decision time range st0 (convenience only,
//           not a hard constraint on validity)
//
// Performance note: Stan's own `sw == 0 && st0 == 0` fast path inside
// `wiener_lpdf()`/`wiener_full_lpdf()` still delegates to the newer 'wiener5'
// algorithm (which also handles sv analytically), which is a structurally
// different and empirically much costlier implementation than the classic
// 4-parameter `wiener_lpdf()`, even when sv = 0. We therefore short-circuit
// ourselves whenever sw and st0 both vanish: to the classic 4-parameter
// density when sv is also 0 (fastest path, identical to the original
// `ddm()`), or to the dedicated 5-parameter (sv-only) form otherwise, which
// is still much cheaper than the general 7-parameter form (the latter falls
// back to adaptive numerical quadrature whenever sw or st0 is nonzero).
real ddm_lpdf(real y, real mu, real bs, real bias,
              real tau, real minrt,
              real sigmadrift, real sigmabias, real sigmatau,
              int dec) {
    real ndt = tau * minrt;
    real sw  = sigmabias * fmin(2 * bias, 2 * (1 - bias));
    real st0 = sigmatau * minrt;

    if (sw == 0 && st0 == 0) {
        if (sigmadrift == 0) {
            // Fastest path: classic 4-parameter Navarro & Fuss density.
            if (dec == 1) {
                return wiener_lpdf(y | bs, ndt, bias, mu);
            } else {
                return wiener_lpdf(y | bs, ndt, 1 - bias, -mu);
            }
        } else {
            // Drift-variability-only: dedicated 5-parameter density.
            if (dec == 1) {
                return wiener_lpdf(y | bs, ndt, bias, mu, sigmadrift);
            } else {
                return wiener_lpdf(y | bs, ndt, 1 - bias, -mu, sigmadrift);
            }
        }
    }

    // General case: full 7-parameter Wiener density (adaptive quadrature).
    if (dec == 1) {
        return wiener_lpdf(y | bs, ndt, bias, mu, sigmadrift, sw, st0);
    } else {
        return wiener_lpdf(y | bs, ndt, 1 - bias, -mu, sigmadrift, sw, st0);
    }
}
"
}


#' @rdname rddm
#' @export
ddm_lpdf_expose <- function() {
  insight::check_if_installed("cmdstanr")

  # Wrap the function Stan block (done normally by brms on model compilation)
  stancode <- paste0(
    "functions {
",
    .ddm_lpdf(),
    "
}"
  )

  mod <- cmdstanr::cmdstan_model(
    cmdstanr::write_stan_file(stancode),
    force_recompile = TRUE
  )
  mod$expose_functions()
  mod$functions$ddm_lpdf
}


#' @rdname rddm
#' @export
ddm_stanvars <- function() {
  brms::stanvar(scode = .ddm_lpdf(), block = "functions")
}


#' @rdname rddm
#' @param link_mu,link_bs,link_bias,link_tau,link_minrt Link functions for the parameters.
#' @param link_sigmadrift,link_sigmabias,link_sigmatau Link functions for the extra
#'   inter-trial variability parameters (fix them to 0, e.g. `sigmadrift = 0`, in the
#'   `brms::bf()` formula to recover the classic 4-parameter DDM).
#' @export
ddm <- function(
  link_mu = "identity",
  link_bs = "softplus",
  link_bias = "logit",
  link_tau = "logit",
  link_minrt = "identity",
  link_sigmadrift = "softplus",
  link_sigmabias = "logit",
  link_sigmatau = "logit"
) {
  brms::custom_family(
    name = "ddm",
    dpars = c(
      "mu",
      "bs",
      "bias",
      "tau",
      "minrt",
      "sigmadrift",
      "sigmabias",
      "sigmatau"
    ),
    links = c(
      link_mu,
      link_bs,
      link_bias,
      link_tau,
      link_minrt,
      link_sigmadrift,
      link_sigmabias,
      link_sigmatau
    ),
    lb = c(NA, 0, 0, 0, 0, 0, 0, 0),
    ub = c(NA, NA, 1, 1, NA, NA, 1, 1),
    type = "real", # Response variable type
    vars = "dec[n]" # Required additional data variable for the decision
  )
}


# brms --------------------------------------------------------------------

#' @rdname rddm
#' @inheritParams rlnr
#' @export
log_lik_ddm <- function(i, prep, ...) {
  # Extract parameters for observation i across all posterior draws (vectors)
  drift <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i) # Vector if varying, scalar if fixed
  sigmadrift <- brms::get_dpar(prep, "sigmadrift", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  sigmatau <- brms::get_dpar(prep, "sigmatau", i = i)

  # Extract observed data for observation i (scalars)
  y <- prep$data$Y[i] # Observed RT
  dec <- prep$data$dec[i] # Observed decision

  # Handle missing data for the observation
  if (is.na(y) || is.na(dec)) {
    # Return vector of NAs, one for each draw
    return(rep(NA_real_, prep$ndraws))
  }

  # Basic check for valid response coding (0 or 1)
  if (!dec %in% c(0, 1)) {
    warning(
      "Response ('dec') must be 0 or 1. Found: ",
      dec,
      " for observation ",
      i,
      ". Returning -Inf log-likelihood."
    )
    # Return -Inf for all draws for this observation
    return(rep(-Inf, prep$ndraws))
  }

  # Calculate ndt (vector, one element per draw)
  # Element-wise multiplication handles if minrt is scalar or vector
  ndt <- tau * minrt

  # Compute log-likelihood using dddm(), which automatically dispatches to
  # the full 7-parameter density (via `rtdists`) whenever sigmadrift,
  # sigmabias, or sigmatau is non-zero, and to the fast 4-parameter
  # brms::dwiener() otherwise.
  ll <- dddm(
    x = y,
    drift = drift,
    bs = bs,
    bias = bias,
    ndt = ndt,
    response = dec,
    log = TRUE,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt,
    ...
  )

  # Return vector of log-likelihoods (one per posterior draw)
  ll
}


#' @rdname rddm
#' @inheritParams rlnr
#' @export
posterior_epred_ddm <- function(prep) {
  # Extract parameters
  mu <- prep$dpars$mu
  bs <- prep$dpars$bs
  bias <- prep$dpars$bias
  tau <- prep$dpars$tau
  minrt <- prep$dpars$minrt
  sigmadrift <- prep$dpars$sigmadrift
  sigmabias <- prep$dpars$sigmabias
  sigmatau <- prep$dpars$sigmatau

  # The closed-form E[RT] formula below is only exact for the classic
  # 4-parameter DDM. Warn (once) when inter-trial variability is present,
  # since there is no simple closed form for E[RT] in that case: users who
  # need exact predictions should use posterior_predict() and average the
  # simulated draws instead.
  if (any(sigmadrift != 0) || any(sigmabias != 0) || any(sigmatau != 0)) {
    warning(
      "posterior_epred() for `ddm()` uses a closed-form approximation that is only exact ",
      "for the classic 4-parameter DDM (sigmadrift = sigmabias = sigmatau = 0). Since ",
      "inter-trial variability parameters are non-zero here, these expected-RT predictions ",
      "are only approximate. For exact predictions, use `posterior_predict()` and average ",
      "the simulated draws instead.",
      call. = FALSE
    )
  }

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Compute expected reaction time (E[RT])
  # Formula adapted from https://doi.org/10.1016/j.jmp.2009.01.006
  epred <- ndt -
    bias / mu +
    bs / mu * (exp(-2 * mu * bias) - 1) / (exp(-2 * mu * bs) - 1)

  # Return expected reaction time
  epred
}


#' @rdname rddm
#' @inheritParams rlnr
#' @export
posterior_predict_ddm <- function(i, prep, ...) {
  # Extract parameters for observation i
  mu <- brms::get_dpar(prep, "mu", i = i)
  bs <- brms::get_dpar(prep, "bs", i = i)
  bias <- brms::get_dpar(prep, "bias", i = i)
  tau <- brms::get_dpar(prep, "tau", i = i)
  minrt <- brms::get_dpar(prep, "minrt", i = i)
  sigmadrift <- brms::get_dpar(prep, "sigmadrift", i = i)
  sigmabias <- brms::get_dpar(prep, "sigmabias", i = i)
  sigmatau <- brms::get_dpar(prep, "sigmatau", i = i)

  # Compute ndt from tau and minrt
  ndt <- tau * minrt

  # Simulate reaction times and responses using rddm(), which automatically
  # dispatches to the full 7-parameter simulation (via `rtdists`) whenever
  # sigmadrift, sigmabias, or sigmatau is non-zero.
  sim_data <- rddm(
    n = prep$ndraws,
    drift = mu,
    bs = bs,
    bias = bias,
    ndt = ndt,
    sigmadrift = sigmadrift,
    sigmabias = sigmabias,
    sigmatau = sigmatau,
    minrt = minrt
  )

  # Return simulated reaction times
  as.matrix(sim_data)
}
