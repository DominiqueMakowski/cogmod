# Stanvars ----------------------------------------------------------------

#' @keywords internal
.conf_sdt_lpdf <- function() {
  "
// Log-likelihood for the reparameterized confidence SDT model.
// Y: observed confidence rating (0, 1, 2, ...).
real conf_sdt_lpdf(real Y, real mu, real c, real tzeroone, real tzerotwo, real toneone, real tonetwo, int dec, int truth) {
  real lower_bound;
  real upper_bound;
  real decision_criterion = -c;
  real mean_x = (truth * 2 - 1) * mu / 2.0;
  real eps = 1e-6; // Small epsilon to avoid numerical issues

  // --- Reparameterization ---
  // The input parameters are positive offsets. Calculate actual criteria.
  real crit_toneone = decision_criterion + toneone;
  real crit_tonetwo = crit_toneone + tonetwo;
  real crit_tzeroone = decision_criterion - tzeroone;
  real crit_tzerotwo = crit_tzeroone - tzerotwo;

  // Input validation is now implicit in the parameterization.

  if (dec == 1) {
    if (Y <= eps) {
      lower_bound = decision_criterion;
      upper_bound = crit_toneone;
    } else if (Y - 1.0 <= eps) {
      lower_bound = crit_toneone;
      upper_bound = crit_tonetwo;
    } else if (Y - 2.0 <= eps) {
      lower_bound = crit_tonetwo;
      upper_bound = positive_infinity();
    } else {
      print(1);
      return negative_infinity(); // Invalid confidence level
    }
  } else if (dec == 0) {
    if (Y <= eps) {
      lower_bound = crit_tzeroone;
      upper_bound = decision_criterion;
    } else if (Y - 1.0 <= eps) {
      lower_bound = crit_tzerotwo;
      upper_bound = crit_tzeroone;
    } else if (Y - 2.0 <= eps) {
      lower_bound = negative_infinity();
      upper_bound = crit_tzerotwo;
    } else {
      print(2);
      return negative_infinity(); // Invalid confidence level
    }
  } else {
    print(3);
    return negative_infinity(); // Invalid response
  }

  // print(\"upper_bound:\", upper_bound);
  // print(\"lower_bound:\", lower_bound);
  // print(\"mean_x:\", mean_x);
  // print(\"upper\", normal_lcdf(upper_bound | mean_x, 1.0));
  // print(\"lower\", normal_lcdf(lower_bound | mean_x, 1.0));
  print(\"logdiffexp\",
        log_diff_exp(normal_lcdf(upper_bound | mean_x, 1.0),
                     normal_lcdf(lower_bound | mean_x, 1.0)));

  // --- Log-likelihood Calculation ---
  return log_diff_exp(normal_lcdf(upper_bound | mean_x, 1.0),
                      normal_lcdf(lower_bound | mean_x, 1.0));
}
"
}

#' @rdname conf_sdt
#' @export
conf_sdt_stanvars <- function() {
  brms::stanvar(scode = .conf_sdt_lpdf(), block = "functions")
}


#' @rdname conf_sdt
#' @param link_mu Link function for the dprime parameter.
#' @param link_c Link function for the c parameter.
#' @param link_tzeroone,link_tzerotwo,link_toneone,link_tonetwo Link functions for the criteria offset parameters.
#' @export
conf_sdt_custom_family <- function(link_mu = "identity", link_c = "identity",
                                   link_tzeroone = "log", link_tzerotwo = "log",
                                   link_toneone = "log", link_tonetwo = "log") {
  brms::custom_family(
    name = "conf_sdt",
    dpars = c("mu", "c", "tzeroone", "tzerotwo", "toneone", "tonetwo"),
    links = c(link_mu, link_c, link_tzeroone, link_tzerotwo, link_toneone, link_tonetwo),
    lb = c(NA, NA, 0, 0, 0, 0), # Lower bounds for positive offsets
    ub = c(NA, NA, NA, NA, NA, NA),
    type = "real",
    vars = c("dec[n]", "vint1[n]")
  )
}

# brms --------------------------------------------------------------------

#' @rdname conf_sdt
#' @export
log_lik_conf_sdt <- function(i, prep) {
  y <- prep$data$Y[i]
  if (is.na(y)) {
    return(NA_real_)
  }

  # Get parameters (offsets for criteria)
  dprime <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  tzeroone_off <- brms::get_dpar(prep, "tzeroone", i = i)
  tzerotwo_off <- brms::get_dpar(prep, "tzerotwo", i = i)
  toneone_off <- brms::get_dpar(prep, "toneone", i = i)
  tonetwo_off <- brms::get_dpar(prep, "tonetwo", i = i)

  truth <- prep$data$truth[i]
  response <- prep$data$dec[i]

  # Reconstruct actual criteria from offsets
  crit_toneone <- -c + toneone_off
  crit_tonetwo <- crit_toneone + tonetwo_off
  crit_tzeroone <- -c - tzeroone_off
  crit_tzerotwo <- crit_tzeroone - tzerotwo_off

  # Call the R density function
  dconf_sdt(
    truth = truth,
    response = response,
    confidence = y,
    dprime = dprime,
    c = c,
    thetazero = cbind(crit_tzeroone, crit_tzerotwo),
    thetaone = cbind(crit_toneone, crit_tonetwo),
    log = TRUE
  )
}


#' @rdname conf_sdt
#' @importFrom brms get_dpar
#' @inheritParams rlnr
#' @export
posterior_predict_conf_sdt <- function(i, prep, ...) {
  # Get parameters (offsets for criteria)
  dprime <- brms::get_dpar(prep, "mu", i = i)
  c <- brms::get_dpar(prep, "c", i = i)
  tzeroone_off <- brms::get_dpar(prep, "tzeroone", i = i)
  tzerotwo_off <- brms::get_dpar(prep, "tzerotwo", i = i)
  toneone_off <- brms::get_dpar(prep, "toneone", i = i)
  tonetwo_off <- brms::get_dpar(prep, "tonetwo", i = i)

  truth <- prep$data$truth[i]

  # Reconstruct actual criteria from offsets
  crit_toneone <- -c + toneone_off
  crit_tonetwo <- crit_toneone + tonetwo_off
  crit_tzeroone <- -c - tzeroone_off
  crit_tzerotwo <- crit_tzeroone - tzerotwo_off

  # Loop over draws because rconf_sdt is not vectorized over the structure of criteria
  n_draws <- length(dprime)
  predictions <- numeric(n_draws)
  for (j in 1:n_draws) {
    sim <- rconf_sdt(
      n = 1,
      dprime = dprime[j],
      c = c[j],
      thetazero = c(crit_tzeroone[j], crit_tzerotwo[j]),
      thetaone = c(crit_toneone[j], crit_tonetwo[j]),
      truth = truth
    )
    if (nrow(sim) > 0) {
      predictions[j] <- sim$confidence
    } else {
      predictions[j] <- NA
    }
  }
  predictions
}


#' @rdname conf_sdt
#' @export
posterior_epred_conf_sdt <- function(prep) {
  # Get parameters (offsets for criteria)
  dprime <- brms::get_dpar(prep, "mu")
  c <- brms::get_dpar(prep, "c")
  tzeroone_off <- brms::get_dpar(prep, "tzeroone")
  tzerotwo_off <- brms::get_dpar(prep, "tzerotwo")
  toneone_off <- brms::get_dpar(prep, "toneone")
  tonetwo_off <- brms::get_dpar(prep, "tonetwo")

  truth <- prep$data$truth
  response <- prep$data$dec

  # Reconstruct actual criteria from offsets
  crit_toneone <- -c + toneone_off
  crit_tonetwo <- crit_toneone + tonetwo_off
  crit_tzeroone <- -c - tzeroone_off
  crit_tzerotwo <- crit_tzeroone - tzerotwo_off

  n_draws <- nrow(dprime)
  n_obs <- ncol(dprime)
  epred <- matrix(NA, nrow = n_draws, ncol = n_obs)
  conf_levels <- 0:2

  for (i in 1:n_obs) {
    probs_k_list <- lapply(conf_levels, function(k) {
      sapply(1:n_draws, function(j) {
        dconf_sdt(
          truth = truth[i],
          response = response[i],
          confidence = k,
          dprime = dprime[j, i],
          c = c[j, i],
          thetazero = c(crit_tzeroone[j, i], crit_tzerotwo[j, i]),
          thetaone = c(crit_toneone[j, i], crit_tonetwo[j, i])
        )
      })
    })
    probs_k <- do.call(cbind, probs_k_list)
    epred[, i] <- probs_k %*% conf_levels
  }

  epred
}
