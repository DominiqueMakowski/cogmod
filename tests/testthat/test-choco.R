context("CHOCO")


test_that("rchoco generates correct proportion of zeros, ones, and threshold values with varying parameters", {
  # Set seed for reproducibility
  set.seed(123)

  # Test parameters
  n_samples <- 30000
  threshold <- 0.5
  tolerance <- 0.015  # Allow for some sampling variability

  # Parameter values to test (using new parameterization)
  p_values <- c(0, 0.5, 1)
  conf_values <- c(0.3, 0.7) # Replaced muleft_values
  confleft_values <- c(-0.5, 0, 0.5) # Replaced mudelta_values
  prec_values <- c(4) # Using a single value for simplicity in this test
  precleft_values <- c(0) # Using a single value for simplicity in this test
  pex_values <- c(0.1, 0.5)
  bex_values <- c(0.2, 0.5, 0.8)
  pmid_values <- c(0, 0.2, 0.5)

  # Loop through all parameter combinations
  for (p in p_values) {
    for (conf in conf_values) {
      for (confleft in confleft_values) {
        for (prec in prec_values) { # Loop added, though only one value
          for (precleft in precleft_values) { # Loop added, though only one value
            for (pex in pex_values) {
              for (bex in bex_values) {
                for (pmid in pmid_values) {
                  # Generate descriptive label for this parameter combination
                  label <- sprintf("p=%.1f, conf=%.1f, confleft=%.1f, prec=%.1f, precleft=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                                   p, conf, confleft, prec, precleft, pex, bex, pmid)

                  # Generate data with current parameter combination (using new names)
                  x <- rchoco(n = n_samples,
                              p = p,
                              conf = conf,
                              confleft = confleft,
                              prec = prec,
                              precleft = precleft,
                              pex = pex,
                              bex = bex,
                              pmid = pmid,
                              threshold = threshold)

                  # Calculate empirical proportions
                  prop_zeros <- sum(abs(x) < 1e-10) / n_samples
                  prop_ones <- sum(abs(x - 1) < 1e-10) / n_samples
                  prop_threshold <- sum(abs(x - threshold) < 1e-10) / n_samples

                  # Expected proportions based on the logic in rchoco/dchoco
                  expected_prop_threshold <- pmid
                  # Calculate effective pex for left (0) and right (1) based on rchoco logic
                  pex_left_eff <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
                  pex_right_eff <- pmin(1, pmax(0, bex * (pex * 2)))
                  # Expected proportions for 0 and 1 depend on p and the effective pex values
                  expected_prop_zeros <- (1 - pmid) * (1 - p) * pex_left_eff
                  expected_prop_ones <- (1 - pmid) * p * pex_right_eff

                  # Test proportions with informative labels
                  expect_equal(prop_threshold, expected_prop_threshold,
                              tolerance = tolerance,
                              label = paste0("Proportion of threshold values with ", label))

                  expect_equal(prop_zeros, expected_prop_zeros,
                              tolerance = tolerance,
                              label = paste0("Proportion of zeros with ", label))

                  expect_equal(prop_ones, expected_prop_ones,
                              tolerance = tolerance,
                              label = paste0("Proportion of ones with ", label))

                  # Also test total proportion of extremes
                  expect_equal(prop_zeros + prop_ones, expected_prop_zeros + expected_prop_ones,
                              tolerance = tolerance, # Use same tolerance for sum
                              label = paste0("Total proportion of extremes with ", label))
                } # pmid
              } # bex
            } # pex
          } # precleft
        } # prec
      } # confleft
    } # conf
  } # p
})



test_that("dchoco matches rchoco empirical distribution and integrates correctly", { # Renamed rbext to rchoco
  # Set seed for reproducibility
  set.seed(456)

  # Test parameters
  n_samples <- 15000 # Increased samples for better empirical estimates
  threshold <- 0.5
  tol_prob <- 0.025 # Tolerance for comparing probabilities (empirical vs theoretical)
  tol_integration <- 0.01 # Tolerance for numerical integration check
  tol_ks_p <- 0.001 # Minimum p-value for KS test to pass

  # Reduced parameter set for faster testing, focusing on variety (using new names)
  p_values <- c(0.2, 0.5, 0.8)
  conf_values <- c(0.3, 0.7) # Replaced muleft_values
  confleft_values <- c(-0.5, 0.5) # Replaced mudelta_values
  prec_values <- c(2, 5) # Replaced phileft_values
  precleft_values <- c(-0.5, 0.5) # Replaced phidelta_values
  pex_values <- c(0, 0.2, 0.6)
  bex_values <- c(0.2, 0.8)
  pmid_values <- c(0, 0.3)


  # Loop through parameter combinations
  for (p in p_values) {
    for (conf in conf_values) {
      for (confleft in confleft_values) {
        for (prec in prec_values) {
          for (precleft in precleft_values) {
            for (pex in pex_values) {
              for (bex in bex_values) {
                for (pmid in pmid_values) {

                  # Generate descriptive label (using new names)
                  label <- sprintf("p=%.1f, conf=%.1f, confL=%.1f, prec=%.1f, precL=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                                   p, conf, confleft, prec, precleft, pex, bex, pmid)

                  # Generate data (using new names)
                  x_sample <- rchoco(n = n_samples, p = p, conf = conf, confleft = confleft,
                                      prec = prec, precleft = precleft, pex = pex,
                                      bex = bex, pmid = pmid, threshold = threshold)

                  # Calculate empirical proportions
                  emp_p0 <- mean(abs(x_sample) < 1e-10)
                  emp_p1 <- mean(abs(x_sample - 1) < 1e-10)
                  emp_pth <- mean(abs(x_sample - threshold) < 1e-10)
                  emp_pleft_cont <- mean(x_sample > 1e-10 & x_sample < threshold - 1e-10)
                  emp_pright_cont <- mean(x_sample > threshold + 1e-10 & x_sample < 1 - 1e-10)

                  # Calculate theoretical probabilities at point masses (using new names)
                  theo_p0 <- dchoco(0, p = p, conf = conf, confleft = confleft, prec = prec,
                                     precleft = precleft, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  theo_p1 <- dchoco(1, p = p, conf = conf, confleft = confleft, prec = prec,
                                     precleft = precleft, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  theo_pth <- dchoco(threshold, p = p, conf = conf, confleft = confleft, prec = prec,
                                      precleft = precleft, pex = pex, bex = bex, pmid = pmid, threshold = threshold)

                  # --- Check Point Masses ---
                  expect_equal(emp_p0, theo_p0, tolerance = tol_prob, label = paste(label, "- P(0)"))
                  expect_equal(emp_p1, theo_p1, tolerance = tol_prob, label = paste(label, "- P(1)"))
                  expect_equal(emp_pth, theo_pth, tolerance = tol_prob, label = paste(label, "- P(threshold)"))

                  # --- Integration Check ---
                  integrand <- function(x_int) {
                    dchoco(x_int, p = p, conf = conf, confleft = confleft, prec = prec,
                            precleft = precleft, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  }

                  # Integrate continuous parts (handle potential errors)
                  integral_left <- tryCatch(
                    stats::integrate(integrand, lower = .Machine$double.eps, upper = threshold - .Machine$double.eps,
                                     subdivisions = 200, stop.on.error = FALSE)$value,
                    error = function(e) 0 # Assume 0 if integration fails
                  )
                  integral_right <- tryCatch(
                    stats::integrate(integrand, lower = threshold + .Machine$double.eps, upper = 1 - .Machine$double.eps,
                                     subdivisions = 200, stop.on.error = FALSE)$value,
                    error = function(e) 0
                  )

                  # Check total probability sums to 1
                  total_prob <- theo_p0 + theo_p1 + theo_pth + integral_left + integral_right
                  expect_equal(total_prob, 1, tolerance = tol_integration, label = paste(label, "- Total Probability Integration"))

                  # --- Mean Check for Continuous Parts ---
                  # Calculate parameters needed for underlying Beta distributions (using new logic)
                  eps_mu <- 1e-9
                  muright <- pmax(eps_mu, pmin(conf, 1 - eps_mu)) # Clamp conf directly
                  phiright <- pmax(eps_mu, prec) # Clamp prec directly

                  logit_conf <- log(conf / (1 - conf))
                  logit_muleft <- logit_conf + confleft
                  muleft <- exp(logit_muleft) / (1 + exp(logit_muleft))
                  muleft_clamped <- pmax(eps_mu, pmin(muleft, 1 - eps_mu)) # Clamp derived muleft
                  phileft <- prec * exp(precleft)
                  phileft_clamped <- pmax(eps_mu, phileft) # Clamp derived phileft

                  # Effective pex values needed to check if continuous part exists
                  pex_left_eff <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
                  pex_right_eff <- pmin(1, pmax(0, bex * (pex * 2)))

                  # Theoretical means of the underlying Beta distributions
                  # Mean = shape1 / (shape1 + shape2) = mu * phi * 2 / (phi * 2) = mu
                  theo_mean_left <- muleft_clamped
                  theo_mean_right <- muright # Use the clamped value

                  # Mean Check for Left Continuous Part
                  if (emp_pleft_cont > 0.01 && pex_left_eff < 1) { # Need sufficient continuous data
                    x_left_cont <- x_sample[x_sample > 1e-10 & x_sample < threshold - 1e-10]
                    if (length(x_left_cont) > 10) {
                      y_raw_left <- 1 - x_left_cont / threshold # Transform back to [0, 1]
                      emp_mean_left <- mean(y_raw_left)
                      expect_equal(emp_mean_left, theo_mean_left, tolerance = 0.05, # Increased tolerance for mean
                                   label = paste(label, "- Mean Check Left Continuous"))
                    }
                  }

                  # Mean Check for Right Continuous Part
                  if (emp_pright_cont > 0.01 && pex_right_eff < 1) { # Need sufficient continuous data
                    x_right_cont <- x_sample[x_sample > threshold + 1e-10 & x_sample < 1 - 1e-10]
                     if (length(x_right_cont) > 10) {
                      y_raw_right <- (x_right_cont - threshold) / (1 - threshold) # Transform back to [0, 1]
                      emp_mean_right <- mean(y_raw_right)
                      expect_equal(emp_mean_right, theo_mean_right, tolerance = 0.05, # Increased tolerance for mean
                                   label = paste(label, "- Mean Check Right Continuous"))
                    }
                  }
                } # pmid
              } # bex
            } # pex
          } # precleft
        } # prec
      } # confleft
    } # conf
  } # p
})



context("CHOCO - brms")

test_that("choco model can recover parameters with brms using variational inference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # --- 1. Set up Simulation ---
  # True parameters using the new parameterization
  true_p <- 0.6 # Corresponds to 'mu' in brms
  true_conf <- 0.7
  true_confleft <- -0.5
  true_prec <- 5
  true_precleft <- 0.2
  true_pex <- 0.1
  true_bex <- 0.7
  n_obs <- 3000

  # Generate synthetic data using new parameter names
  set.seed(456)
  df <- data.frame(
    y = rchoco(n = n_obs, p = true_p, conf = true_conf, confleft = true_confleft,
                prec = true_prec, precleft = true_precleft, pex = true_pex,
                bex = true_bex, pmid = 0, threshold = 0.5)
  )

  # --- 2. Define and Fit brms Model ---
  # Formula using new parameter names
  f <- brms::bf(
    y ~ 1, # Intercept for 'mu' (p)
    conf ~ 1,
    confleft ~ 1,
    prec ~ 1,
    precleft ~ 1,
    pex ~ 1,
    bex ~ 1,
    pmid = 0, # Fixed pmid
    family = choco() # Use the custom family
  )

  # Fit model using Pathfinder VI
  m <- brms::brm(
    formula = f,
    data = df,
    stanvars = choco_stanvars(), # Include Stan functions
    init =  0,
    seed = 456,
    refresh = 0,
    algorithm = "pathfinder",
    backend = "cmdstanr"
  )

  # --- 3. Check Parameter Recovery ---
  post_summary <- brms::posterior_summary(m, probs = c(0.05, 0.95))

  # Apply inverse-link functions to recover parameters using new names
  post_p <- brms::inv_logit_scaled(post_summary["b_Intercept", "Estimate"])
  post_conf <- brms::inv_logit_scaled(post_summary["b_conf_Intercept", "Estimate"])
  post_confleft <- post_summary["b_confleft_Intercept", "Estimate"] # identity link
  post_prec <- log(1 + exp(post_summary["b_prec_Intercept", "Estimate"])) # softplus
  post_precleft <- post_summary["b_precleft_Intercept", "Estimate"] # identity link
  post_pex <- brms::inv_logit_scaled(post_summary["b_pex_Intercept", "Estimate"])
  post_bex <- brms::inv_logit_scaled(post_summary["b_bex_Intercept", "Estimate"])

  # Check if posterior means are close to true values
  expect_equal(post_p,  true_p, tolerance = 0.1,
               label = sprintf("Posterior mean of p (%.3f) is close to true p (%.3f)", post_p, true_p))
  expect_equal(post_conf, true_conf, tolerance = 0.1,
               label = sprintf("Posterior mean of conf (%.3f) is close to true conf (%.3f)", post_conf, true_conf))
  expect_equal(post_confleft, true_confleft, tolerance = 0.2, # Allow slightly higher tolerance for identity links
               label = sprintf("Posterior mean of confleft (%.3f) is close to true confleft (%.3f)", post_confleft, true_confleft))
  expect_equal(post_prec, true_prec, tolerance = 1, # Allow higher tolerance for precision
               label = sprintf("Posterior mean of prec (%.3f) is close to true prec (%.3f)", post_prec, true_prec))
  expect_equal(post_precleft, true_precleft, tolerance = 0.3, # Allow slightly higher tolerance for identity links
               label = sprintf("Posterior mean of precleft (%.3f) is close to true precleft (%.3f)", post_precleft, true_precleft))
  expect_equal(post_pex, true_pex, tolerance = 0.1,
               label = sprintf("Posterior mean of pex (%.3f) is close to true pex (%.3f)", post_pex, true_pex))
  expect_equal(post_bex, true_bex, tolerance = 0.1,
               label = sprintf("Posterior mean of bex (%.3f) is close to true bex (%.3f)", post_bex, true_bex))


  # --- 4. Test Post-processing Functions ---
  n_pred_draws <- 10
  pred <- brms::posterior_predict(m, ndraws = n_pred_draws)
  expect_equal(nrow(pred), n_pred_draws)
  expect_equal(ncol(pred), n_obs)
  expect_true(all(pred >= 0 & pred <= 1), "Posterior predictions outside [0, 1]")
  expect_false(any(is.na(pred)), "NA values in posterior predictions")

  n_ll_draws <- 5
  ll <- brms::log_lik(m, ndraws = n_ll_draws)
  expect_equal(nrow(ll), n_ll_draws)
  expect_equal(ncol(ll), n_obs)
  expect_true(all(is.finite(ll)), "Non-finite values found in log-likelihood")
})




test_that("Stan choco_lpdf matches R dchoco function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  choco_lpdf_stan <- choco_lpdf_expose()

  # --- Define parameter grids using new names ---
  y_values <- c(0, 0.01, 0.49, 0.5, 0.51, 0.99, 1) # Test around threshold
  p_values <- c(0.1, 0.5, 0.9) # Corresponds to 'mu' in Stan function
  conf_values <- c(0.2, 0.7)
  confleft_values <- c(-0.5, 0.5)
  prec_values <- c(2, 6)
  precleft_values <- c(-0.3, 0.3)
  pex_values <- c(0, 0.2, 0.8, 1)
  bex_values <- c(0, 0.5, 1)
  pmid_values <- c(0, 0.1, 0.6, 1)

  # --- Loop through parameter combinations ---
  for (p in p_values) {
    for (conf in conf_values) {
      for (confleft in confleft_values) {
        for (prec in prec_values) {
          for (precleft in precleft_values) {
            for (pex in pex_values) {
              for (bex in bex_values) {
                for (pmid in pmid_values) {
                  # Test across different y values
                  for (y in y_values) {
                    label <- sprintf(
                      "y=%.2f, p=%.1f, conf=%.1f, confleft=%.1f, prec=%.1f, precleft=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                      y, p, conf, confleft, prec, precleft, pex, bex, pmid
                    )

                    # Calculate log-density using Stan function (pass new parameters)
                    # Order: y, mu(p), conf, confleft, prec, precleft, pex, bex, pmid
                    s <- capture.output(stan_log_lik <- choco_lpdf_stan(y, p, conf, confleft, prec, precleft, pex, bex, pmid))

                    # Calculate log-density using R function (pass new parameters)
                    r_log_lik <- dchoco(y, p = p, conf = conf, confleft = confleft, prec = prec,
                                        precleft = precleft, pex = pex, bex = bex, pmid = pmid,
                                        threshold = 0.5, log = TRUE)

                    # Compare log-likelihoods
                    expect_equal(stan_log_lik, r_log_lik,
                      tolerance = 1e-6,
                      label = paste("Log-likelihood comparison:", label)
                    )
                  } # y
                } # pmid
              } # bex
            } # pex
          } # precleft
        } # prec
      } # confleft
    } # conf
  } # p

  # --- Test invalid parameter handling using new names ---
  valid_params <- list(x=0.2, p=0.5, conf=0.7, confleft=0, prec=4, precleft=0, pex=0.1, bex=0.5, pmid=0.1, threshold=0.5) # Added threshold

  # R function errors (changed from expect_warning)
  expect_error(do.call(dchoco, c(valid_params[-2], list(p=-0.1))), "p must be between 0 and 1")
  expect_error(do.call(dchoco, c(valid_params[-3], list(conf=1.1))), "conf must be between 0 and 1 \\(exclusive\\)") # Escaped parentheses
  expect_error(do.call(dchoco, c(valid_params[-5], list(prec=-1))), "prec must be positive")

  # Stan function returns -Inf (these remain the same)
  expect_equal(choco_lpdf_stan(y=0.2, mu=-0.1, conf=0.7, confleft=0, prec=4, precleft=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid p (mu)")
  expect_equal(choco_lpdf_stan(y=0.2, mu=0.5, conf=1.1, confleft=0, prec=4, precleft=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid conf")
  expect_equal(choco_lpdf_stan(y=0.2, mu=0.5, conf=0.7, confleft=0, prec=-1, precleft=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid prec")
  expect_equal(choco_lpdf_stan(y=-0.1, mu=0.5, conf=0.7, confleft=0, prec=4, precleft=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid y (<0)")
  expect_equal(choco_lpdf_stan(y=1.1, mu=0.5, conf=0.7, confleft=0, prec=4, precleft=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid y (>1)")

})
