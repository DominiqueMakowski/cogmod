context("CHOCO")


test_that("rchoco generates correct proportion of zeros, ones, and threshold values with varying parameters", {
  # Set seed for reproducibility
  set.seed(123)

  # Test parameters
  n_samples <- 30000
  threshold <- 0.5
  tolerance <- 0.015  # Allow for some sampling variability

  # Parameter values to test
  p_values <- c(0, 0.5, 1)
  # p_values <- c(0.5)
  pex_values <- c(0.1, 0.5)
  muleft_values <- c(0.3, 0.7)
  mudelta_values <- c(-0.5, 0, 0.5)
  bex_values <- c(0.2, 0.5, 0.8)
  pmid_values <- c(0, 0.2, 0.5)

  # Loop through all parameter combinations
  for (p in p_values) {
    for (pex in pex_values) {
      for (muleft in muleft_values) {
        for (mudelta in mudelta_values) {
          for (bex in bex_values) {
            for (pmid in pmid_values) {
              # Generate descriptive label for this parameter combination
              label <- sprintf("p=%.1f, muleft=%.1f, mudelta=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                               p, muleft, mudelta, pex, bex, pmid)

              # Generate data with current parameter combination
              x <- rchoco(n = n_samples,
                          p = p,
                          muleft = muleft,
                          mudelta = mudelta,
                          phileft = 4,
                          phidelta = 0,
                          pex = pex,
                          bex = bex,
                          pmid = pmid,
                          threshold = threshold)

              # Visualize
              # plot(hist(x, breaks = seq(-0.05, 1.05, by = 0.01)))

              # Calculate empirical proportions
              prop_zeros <- sum(abs(x) < 1e-10) / n_samples
              prop_ones <- sum(abs(x - 1) < 1e-10) / n_samples
              prop_threshold <- sum(abs(x - threshold) < 1e-10) / n_samples

              # Expected proportions
              expected_prop_threshold <- pmid
              expected_prop_extremes <- (1 - pmid) * pex
              expected_prop_zeros <- expected_prop_extremes * (1 - bex)
              expected_prop_ones <- expected_prop_extremes * bex

              # Scale by p
              expected_prop_zeros <- 2 * (expected_prop_zeros * (1 - p))
              expected_prop_ones <- 2 * (expected_prop_ones * p)

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
                          tolerance = tolerance,
                          label = paste0("Total proportion of extremes with ", label))
            }
          }
        }
      }
    }
  }
})



test_that("dchoco matches rbext empirical distribution and integrates correctly", {
  # Set seed for reproducibility
  set.seed(456)

  # Test parameters
  n_samples <- 15000 # Increased samples for better empirical estimates
  threshold <- 0.5
  tol_prob <- 0.025 # Tolerance for comparing probabilities (empirical vs theoretical)
  tol_integration <- 0.01 # Tolerance for numerical integration check
  tol_ks_p <- 0.001 # Minimum p-value for KS test to pass

  # Reduced parameter set for faster testing, focusing on variety
  p_values <- c(0.2, 0.5, 0.8)
  pex_values <- c(0, 0.2, 0.6)
  muleft_values <- c(0.3, 0.7)
  mudelta_values <- c(-0.5, 0.5)
  bex_values <- c(0.2, 0.8)
  pmid_values <- c(0, 0.3)
  phileft_values <- c(2, 5) # Test different precisions
  phidelta_values <- c(-0.5, 0.5) # Test different precision changes

  # Loop through parameter combinations
  for (p in p_values) {
    for (pex in pex_values) {
      for (muleft in muleft_values) {
        for (mudelta in mudelta_values) {
          for (bex in bex_values) {
            for (pmid in pmid_values) {
              for (phileft in phileft_values) {
                for (phidelta in phidelta_values) {

                  # Generate descriptive label
                  label <- sprintf("p=%.1f, muL=%.1f, muD=%.1f, phiL=%.1f, phiD=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                                   p, muleft, mudelta, phileft, phidelta, pex, bex, pmid)

                  # Generate data
                  x_sample <- rchoco(n = n_samples, p = p, muleft = muleft, mudelta = mudelta,
                                      phileft = phileft, phidelta = phidelta, pex = pex,
                                      bex = bex, pmid = pmid, threshold = threshold)

                  # Calculate empirical proportions
                  emp_p0 <- mean(abs(x_sample) < 1e-10)
                  emp_p1 <- mean(abs(x_sample - 1) < 1e-10)
                  emp_pth <- mean(abs(x_sample - threshold) < 1e-10)
                  emp_pleft_cont <- mean(x_sample > 1e-10 & x_sample < threshold - 1e-10)
                  emp_pright_cont <- mean(x_sample > threshold + 1e-10 & x_sample < 1 - 1e-10)

                  # Calculate theoretical probabilities at point masses
                  theo_p0 <- dchoco(0, p = p, muleft = muleft, mudelta = mudelta, phileft = phileft,
                                     phidelta = phidelta, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  theo_p1 <- dchoco(1, p = p, muleft = muleft, mudelta = mudelta, phileft = phileft,
                                     phidelta = phidelta, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  theo_pth <- dchoco(threshold, p = p, muleft = muleft, mudelta = mudelta, phileft = phileft,
                                      phidelta = phidelta, pex = pex, bex = bex, pmid = pmid, threshold = threshold)

                  # --- Check Point Masses ---
                  expect_equal(emp_p0, theo_p0, tolerance = tol_prob, label = paste(label, "- P(0)"))
                  expect_equal(emp_p1, theo_p1, tolerance = tol_prob, label = paste(label, "- P(1)"))
                  expect_equal(emp_pth, theo_pth, tolerance = tol_prob, label = paste(label, "- P(threshold)"))

                  # --- Integration Check ---
                  integrand <- function(x_int) {
                    dchoco(x_int, p = p, muleft = muleft, mudelta = mudelta, phileft = phileft,
                            phidelta = phidelta, pex = pex, bex = bex, pmid = pmid, threshold = threshold)
                  }

                  # Integrate continuous parts (handle potential errors)
                  integral_left <- tryCatch(
                    stats::integrate(integrand, lower = .Machine$double.eps, upper = threshold - .Machine$double.eps,
                                     subdivisions = 200)$value,
                    error = function(e) 0 # Assume 0 if integration fails (e.g., no density there)
                  )
                  integral_right <- tryCatch(
                    stats::integrate(integrand, lower = threshold + .Machine$double.eps, upper = 1 - .Machine$double.eps,
                                     subdivisions = 200)$value,
                    error = function(e) 0
                  )

                  # Check total probability sums to 1
                  total_prob <- theo_p0 + theo_p1 + theo_pth + integral_left + integral_right
                  expect_equal(total_prob, 1, tolerance = tol_integration, label = paste(label, "- Total Probability Integration"))

                  # --- KS Test for Continuous Parts ---
                  # Calculate parameters needed for underlying Beta distributions
                  logit_muleft <- log(muleft / (1 - muleft))
                  logit_muright <- logit_muleft + mudelta
                  muright <- exp(logit_muright) / (1 + exp(logit_muright))
                  phiright <- phileft * exp(phidelta)
                  eps_mu <- 1e-9
                  muright <- pmax(eps_mu, pmin(muright, 1 - eps_mu))
                  muleft_clamped <- pmax(eps_mu, pmin(muleft, 1 - eps_mu))
                  pex_left_eff <- pmin(1, pmax(0, (1 - bex) * (pex * 2)))
                  pex_right_eff <- pmin(1, pmax(0, bex * (pex * 2)))
                  shape1_left <- muleft_clamped * phileft * 2
                  shape2_left <- (1 - muleft_clamped) * phileft * 2
                  shape1_right <- muright * phiright * 2
                  shape2_right <- (1 - muright) * phiright * 2

                  # KS Test for Left Continuous Part
                  if (emp_pleft_cont > 0.01 && pex_left_eff < 1) { # Need enough data and non-zero beta density
                    x_left_cont <- x_sample[x_sample > 1e-10 & x_sample < threshold - 1e-10]
                    if (length(x_left_cont) > 10) {
                      y_raw_left <- 1 - x_left_cont / threshold # Transform back to [0, 1]
                      y_raw_left_jittered <- y_raw_left + runif(length(y_raw_left), -1e-9, 1e-9) # Add jitter
                      y_raw_left_jittered <- pmax(1e-12, pmin(1 - 1e-12, y_raw_left_jittered)) # Clamp
                      ks_left <- stats::ks.test(y_raw_left_jittered, "pbeta", shape1 = shape1_left, shape2 = shape2_left)
                      expect_gt(ks_left$p.value, tol_ks_p,
                                label = paste(label, "- KS Test Left Continuous (p=", format(ks_left$p.value, digits=3), ")"))
                    }
                  }

                  # KS Test for Right Continuous Part
                  if (emp_pright_cont > 0.01 && pex_right_eff < 1) { # Need enough data and non-zero beta density
                    x_right_cont <- x_sample[x_sample > threshold + 1e-10 & x_sample < 1 - 1e-10]
                     if (length(x_right_cont) > 10) {
                      y_raw_right <- (x_right_cont - threshold) / (1 - threshold) # Transform back to [0, 1]
                      y_raw_right_jittered <- y_raw_right + runif(length(y_raw_right), -1e-9, 1e-9) # Add jitter
                      y_raw_right_jittered <- pmax(1e-12, pmin(1 - 1e-12, y_raw_right_jittered)) # Clamp
                      ks_right <- stats::ks.test(y_raw_right_jittered, "pbeta", shape1 = shape1_right, shape2 = shape2_right)
                      expect_gt(ks_right$p.value, tol_ks_p,
                                label = paste(label, "- KS Test Right Continuous (p=", format(ks_right$p.value, digits=3), ")"))
                    }
                  }
                } # phidelta
              } # phileft
            } # pmid
          } # bex
        } # mudelta
      } # muleft
    } # pex
  } # p
})



context("CHOCO - brms")

test_that("choco model can recover parameters with brms using variational inference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # --- 1. Set up Simulation ---
  true_mu <- 0.6 # Renamed from true_p for clarity in brms context
  true_muleft <- 0.3
  true_mudelta <- 0.5
  true_phileft <- 5
  true_phidelta <- -0.2
  true_pex <- 0.1
  true_bex <- 0.7
  n_obs <- 3000

  # Generate synthetic data
  set.seed(456)
  df <- data.frame(
    y = rchoco(n = n_obs, p = true_mu, muleft = true_muleft, mudelta = true_mudelta, # Pass true_mu to 'p' arg of rchoco
                phileft = true_phileft, phidelta = true_phidelta, pex = true_pex,
                bex = true_bex, pmid = 0, threshold = 0.5)
  )

  # --- 2. Define and Fit brms Model ---
  # Formula for intercept-only model
  # 'y ~ 1' models the intercept for the 'mu' parameter by default
  f <- brms::bf(
    y ~ 1,
    muleft ~ 1,
    mudelta ~ 1,
    phileft ~ 1,
    phidelta ~ 1,
    pex ~ 1,
    bex ~ 1,
    pmid = 0,
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
  post_summary <- brms::posterior_summary(m, probs = c(0.05, 0.95)) # Keep CI for reference if needed

  # Apply inverse-link functions to recover parameters
  post_mu <- brms::inv_logit_scaled(post_summary["b_Intercept", "Estimate"])
  post_muleft <- brms::inv_logit_scaled(post_summary["b_muleft_Intercept", "Estimate"])
  post_mudelta <- post_summary["b_mudelta_Intercept", "Estimate"]
  post_phileft <- log(1 + exp(post_summary["b_phileft_Intercept", "Estimate"])) # softplus
  post_phidelta <- post_summary["b_phidelta_Intercept", "Estimate"]
  post_pex <- brms::inv_logit_scaled(post_summary["b_pex_Intercept", "Estimate"])
  post_bex <- brms::inv_logit_scaled(post_summary["b_bex_Intercept", "Estimate"])

  expect_equal(post_mu,  true_mu, tolerance = 0.1,
               label = sprintf("Posterior mean of mu (%.3f) is close to true mu (%.3f)", post_mu, true_mu))
  expect_equal(post_muleft, true_muleft, tolerance = 0.1,
               label = sprintf("Posterior mean of muleft (%.3f) is close to true muleft (%.3f)", post_muleft, true_muleft))
  expect_equal(post_mudelta, true_mudelta, tolerance = 0.1,
               label = sprintf("Posterior mean of mudelta (%.3f) is close to true mudelta (%.3f)", post_mudelta, true_mudelta))
  expect_equal(post_phileft, true_phileft, tolerance = 1,
               label = sprintf("Posterior mean of phileft (%.3f) is close to true phileft (%.3f)", post_phileft, true_phileft))
  expect_equal(post_phidelta, true_phidelta, tolerance = 0.3,
               label = sprintf("Posterior mean of phidelta (%.3f) is close to true phidelta (%.3f)", post_phidelta, true_phidelta))
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
  # This assumes the Stan function hardcodes threshold = 0.5
  # If the Stan function takes threshold as an argument, adjust the call below
  choco_lpdf_stan <- choco_lpdf_expose()

  # --- Define parameter grids ---
  y_values <- c(0, 0.01, 0.49, 0.5, 0.51, 0.99, 1) # Test around threshold
  p_values <- c(0.1, 0.5, 0.9)
  muleft_values <- c(0.2, 0.7)
  mudelta_values <- c(-0.5, 0.5)
  phileft_values <- c(2, 6)
  phidelta_values <- c(-0.3, 0.3)
  pex_values <- c(0, 0.2, 0.8, 1)
  bex_values <- c(0, 0.5, 1)
  pmid_values <- c(0, 0.1, 0.6, 1)

  # --- Loop through parameter combinations ---
  for (p in p_values) {
    for (muleft in muleft_values) {
      for (mudelta in mudelta_values) {
        for (phileft in phileft_values) {
          for (phidelta in phidelta_values) {
            for (pex in pex_values) {
              for (bex in bex_values) {
                for (pmid in pmid_values) {
                  # Test across different y values
                  for (y in y_values) {
                    label <- sprintf(
                      "y=%.2f, p=%.1f, muleft=%.1f, mudelta=%.1f, phileft=%.1f, phidelta=%.1f, pex=%.1f, bex=%.1f, pmid=%.1f",
                      y, p, muleft, mudelta, phileft, phidelta, pex, bex, pmid
                    )

                    # Calculate log-density using Stan function
                    # Note: The exposed Stan function here doesn't take threshold
                    stan_log_lik <- choco_lpdf_stan(y, p, muleft, mudelta, phileft, phidelta, pex, bex, pmid)

                    # Calculate log-density using R function
                    r_log_lik <- dchoco(y, p, muleft, mudelta, phileft, phidelta, pex, bex, pmid,
                      threshold = 0.5, log = TRUE
                    )

                    # Visualize
                    # plot(seq(0, 1, length.out = 100),
                    #      dchoco(seq(0, 1, length.out = 100), p, muleft, mudelta, phileft, phidelta, pex, bex, pmid),
                    #      type = "l", main = label)

                    # Compare log-likelihoods
                    expect_equal(stan_log_lik, r_log_lik,
                      tolerance = 1e-6,
                      label = paste("Log-likelihood comparison:", label)
                    )
                  } # y
                } # pmid
              } # bex
            } # pex
          } # phidelta
        } # phileft
      } # mudelta
    } # muleft
  } # p

  # --- Test invalid parameter handling ---
  valid_params <- list(x=0.2, p=0.5, muleft=0.5, mudelta=0, phileft=4, phidelta=0, pex=0.1, bex=0.5, pmid=0.1)

  # R function errors
  expect_warning(do.call(dchoco, c(valid_params[-2], list(p=-0.1))), "p must be between 0 and 1")
  expect_warning(do.call(dchoco, c(valid_params[-3], list(muleft=1.1))), "muleft must be between 0 and 1")
  expect_warning(do.call(dchoco, c(valid_params[-5], list(phileft=-1))), "phileft must be positive")

  # Stan function returns -Inf
  expect_equal(choco_lpdf_stan(y=0.2, mu=-0.1, muleft=0.5, mudelta=0, phileft=4, phidelta=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid p")
  expect_equal(choco_lpdf_stan(y=0.2, mu=0.5, muleft=1.1, mudelta=0, phileft=4, phidelta=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid muleft")
  expect_equal(choco_lpdf_stan(y=0.2, mu=0.5, muleft=0.5, mudelta=0, phileft=-1, phidelta=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid phileft")
  expect_equal(choco_lpdf_stan(y=-0.1, mu=0.5, muleft=0.5, mudelta=0, phileft=4, phidelta=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid y (<0)")
  expect_equal(choco_lpdf_stan(y=1.1, mu=0.5, muleft=0.5, mudelta=0, phileft=4, phidelta=0, pex=0.1, bex=0.5, pmid=0.1), -Inf, label="Stan invalid y (>1)")

})


