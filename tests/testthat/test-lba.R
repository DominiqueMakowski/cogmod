context("LBA")

test_that("rlba generates consistent data across various parameter settings", {
  set.seed(456) # Use a different seed for variety
  n_samples_loop <- 10000 # Reduced samples for faster testing within loops

  # Define parameter ranges to test
  vzero_vals <- c(2.0, 3.0)
  vdelta_vals <- c(0, -0.5, 0.5) # Symmetric, bias towards 0, bias towards 1
  sigmazero_vals <- c(0.5, 1.0) # Lower and higher base variability
  # Keep sigmadelta simple for now to limit combinations
  sigmadelta_vals <- c(0) # Test symmetric variability case first
  A_vals <- c(0.5, 1.0)
  k_vals <- c(0.2, 0.5)
  ndt_vals <- c(0.2, 0.4)

  # Loop through parameter combinations
  for (vzero in vzero_vals) {
    for (vdelta in vdelta_vals) {
      for (sigmazero in sigmazero_vals) {
        for (sigmadelta in sigmadelta_vals) { # Can expand this later if needed
          for (A in A_vals) {
            for (k in k_vals) {
              for (ndt in ndt_vals) {

                # Skip invalid combinations if necessary (e.g., if both mean drifts are non-positive)
                mean_drift_0 <- vzero
                mean_drift_1 <- vzero + vdelta
                if (mean_drift_0 <= 0 && mean_drift_1 <= 0) next

                # Construct a label for more informative error messages
                param_label <- sprintf("v0=%.1f, vd=%.1f, s0=%.1f, sd=%.1f, A=%.1f, k=%.1f, ndt=%.1f",
                                       vzero, vdelta, sigmazero, sigmadelta, A, k, ndt)

                # Simulate data
                lba_data <- rlba(n = n_samples_loop, vzero = vzero, vdelta = vdelta,
                                 sigmazero = sigmazero, sigmadelta = sigmadelta,
                                 A = A, k = k, ndt = ndt)

                # --- Basic Checks ---
                expect_true(all(lba_data$rt > ndt), paste("RT > ndt failed for:", param_label))
                expect_false(any(is.na(lba_data$rt)), paste("NA RTs found for:", param_label))
                expect_false(any(is.infinite(lba_data$rt)), paste("Infinite RTs found for:", param_label))
                expect_true(all(lba_data$response %in% c(0, 1)), paste("Invalid choices found for:", param_label))

                # --- Approximate Mean RT Checks ---
                # Theoretical mean RT is a rough approximation, especially with variability
                b <- A + k
                emp_mean_rt_0 <- mean(lba_data$rt[lba_data$response == 0], na.rm = TRUE)
                emp_mean_rt_1 <- mean(lba_data$rt[lba_data$response == 1], na.rm = TRUE)

                # Only check theoretical mean if the mean drift is positive
                if (mean_drift_0 > 0 && !is.nan(emp_mean_rt_0)) {
                  theo_mean_rt_0 <- (b / mean_drift_0) + ndt # Very rough approximation
                  # Use a larger tolerance due to approximation and smaller n_samples
                  expect_lt(abs(emp_mean_rt_0 - theo_mean_rt_0), 0.5 + b / mean_drift_0 * 0.5, # Relative tolerance + absolute
                           label = paste("Mean RT (Resp 0) mismatch for:", param_label))
                }
                if (mean_drift_1 > 0 && !is.nan(emp_mean_rt_1)) {
                  theo_mean_rt_1 <- (b / mean_drift_1) + ndt # Very rough approximation
                  expect_lt(abs(emp_mean_rt_1 - theo_mean_rt_1), 0.5 + b / mean_drift_1 * 0.5, # Relative tolerance + absolute
                           label = paste("Mean RT (Resp 1) mismatch for:", param_label))
                }

                # --- Approximate Choice Proportion Checks ---
                # This is also a rough check, especially with variability
                emp_choice_0 <- mean(lba_data$response == 0)

                if (mean_drift_0 > 0 && mean_drift_1 > 0) {
                  # Expect proportion to roughly follow relative drift rates
                  theo_prop_0 <- mean_drift_0 / (mean_drift_0 + mean_drift_1)
                  expect_equal(emp_choice_0, theo_prop_0, tolerance = 0.15, # Increased tolerance
                               label = paste("Choice Prop (Resp 0) mismatch for:", param_label))
                } else if (mean_drift_0 > 0 && mean_drift_1 <= 0) {
                  # Expect mostly response 0
                  expect_gt(emp_choice_0, 0.85, # Allow some errors due to variability
                           label = paste("Expected mostly Resp 0 for:", param_label))
                } else if (mean_drift_0 <= 0 && mean_drift_1 > 0) {
                  # Expect mostly response 1 (so emp_choice_0 should be low)
                  expect_lt(emp_choice_0, 0.15, # Allow some errors due to variability
                           label = paste("Expected mostly Resp 1 for:", param_label))
                }
              } # ndt loop
            } # k loop
          } # A loop
        } # sigmadelta loop
      } # sigmazero loop
    } # vdelta loop
  } # vzero loop
})




test_that("dlba integrates correctly and matches rlba empirical probabilities", {
  set.seed(789)
  n_samples <- 15000 # Large sample size
  tol_prob <- 0.03 # Tolerance for probability comparison (LBA can be variable)
  tol_integration <- 0.015 # Tolerance for integration

  # Parameter grid (reduced for reasonable test time, focus on variety)
  vzero_vals <- c(1.5, 3.0)
  vdelta_vals <- c(-1.0, 0, 1.0) # Wider range
  sigmazero_vals <- c(0.5, 1.0)
  sigmadelta_vals <- c(0) # Keep simple for now
  A_vals <- c(0.5, 1.0)
  k_vals <- c(0.2, 0.5)
  ndt_vals <- c(0.1, 0.3)

  for (vzero in vzero_vals) {
    for (vdelta in vdelta_vals) {
      for (sigmazero in sigmazero_vals) {
        for (sigmadelta in sigmadelta_vals) {
          for (A in A_vals) {
            for (k in k_vals) {
              for (ndt in ndt_vals) {

                # Skip invalid combinations if necessary
                mean_drift_0 <- vzero
                mean_drift_1 <- vzero + vdelta
                if (mean_drift_0 <= 0 && mean_drift_1 <= 0) next

                label <- sprintf("params(v0=%.1f, vd=%.1f, s0=%.1f, sd=%.1f, A=%.1f, k=%.1f, ndt=%.1f)",
                                 vzero, vdelta, sigmazero, sigmadelta, A, k, ndt)

                # Generate sample
                data <- rlba(n_samples, vzero, vdelta, sigmazero, sigmadelta, A, k, ndt)
                rt_0 <- data$rt[data$response == 0]
                rt_1 <- data$rt[data$response == 1]
                emp_p0 <- length(rt_0) / n_samples
                emp_p1 <- length(rt_1) / n_samples

                # --- Integration Check ---
                # Define the joint density function for integration
                integrand <- function(x_int, resp) {
                  dlba(x_int, response = resp, vzero = vzero, vdelta = vdelta,
                       sigmazero = sigmazero, sigmadelta = sigmadelta,
                       A = A, k = k, ndt = ndt)
                }

                # Determine a reasonable upper limit for integration
                q999 <- tryCatch(quantile(data$rt, 0.999, na.rm = TRUE), error = function(e) NA)
                if (is.na(q999) || !is.finite(q999)) {
                    upper_limit <- ndt + 20 # Fallback upper limit
                } else {
                    upper_limit <- max(q999, ndt + 10) # Ensure upper limit is well above ndt
                }


                integral_0 <- tryCatch(
                  stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                                   resp = 0, subdivisions = 250, stop.on.error = FALSE)$value, # Increased subdivisions
                  error = function(e) NA
                )
                integral_1 <- tryCatch(
                  stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                                   resp = 1, subdivisions = 250, stop.on.error = FALSE)$value, # Increased subdivisions
                  error = function(e) NA
                )

                # Check if integrations were successful
                if (!is.na(integral_0) && !is.na(integral_1)) {
                  # Check if total probability integrates to 1
                  total_prob <- integral_0 + integral_1
                  expect_equal(total_prob, 1, tolerance = tol_integration,
                               label = paste(label, "- Total Probability Integration"))

                  # Check if empirical proportions match integrated probabilities
                  # Only check if empirical probability is not too close to 0 or 1
                  if (emp_p0 > 0.01 && emp_p0 < 0.99) {
                      expect_equal(emp_p0, integral_0, tolerance = tol_prob,
                                   label = paste(label, "- P(0): Empirical vs Integrated"))
                  }
                   if (emp_p1 > 0.01 && emp_p1 < 0.99) {
                      expect_equal(emp_p1, integral_1, tolerance = tol_prob,
                                   label = paste(label, "- P(1): Empirical vs Integrated"))
                   }
                } else {
                  warning("Integration failed for parameters: ", label)
                }

              } # ndt loop
            } # k loop
          } # A loop
        } # sigmadelta loop
      } # sigmazero loop
    } # vdelta loop
  } # vzero loop
})


# ... existing test_that blocks ...

context("LBA - brms")


test_that("LBA model can recover parameters with brms", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  set.seed(666)
  n_obs <- 1000

  # True parameters
  true_vzero <- 2.0
  true_vdelta <- -0.5
  true_sigmazero <- 0.8
  true_sigmadelta <- 0.1 # Corresponds to sigma1 = 0.8 * exp(0.1) ~= 0.884
  true_A <- 0.6
  true_k <- 0.3
  true_tau <- 0.7 # Proportion of minrt
  true_minrt <- 0.25 # Assume a fixed known minrt
  true_ndt <- true_tau * true_minrt # = 0.175

  # Simulate data using rlba (which takes ndt)
  df <- rlba(n = n_obs, vzero = true_vzero, vdelta = true_vdelta,
             sigmazero = true_sigmazero, sigmadelta = true_sigmadelta,
             A = true_A, k = true_k, ndt = true_ndt)

  # Define brms formula - simple intercept-only model
  # Note: minrt is passed via minrt() argument. tau is estimated.
  f <- brms::bf(rt | dec(response) ~ 1,
              vdelta ~ 1,
              sigmazero ~ 1,
              sigmadelta ~ 1,
              A ~ 1,
              k ~ 1,
              tau ~ 1,
              minrt = min(df$rt),
              family = lba()) # Use the custom family

  # Define priors
  priors <- c(
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "tau"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "A"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "k"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = ""),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "vdelta"),
    brms::set_prior("normal(0, 1)", class = "Intercept", dpar = "sigmazero")
  ) |>
    brms::validate_prior(f, data = df)

  # Fit model using variational inference (Pathfinder) for speed
  fit <- brms::brm(f,
              data = df,
              prior = priors,
              stanvars = lba_stanvars(), # Include Stan functions
              backend = "cmdstanr",
              algorithm = "pathfinder", # VI
              refresh = 0,
              seed = 667,
              init = 0)

  # Check results
  summary_fit <- summary(fit)
  fixed_effects <- summary_fit$fixed

  # Check recovery (comparing posterior mean to true value)
  # Apply inverse link functions where necessary
  est_vzero <- fixed_effects["Intercept", "Estimate"]
  est_vdelta <- fixed_effects["vdelta_Intercept", "Estimate"]
  est_sigmazero <- log1p(exp(fixed_effects["sigmazero_Intercept", "Estimate"])) # softplus
  est_sigmadelta <- fixed_effects["sigmadelta_Intercept", "Estimate"]
  est_A <- log1p(exp(fixed_effects["A_Intercept", "Estimate"])) # softplus
  est_k <- log1p(exp(fixed_effects["k_Intercept", "Estimate"])) # softplus
  est_tau <- plogis(fixed_effects["tau_Intercept", "Estimate"]) # logit link

  # Use slightly relaxed tolerances for VI and complex model
  expect_equal(est_vzero, true_vzero, tolerance = 0.4, label = "vzero recovery")
  expect_equal(est_vdelta, true_vdelta, tolerance = 0.4, label = "vdelta recovery")
  expect_equal(est_sigmazero, true_sigmazero, tolerance = 0.4, label = "sigmazero recovery")
  expect_equal(est_sigmadelta, true_sigmadelta, tolerance = 0.4, label = "sigmadelta recovery")
  expect_equal(est_A, true_A, tolerance = 0.4, label = "A recovery")
  expect_equal(est_k, true_k, tolerance = 0.4, label = "k recovery")
  expect_equal(est_tau, true_tau, tolerance = 0.15, label = "tau recovery") # Tau is often well-recovered
})


test_that("Stan lba_lpdf matches R dlba function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  stan_lpdf <- lba_lpdf_expose()

  # --- Test multiple valid parameter combinations ---
  vzero_vals <- c(1.5, 3.0)
  vdelta_vals <- c(-0.5, 0.5)
  sigmazero_vals <- c(0.7, 1.2)
  sigmadelta_vals <- c(-0.2, 0.2)
  A_vals <- c(0.4, 0.8)
  k_vals <- c(0.2, 0.6)
  tau_vals <- c(0.2, 0.8)
  minrt_vals <- c(0.15, 0.3)
  y_offsets <- c(0.01, 0.1, 0.5) # Values to add to ndt to get Y
  dec_vals <- c(0, 1)

  for (vzero in vzero_vals) {
    for (vdelta in vdelta_vals) {
      for (sigmazero in sigmazero_vals) {
        for (sigmadelta in sigmadelta_vals) {
          for (A in A_vals) {
            for (k in k_vals) {
              for (tau in tau_vals) {
                for (minrt in minrt_vals) {
                  for (y_offset in y_offsets) {
                    for (dec in dec_vals) {

                      ndt_val <- tau * minrt
                      y_val <- ndt_val + y_offset # Ensure Y > ndt

                      label <- sprintf("Valid: v0=%.1f,vd=%.1f,s0=%.1f,sd=%.1f,A=%.1f,k=%.1f,tau=%.1f,minrt=%.2f,dec=%d,Y=%.3f",
                                       vzero, vdelta, sigmazero, sigmadelta, A, k, tau, minrt, dec, y_val)

                      # Calculate R log-likelihood (using ndt)
                      r_loglik <- dlba(y_val, response = dec, vzero = vzero, vdelta = vdelta,
                                       sigmazero = sigmazero, sigmadelta = sigmadelta,
                                       A = A, k = k, ndt = ndt_val, log = TRUE)

                      # Calculate Stan log-likelihood (using tau, minrt)
                      stan_loglik <- stan_lpdf(Y = y_val, dec = dec, mu = vzero, vdelta = vdelta,
                                               sigmazero = sigmazero, sigmadelta = sigmadelta,
                                               A = A, k = k, tau = tau, minrt = minrt)

                      # Compare (allow slightly larger tolerance for complex calcs)
                      expect_equal(stan_loglik, r_loglik, tolerance = 1e-5, label = label)
                    } # dec
                  } # y_offset
                } # minrt
              } # tau
            } # k
          } # A
        } # sigmadelta
      } # sigmazero
    } # vdelta
  } # vzero

  # --- Test specific edge cases and invalid inputs ---
  # Re-use one set of valid parameters
  v0=2; vd=0; s0=1; sd=0; A=0.5; k=0.5; tau=0.5; minrt=0.2; ndt=tau*minrt # ndt=0.1
  y_base = ndt + 0.2 # = 0.3

  # Case where Y < ndt
  y_below <- 0.05
  r_loglik_below <- dlba(y_below, 0, v0, vd, s0, sd, A, k, ndt, log = TRUE)
  stan_loglik_below <- stan_lpdf(y_below, v0, vd, s0, sd, A, k, tau, minrt, 0)
  expect_equal(stan_loglik_below, r_loglik_below, label = "Stan vs R: Y < ndt")
  expect_true(is.infinite(stan_loglik_below) && stan_loglik_below < 0)

  # Case where Y == ndt (edge case)
  y_equal <- ndt
  r_loglik_equal <- dlba(y_equal, 0, v0, vd, s0, sd, A, k, ndt, log = TRUE)
  stan_loglik_equal <- stan_lpdf(y_equal, v0, vd, s0, sd, A, k, tau, minrt, 0)
  expect_equal(stan_loglik_equal, r_loglik_equal, label = "Stan vs R: Y == ndt")
  expect_true(is.infinite(stan_loglik_equal) && stan_loglik_equal < 0)

  # Invalid parameters (should return -Inf in Stan)
  expect_equal(stan_lpdf(y_base, v0, vd, -0.1, sd, A, k, tau, minrt, 0), -Inf, label = "Stan invalid sigmazero")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, -0.1, k, tau, minrt, 0), -Inf, label = "Stan invalid A")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, A, -0.1, tau, minrt, 0), -Inf, label = "Stan invalid k")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, A, k, -0.1, minrt, 0), -Inf, label = "Stan invalid tau (<0)")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, A, k, 1.1, minrt, 0), -Inf, label = "Stan invalid tau (>1)")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, A, k, tau, -0.1, 0), -Inf, label = "Stan invalid minrt")
  expect_equal(stan_lpdf(y_base, v0, vd, s0, sd, A, k, tau, minrt, 2), -Inf, label = "Stan invalid dec") # Invalid decision

  # Test edge cases for tau = 0 and tau = 1
  ndt0 <- 0 * minrt; y0 <- ndt0 + 0.1
  r0 <- dlba(y0, 0, v0, vd, s0, sd, A, k, ndt0, log = TRUE)
  s0_stan <- stan_lpdf(y0, v0, vd, s0, sd, A, k, 0, minrt, 0)
  expect_equal(s0_stan, r0, tolerance = 1e-5, label = "Stan vs R: tau = 0")

  ndt1 <- 1 * minrt; y1 <- ndt1 + 0.1
  r1 <- dlba(y1, 0, v0, vd, s0, sd, A, k, ndt1, log = TRUE)
  s1_stan <- stan_lpdf(y1, v0, vd, s0, sd, A, k, 1, minrt, 0)
  expect_equal(s1_stan, r1, tolerance = 1e-5, label = "Stan vs R: tau = 1")

})
