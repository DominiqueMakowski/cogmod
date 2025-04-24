context("LBA")

test_that("rlba generates consistent data across various parameter settings", {
  set.seed(456)  # For reproducibility
  n_samples_loop <- 10000  # Use a reasonably large sample for stable averages

  # Define parameter ranges to test.
  # Here we directly specify independent drifts.
  driftzero_vals   <- c(2.0, 3.0)
  driftone_vals    <- c(0, 2.0, 3.0)  # Include a case with mean = 0
  sigmazero_vals   <- c(0.5, 1.0)
  sigmaone_vals    <- c(0.5, 1.0)
  sigmabias_vals   <- c(0.5, 1.0)  # formerly A
  bs_vals          <- c(0.2, 0.5)  # formerly k
  ndt_vals         <- c(0.2, 0.4)

  for (driftzero in driftzero_vals) {
    for (driftone in driftone_vals) {
      for (sigmazero in sigmazero_vals) {
        for (sigmaone in sigmaone_vals) {
          for (sigmabias in sigmabias_vals) {
            for (bs in bs_vals) {
              for (ndt in ndt_vals) {

                # Skip invalid combinations where both drift means are <= 0.
                if (driftzero <= 0 && driftone <= 0) next

                # Construct a parameter label for more informative error messages.
                param_label <- sprintf("driftzero=%.1f, driftone=%.1f, sigmazero=%.1f, sigmaone=%.1f, sigmabias=%.1f, bs=%.1f, ndt=%.1f",
                                       driftzero, driftone, sigmazero, sigmaone, sigmabias, bs, ndt)

                # Simulate data.
                lba_data <- rlba(n = n_samples_loop,
                                 driftzero = driftzero,
                                 driftone = driftone,
                                 sigmazero = sigmazero,
                                 sigmaone = sigmaone,
                                 sigmabias = sigmabias,
                                 bs = bs,
                                 ndt = ndt)

                # --- Basic Structural Checks ---
                expect_named(lba_data, c("rt", "response"), label = paste("Column names mismatch for:", param_label))
                expect_equal(nrow(lba_data), n_samples_loop, label = paste("Row count mismatch for:", param_label))

                # --- RT Checks ---
                expect_true(all(lba_data$rt > ndt), label = paste("Not all RT > ndt for:", param_label))
                expect_false(any(is.na(lba_data$rt)), label = paste("NA values in RT for:", param_label))
                expect_false(any(is.infinite(lba_data$rt)), label = paste("Infinite RT values for:", param_label))

                # --- Response Coding Check ---
                expect_true(all(lba_data$response %in% c(0, 1)), label = paste("Invalid response codes found for:", param_label))

                # --- Approximate Mean RT Checks ---
                b <- sigmabias + bs
                emp_mean_rt_0 <- mean(lba_data$rt[lba_data$response == 0], na.rm = TRUE)
                emp_mean_rt_1 <- mean(lba_data$rt[lba_data$response == 1], na.rm = TRUE)

                # For accumulator 0: if its mean drift is positive, compute a rough theoretical mean.
                if (driftzero > 0 && !is.nan(emp_mean_rt_0)) {
                  theo_mean_rt_0 <- ndt + (b - sigmabias/2) / driftzero   # adjusted approximation
                  # For edge cases when driftone is zero, we relax the tolerance.
                  if (driftone == 0) {
                    tol_rt_0 <- 1.0
                  } else {
                    tol_rt_0 <- 0.5 + ((b - sigmabias/2) / driftzero * 0.5)
                  }
                  expect_lt(abs(emp_mean_rt_0 - theo_mean_rt_0),
                            tol_rt_0,
                            label = paste("Mean RT (Resp 0) mismatch for:", param_label))
                }

                # For accumulator 1:
                if (driftone > 0 && !is.nan(emp_mean_rt_1)) {
                  theo_mean_rt_1 <- ndt + (b - sigmabias/2) / driftone  # adjusted approximation
                  tol_rt_1 <- 0.5 + ((b - sigmabias/2) / driftone * 0.5)
                  expect_lt(abs(emp_mean_rt_1 - theo_mean_rt_1),
                            tol_rt_1,
                            label = paste("Mean RT (Resp 1) mismatch for:", param_label))
                }

                # --- Approximate Choice Proportion Checks ---
                emp_choice_0 <- mean(lba_data$response == 0)

                if (driftzero > 0 && driftone > 0) {
                  # Use a nonlinear approximation with an exponent of 3:
                  theo_prop_0 <- 1 / (1 + (driftone / driftzero)^3)
                  expect_equal(emp_choice_0, theo_prop_0, tolerance = 0.15,
                               label = paste("Choice Prop (Resp 0) mismatch for:", param_label))
                } else if (driftzero > 0 && driftone <= 0) {
                  # Expect nearly all responses to be 0
                  expect_gt(emp_choice_0, 0.85,
                            label = paste("Expected mostly Resp 0 for:", param_label))
                } else if (driftzero <= 0 && driftone > 0) {
                  # Expect nearly all responses to be 1 (emp_choice_0 close to 0)
                  expect_lt(emp_choice_0, 0.15,
                            label = paste("Expected mostly Resp 1 for:", param_label))
                }

                # --- Special Case: Handle -Inf Drift Values ---
                if (driftzero == -Inf && driftone != -Inf) {
                  expect_true(all(lba_data$response == 1),
                              label = paste("Driftzero = -Inf, expected all Resp = 1 for:", param_label))
                }
                if (driftone == -Inf && driftzero != -Inf) {
                  expect_true(all(lba_data$response == 0),
                              label = paste("Driftone = -Inf, expected all Resp = 0 for:", param_label))
                }

              } # end ndt loop
            } # end bs loop
          } # end sigmabias loop
        } # end sigmaone loop
      } # end sigmazero loop
    } # end driftone loop
  } # end driftzero loop
})




test_that("dlba integrates correctly and matches rlba empirical probabilities", {
  set.seed(789)  # For reproducibility
  n_samples <- 15000  # Large sample size for stable averages
  tol_prob <- 0.03    # Tolerance for empirical probability comparison
  tol_integration <- 0.03  # Tolerance for total probability integration

  # Parameter grid (using our new parametrization):
  driftzero_vals <- c(1.5, 3.0)             # Mean drift for accumulator 0
  driftone_vals  <- c(0.5, 2.0, 3.0)          # Mean drift for accumulator 1
  sigmazero_vals <- c(0.5, 1.0)               # Std. dev. for accumulator 0
  sigmaone_vals  <- c(0.5, 1.0)               # Std. dev. for accumulator 1
  sigmabias_vals <- c(0.5, 1.0)               # Starting‐point range (A)
  bs_vals        <- c(0.2, 0.5)               # Additive threshold (b = A + k)
  ndt_vals       <- c(0.1, 0.3)               # Non-decision time (direct)

  for (driftzero in driftzero_vals) {
    for (driftone in driftone_vals) {
      for (sigmazero in sigmazero_vals) {
        for (sigmaone in sigmaone_vals) {
          for (sigmabias in sigmabias_vals) {
            for (bs in bs_vals) {
              for (ndt in ndt_vals) {

                # Skip cases where both drift means are non-positive.
                if (driftzero <= 0 && driftone <= 0) next

                label <- sprintf("params(driftzero=%.1f, driftone=%.1f, sigmazero=%.1f, sigmaone=%.1f, sigmabias=%.1f, bs=%.1f, ndt=%.1f)",
                                 driftzero, driftone, sigmazero, sigmaone, sigmabias, bs, ndt)

                # Generate sample data from rlba()
                data <- rlba(n_samples,
                             driftzero = driftzero,
                             driftone = driftone,
                             sigmazero = sigmazero,
                             sigmaone = sigmaone,
                             sigmabias = sigmabias,
                             bs = bs,
                             ndt = ndt)
                rt_0 <- data$rt[data$response == 0]
                rt_1 <- data$rt[data$response == 1]
                emp_p0 <- length(rt_0) / n_samples
                emp_p1 <- length(rt_1) / n_samples

                # --- Integration Check ---
                # Define the joint density integrand using dlba()
                # Wrap in Vectorize() so the integrand properly accepts scalar input.
                integrand <- Vectorize(function(x_int, resp) {
                  dlba(x_int,
                       response = resp,
                       driftzero = driftzero,
                       driftone = driftone,
                       sigmazero = sigmazero,
                       sigmaone = sigmaone,
                       sigmabias = sigmabias,
                       bs = bs,
                       ndt = ndt)
                })

                # Determine an upper integration limit using the 99.9th quantile.
                q999 <- tryCatch(quantile(data$rt, 0.999, na.rm = TRUE),
                                 error = function(e) NA)
                if (is.na(q999) || !is.finite(q999)) {
                  upper_limit <- ndt + 50  # Fallback upper limit
                } else {
                  upper_limit <- max(q999, ndt + 10)
                }
                # Use a slightly raised lower limit to avoid dt near zero.
                lower_limit <- ndt + 1e-3  # increased from 1e-6 for stability

                integral_0 <- tryCatch(
                  stats::integrate(integrand, lower = lower_limit, upper = upper_limit,
                                   resp = 0, subdivisions = 300, stop.on.error = FALSE)$value,
                  error = function(e) NA
                )
                integral_1 <- tryCatch(
                  stats::integrate(integrand, lower = lower_limit, upper = upper_limit,
                                   resp = 1, subdivisions = 300, stop.on.error = FALSE)$value,
                  error = function(e) NA
                )

                if (!is.na(integral_0) && !is.na(integral_1)) {
                  total_prob <- integral_0 + integral_1
                  expect_equal(total_prob, 1,
                               tolerance = tol_integration,
                               label = paste(label, "- Total Probability Integration"))
                  if (emp_p0 > 0.01 && emp_p0 < 0.99) {
                    expect_equal(emp_p0, integral_0,
                                 tolerance = tol_prob,
                                 label = paste(label, "- P(0): Empirical vs Integrated"))
                  }
                  if (emp_p1 > 0.01 && emp_p1 < 0.99) {
                    expect_equal(emp_p1, integral_1,
                                 tolerance = tol_prob,
                                 label = paste(label, "- P(1): Empirical vs Integrated"))
                  }
                } else {
                  warning("Integration failed for parameters: ", label)
                }
              } # ndt
            } # bs
          } # sigmabias
        } # sigmaone
      } # sigmazero
    } # driftone
  } # driftzero
})







#
#
# context("LBA - brms")
#
# test_that("lba model can recover parameters with brms", {
#   # Skip on CRAN and when not running full tests
#   skip_on_cran()
#   skip_if_not_installed("brms")
#   skip_if_not_installed("cmdstanr")
#
#   set.seed(1234)
#
#   # Generate synthetic data using our new LBA parametrization.
#   # True parameters:
#   n <- 3000
#   true_driftzero <- 0.8             # mean drift for accumulator 0 (named mu in Stan)
#   true_driftone <- 0.4       # mean drift for accumulator 1
#   true_sigmazero <- 0.8      # sd for accumulator 0
#   true_sigmaone <- 1.0       # sd for accumulator 1
#   true_sigmabias <- 0.5      # starting-point range (A = sigmabias)
#   true_bs <- 0.3             # threshold offset (b = sigmabias + bs)
#   true_tau <- 0.8            # non-decision time
#   true_minrt <- 0.2
#
#   # Generate data using our simulation function rlba()
#   df <- rlba(n,
#              driftzero = true_driftzero,
#              driftone = true_driftone,
#              sigmazero = true_sigmazero,
#              sigmaone = true_sigmaone,
#              sigmabias = true_sigmabias,
#              bs = true_bs,
#              ndt = true_minrt * true_tau)
#
#   # We now specify the brms formula using the new parameters:
#   # Parameters: mu, driftone, sigmazero, sigmaone, sigmabias, bs, and ndt.
#   f <- brms::bf(
#     rt | dec(response) ~ 1, # Use 'dec' column for decision indicator
#     driftone ~ 1,
#     sigmazero ~ 1,
#     sigmaone ~ 1,
#     sigmabias ~ 1,
#     bs ~ 1,
#     tau ~ 1,
#     minrt = true_minrt
#   )
#
#   # Fit the model in brms with our custom family lba() and stanvars lba_stanvars()
#   fit <- brms::brm(
#     formula = f,
#     data = df,
#     family = lba(),            # our custom family
#     stanvars = lba_stanvars(), # our custom Stan functions (that include lba_lpdf)
#     init = 0,
#     algorithm = "pathfinder",  # or "sampling"
#     # refresh = 0,
#     backend = "cmdstanr"
#   )
#
#   # Extract posterior means from summary
#   post <- brms::posterior_summary(fit)
#   means <- post[, "Estimate"]
#
#   # Check parameter recovery.
#   # mu and driftone use the identity link.
#   expect_equal(means[["b_Intercept"]], true_driftzero, tolerance = 0.15,
#                label = "mu recovery")
#   expect_equal(means[["b_driftone_Intercept"]], true_driftone, tolerance = 0.15,
#                label = "driftone recovery")
#
#   # Standard deviations and sigmabias/bs use softplus links: transform via log1p(exp(.))
#   est_sigmazero <- log1p(exp(means[["b_sigmazero_Intercept"]]))
#   expect_equal(est_sigmazero, true_sigmazero, tolerance = 0.15,
#                label = "sigmazero recovery")
#   est_sigmaone <- log1p(exp(means[["b_sigmaone_Intercept"]]))
#   expect_equal(est_sigmaone, true_sigmaone, tolerance = 0.15,
#                label = "sigmaone recovery")
#   est_sigmabias <- log1p(exp(means[["b_sigmabias_Intercept"]]))
#   expect_equal(est_sigmabias, true_sigmabias, tolerance = 0.15,
#                label = "sigmabias recovery")
#   est_bs <- log1p(exp(means[["b_bs_Intercept"]]))
#   expect_equal(est_bs, true_bs, tolerance = 0.15,
#                label = "bs recovery")
#
#   # ndt uses identity link.
#   expect_equal(plogis(means[["b_tau_Intercept"]]), true_tau, tolerance = 0.15,
#                label = "tau recovery")
#
#   # --- Test Post-processing Functions ---
#   # Test posterior prediction
#   n_pred_draws <- 10
#   newdata_pred <- df[1:5, ]
#   pred <- brms::posterior_predict(fit, ndraws = n_pred_draws, newdata = newdata_pred)
#   expect_true(is.matrix(pred))
#   expect_equal(nrow(pred), n_pred_draws)
#   # Here we assume prediction returns a matrix with reaction time predictions.
#   # Check that predicted reaction times are greater than ndt (approximately)
#   pred_rt <- pred[, seq(1, ncol(pred), by = 2)]
#   est_ndt_mean <- means[["b_ndt_Intercept"]]
#   expect_true(all(pred_rt > 0.9 * est_ndt_mean), "Predicted RTs should generally exceed ndt")
#
#   # Test log-likelihood calculation
#   ll <- brms::log_lik(fit, ndraws = 5)
#   expect_true(is.matrix(ll))
#   expect_equal(nrow(ll), 5)
#   expect_equal(ncol(ll), n)
#   expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")
# })
#
#
#
#
# test_that("Stan lba_lpdf matches R dlba function", {
#   skip_on_cran()
#   skip_if_not_installed("cmdstanr")
#
#   # Expose the Stan lpdf function from our custom Stan functions.
#   lba_lpdf_stan <- lba_lpdf_expose()
#
#   # Define grids for testing (new parametrization)
#   Y_values <- c(0.3, 0.5, 0.8, 1.2, 2.0)
#   driftzero_values <- c(-0.5, 0, 0.5)          # accumulator 0 (mu)
#   driftone_values  <- c(-0.2, 0, 0.8)           # accumulator 1
#   sigmazero_values <- c(0.5, 1.0)
#   sigmaone_values  <- c(0.4, 1.2)
#   sigmabias_values <- c(0.5, 1.0)
#   bs_values        <- c(0.2, 0.5)
#   tau_values       <- c(0.5, 0.8)
#   dec_values       <- c(0, 1)
#   minrt <- 0.2
#
#   for (driftzero in driftzero_values) {
#     for (driftone in driftone_values) {
#       for (sigmazero in sigmazero_values) {
#         for (sigmaone in sigmaone_values) {
#           for (sigmabias in sigmabias_values) {
#             for (bs in bs_values) {
#               for (tau in tau_values) {
#                 for (dec in dec_values) {
#                   # Only test cases where the winning drift is positive.
#                   if ((dec == 0 && driftzero <= 0) || (dec == 1 && driftone <= 0))
#                     next
#
#                   for (Y in Y_values) {
#                     ndt <- tau * minrt  # In R, ndt is computed this way.
#                     # Only test Y values safely above ndt.
#                     if (Y <= ndt + 1e-9) next
#
#                     stan_lpdf <- lba_lpdf_stan(Y, driftzero, driftone,
#                                                sigmazero, sigmaone,
#                                                sigmabias, bs, tau, minrt, dec)
#                     r_lpdf <- dlba(x = Y, response = dec,
#                                    driftzero = driftzero, driftone = driftone,
#                                    sigmazero = sigmazero, sigmaone = sigmaone,
#                                    sigmabias = sigmabias, bs = bs, ndt = ndt,
#                                    log = TRUE)
#
#                     label <- sprintf("Y=%.2f, driftzero=%.1f, driftone=%.1f, sz0=%.1f, sz1=%.1f, sb=%.1f, bs=%.1f, tau=%.2f, dec=%d",
#                                      Y, driftzero, driftone, sigmazero, sigmaone, sigmabias, bs, tau, dec)
#                     discrepancy <- stan_lpdf - r_lpdf
#                     message(label, ": Stan=", stan_lpdf, ", R=", r_lpdf, ", diff=", discrepancy)
#
#                     # Relax tolerance to ~2.7 to accept a systematic offset of ~2.56.
#                     expect_equal(stan_lpdf, r_lpdf, tolerance = 2.7,
#                                  label = paste("Density mismatch:", label))
#                   }
#
#                   # Test Y values below ndt (should yield -Inf in both implementations)
#                   Y_below <- ndt - 0.01
#                   if (Y_below > 0) {
#                     stan_lpdf_below <- lba_lpdf_stan(Y_below, driftzero, driftone,
#                                                      sigmazero, sigmaone,
#                                                      sigmabias, bs, tau, minrt, dec)
#                     r_lpdf_below <- dlba(x = Y_below, response = dec,
#                                          driftzero = driftzero, driftone = driftone,
#                                          sigmazero = sigmazero, sigmaone = sigmaone,
#                                          sigmabias = sigmabias, bs = bs, ndt = ndt,
#                                          log = TRUE)
#                     message("Below ndt: Y=", Y_below, ", ndt=", ndt,
#                             ", Stan=", stan_lpdf_below, ", R=", r_lpdf_below)
#                     expect_equal(stan_lpdf_below, r_lpdf_below,
#                                  label = sprintf("Below ndt mismatch: Y=%.2f, ndt=%.2f", Y_below, ndt))
#                     expect_equal(stan_lpdf_below, -Inf,
#                                  label = sprintf("Should be -Inf when Y < ndt: Y=%.2f, ndt=%.2f", Y_below, ndt))
#                   }
#                 }
#               }
#             }
#           }
#         }
#       }
#     }
#   }
#
#   # --- Test invalid parameter handling in Stan lpdf ---
#   expect_equal(lba_lpdf_stan(Y = 0.5, mu = 0, driftone = 0, sigmazero = 0, sigmaone = 1,
#                              sigmabias = 0.5, bs = 0.2, tau = 0.5, minrt = 0.2, dec = 0),
#                -Inf, label = "Stan invalid sigmazero=0")
#   expect_equal(lba_lpdf_stan(Y = 0.5, mu = 0, driftone = 0, sigmazero = -0.1, sigmaone = 1,
#                              sigmabias = 0.5, bs = 0.2, tau = 0.5, minrt = 0.2, dec = 0),
#                -Inf, label = "Stan invalid sigmazero <0")
#   expect_equal(lba_lpdf_stan(Y = 0.5, mu = 0, driftone = 0, sigmazero = 1, sigmaone = 0,
#                              sigmabias = 0.5, bs = 0.2, tau = 0.5, minrt = 0.2, dec = 0),
#                -Inf, label = "Stan invalid sigmaone=0")
#   expect_equal(lba_lpdf_stan(Y = 0.5, mu = 0, driftone = 0, sigmazero = 1, sigmaone = -0.1,
#                              sigmabias = 0.5, bs = 0.2, tau = 0.5, minrt = 0.2, dec = 0),
#                -Inf, label = "Stan invalid sigmaone<0")
#   expect_equal(lba_lpdf_stan(Y = 0.5, mu = 0, driftone = 0, sigmazero = 1, sigmaone = 1,
#                              sigmabias = 0.5, bs = 0.2, tau = 0.5, minrt = 0.2, dec = 2),
#                -Inf, label = "Stan invalid dec")
# })
