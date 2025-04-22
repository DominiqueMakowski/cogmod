context("DDM - rddm and dddm")


test_that("rddm and dddm produces correct output shape", {
  n <- 1000
  drift <- 0.2
  bs <- 1
  bias <- 0.5
  ndt <- 0.3

  # Run rddm
  sim_data <- rddm(n, drift = drift, bs = bs, bias = bias, ndt = ndt)

  # Check that the output is a data frame
  expect_s3_class(sim_data, "data.frame")

  # Check that the data frame has the correct number of rows
  expect_equal(nrow(sim_data), n)

  # Check that the data frame has the correct columns
  expect_named(sim_data, c("rt", "response"))

  # DDM
  x <- c(0.5, 0.7, 1.2)
  drift <- 0.2
  bs <- 1
  bias <- 0.5
  ndt <- 0.3
  resp <- c(0, 1, 0)

  # Run dddm
  densities <- dddm(x = x, drift = drift, bs = bs, bias = bias, ndt = ndt, response = resp)

  # Check that the output is a numeric vector
  expect_type(densities, "double")

  # Check that the output has the correct length
  expect_equal(length(densities), length(x))
})


test_that("dddm matches densities of rddm simulated data", {
  set.seed(123)

  # Parameters for the test
  n <- 1000
  drift <- 0.5
  bs <- 1.2
  bias <- 0.6
  ndt <- 0.3

  # Simulate data using rddm
  sim_data <- rddm(n, drift = drift, bs = bs, bias = bias, ndt = ndt)

  # Extract reaction times and responses
  rt <- sim_data$rt
  response <- sim_data$response

  # Compute densities using dddm
  densities <- dddm(x = rt, drift = drift, bs = bs, bias = bias, ndt = ndt, response = response, log = FALSE)

  # Check that densities are finite and non-negative
  expect_true(all(is.finite(densities)), "All densities should be finite")
  expect_true(all(densities >= 0), "All densities should be non-negative")

  # Define plausible RTs as those within a reasonable central range, excluding potential tails
  # Original upper bound (ndt + bs / abs(drift)) might include all samples for some seeds/params.
  # Use a tighter upper bound to better isolate the central tendency.
  plausible_upper_bound <- ndt + bs / (2 * abs(drift)) # e.g., ndt + half the mean decision time approx.
  plausible_rt_indices <- rt > ndt & rt < plausible_upper_bound
  plausible_densities <- densities[plausible_rt_indices]

  # Check that plausible RTs have higher average density than the overall average density
  expect_gt(mean(plausible_densities), mean(densities),
            "Mean density of plausible RTs should be higher than overall mean density")
})

context("DDM - brms")

test_that("DDM model can recover parameters with brms", {
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  set.seed(1234)
  n_obs <- 3000

  # True parameters
  true_nu <- 0.5
  true_bs <- 1.2
  true_bias <- 0.6
  true_tau <- 0.8 # Proportion of minrt
  true_minrt <- 0.3
  true_ndt <- true_tau * true_minrt
  true_minrt <- true_minrt

  # Simulate data
  df <- rddm(n_obs, drift = true_nu, bs = true_bs, bias = true_bias, ndt = true_ndt)

  # Define brms formula
  f <- brms::bf(
    rt | dec(response) ~ 1,
    bs ~ 1,
    bias ~ 1,
    tau ~ 1,
    minrt = true_minrt
  )

  # Fit model using brms
  fit <- brms::brm(
    formula = f,
    data = df,
    family = ddm(),
    stanvars = ddm_stanvars(),
    algorithm = "pathfinder",  # Use VI
    backend = "cmdstanr",
    refresh = 0,
    init = 0
  )

  # Extract posterior means from summary
  post <- brms::posterior_summary(fit)
  means <- post[, "Estimate"]

  # Check parameter recovery
  expect_equal(means[["b_Intercept"]], true_nu, tolerance = 0.1, label = "drift recovery")
  expect_equal(log1p(exp(means[["b_bs_Intercept"]])), true_bs, tolerance = 0.1, label = "bs recovery")
  expect_equal(brms::inv_logit_scaled(means[["b_bias_Intercept"]]), true_bias, tolerance = 0.1, label = "bias recovery")
  expect_equal(brms::inv_logit_scaled(means[["b_tau_Intercept"]]), true_tau, tolerance = 0.1, label = "tau recovery")

  # Check derived ndt recovery
  est_tau_mean <- brms::inv_logit_scaled(means[["b_tau_Intercept"]])
  est_ndt_mean <- est_tau_mean * min(df$rt)
  expect_equal(est_ndt_mean, true_ndt, tolerance = 0.2, label = "Derived NDT recovery")

  # --- Test Post-processing Functions ---
  # Test posterior prediction dimensions
  n_pred_draws <- 10
  n_pred_obs <- 5
  # Ensure newdata has the 'response' column if needed by posterior_predict_ddm
  # (It's not strictly needed for prediction itself, but good practice)
  newdata_pred <- df[1:n_pred_obs, , drop = FALSE]

  # Generate predictions
  pred <- brms::posterior_predict(fit, ndraws = n_pred_draws, newdata = newdata_pred)

  # Check dimensions: posterior_predict for wiener returns draws x (obs*2) matrix
  # Columns alternate between RT (q) and response (resp)
  expect_true(is.matrix(pred), label = "posterior_predict output should be a matrix")
  expect_equal(nrow(pred), n_pred_draws, label = "posterior_predict rows should match ndraws")
  expect_equal(ncol(pred), n_pred_obs * 2, label = "posterior_predict columns should be 2 * nrow(newdata)")

  # Optional: Check predicted RTs are plausible (e.g., > estimated NDT)
  pred_rt_indices <- seq(1, ncol(pred), by = 2) # Indices for RT columns
  pred_rt <- pred[, pred_rt_indices]
  expect_true(all(pred_rt > est_ndt_mean * 0.9, na.rm = TRUE), # Allow some buffer
              label = "Predicted RTs should generally be > estimated NDT")

  # Optional: Check predicted responses are plausible (0 or 1)
  pred_resp_indices <- seq(2, ncol(pred), by = 2) # Indices for response columns
  pred_resp <- pred[, pred_resp_indices]
  expect_true(all(pred_resp %in% c(0, 1)), label = "Predicted responses should be 0 or 1")

  # Test log-likelihood dimensions
  ll <- brms::log_lik(fit, ndraws = 5) # Use a small number of draws for testing
  expect_true(is.matrix(ll), label = "log_lik output should be a matrix")
  expect_equal(nrow(ll), 5, label = "log_lik rows should match ndraws")
  expect_equal(ncol(ll), n_obs, label = "log_lik columns should match number of original observations")
  expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")
})


test_that("Stan DDM lpdf matches R dddm function", {
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  ddm_lpdf_stan <- ddm_lpdf_expose()

  # Define parameter grids for testing
  Y_values <- c(0.5, 0.7, 1.0, 1.5)
  drift_values <- c(0.2, 0.5)
  bs_values <- c(1.0, 1.5)
  bias_values <- c(0.3, 0.6)
  tau_values <- c(0.5, 0.8)
  minrt_values <- c(0.2, 0.3)
  dec_values <- c(0, 1)

  for (drift in drift_values) {
    for (bs in bs_values) {
      for (bias in bias_values) {
        for (tau in tau_values) {
          for (minrt in minrt_values) {
            ndt <- tau * minrt

            for (dec in dec_values) {
              for (Y in Y_values) {
                if (Y <= ndt) next # Skip invalid cases where Y <= ndt

                # Calculate lpdf using Stan function
                stan_lpdf <- ddm_lpdf_stan(Y, drift, bs, bias, tau, minrt, dec)

                # Calculate lpdf using R function
                r_lpdf <- dddm(x = Y, drift = drift, bs = bs, bias = bias, ndt = ndt, response = dec, log = TRUE)

                # Compare results
                label <- sprintf("Y=%.2f, drift=%.2f, bs=%.2f, bias=%.2f, tau=%.2f, minrt=%.2f, dec=%d",
                                 Y, drift, bs, bias, tau, minrt, dec)
                expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-6, label = label)
              }
            }
          }
        }
      }
    }
  }
})
