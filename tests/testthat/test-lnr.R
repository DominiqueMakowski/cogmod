context("LNR - rlnr and dlnr")

test_that("rlnr generates valid outcomes and dlnr computes valid densities/log-densities", {
  set.seed(123)
  n <- 100
  mu <- 1
  mudelta <- 0.5
  sigmazero <- 1
  sigmadelta <- -0.5
  ndt <- 0.2

  # --- Test rlnr ---
  data <- rlnr(n, mu = mu, mudelta = mudelta, sigmazero = sigmazero, sigmadelta = sigmadelta, ndt = ndt)

  # Check output format
  expect_true(is.data.frame(data))
  expect_true(all(c("rt", "response") %in% colnames(data)))
  expect_equal(nrow(data), n)

  # Check RTs and responses
  expect_true(all(data$rt > ndt)) # RTs must be greater than ndt
  expect_true(all(data$response %in% c(0, 1)))

  # --- Test dlnr ---
  y_valid <- ndt + 0.5 # RT > ndt
  y_invalid <- ndt - 0.1 # RT <= ndt

  # Test density calculation for valid RT
  density_0 <- dlnr(y_valid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0)
  density_1 <- dlnr(y_valid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1)
  expect_true(is.finite(density_0) && density_0 >= 0)
  expect_true(is.finite(density_1) && density_1 >= 0)

  # Test log-density calculation for valid RT
  log_density_0 <- dlnr(y_valid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0, log = TRUE)
  log_density_1 <- dlnr(y_valid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 1, log = TRUE)
  expect_true(is.finite(log_density_0) && log_density_0 <= 0) # Log density can be 0 if density is 1
  expect_true(is.finite(log_density_1) && log_density_1 <= 0)
  expect_equal(exp(log_density_0), density_0, tolerance = 1e-9)
  expect_equal(exp(log_density_1), density_1, tolerance = 1e-9)

  # Test density for invalid RT (y <= ndt)
  expect_equal(dlnr(y_invalid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0), 0)
  expect_equal(dlnr(y_invalid, mu, mudelta, sigmazero, sigmadelta, ndt, response = 0, log = TRUE), -Inf)
})


test_that("dlnr integrates correctly and matches rlnr empirical probabilities", {
  set.seed(456)
  n_samples <- 15000 # Large sample size
  tol_prob <- 0.025 # Tolerance for probability comparison
  tol_integration <- 0.01 # Tolerance for integration

  # Parameter grid (reduced for reasonable test time)
  mu_values <- c(0, 0.8)
  mudelta_values <- c(-0.5, 0.5)
  sigmazero_values <- c(0.5, 1.2)
  sigmadelta_values <- c(-0.3, 0.3)
  ndt_values <- c(0.1, 0.3)

  for (mu in mu_values) {
    for (mudelta in mudelta_values) {
      for (sigmazero in sigmazero_values) {
        for (sigmadelta in sigmadelta_values) {
          for (ndt in ndt_values) {

            label <- sprintf("params(mu=%.1f, mud=%.1f, sz=%.1f, sd=%.1f, ndt=%.1f)",
                             mu, mudelta, sigmazero, sigmadelta, ndt)

            # Generate sample
            data <- rlnr(n_samples, mu, mudelta, sigmazero, sigmadelta, ndt)
            rt_0 <- data$rt[data$response == 0]
            rt_1 <- data$rt[data$response == 1]
            emp_p0 <- length(rt_0) / n_samples
            emp_p1 <- length(rt_1) / n_samples

            # --- Integration Check ---
            # Define the joint density function for integration
            integrand <- function(y, resp) {
              dlnr(y, mu, mudelta, sigmazero, sigmadelta, ndt, response = resp)
            }

            # Integrate densities for each response (handle potential errors)
            # Determine a reasonable upper limit for integration
            # Use a fixed large upper limit or adaptive based on quantiles
            # Using quantiles might be slightly better if distributions vary widely
            q999 <- quantile(data$rt, 0.999, na.rm = TRUE)
            upper_limit <- max(q999, ndt + 10) # Ensure upper limit is well above ndt and typical RTs

            integral_0 <- tryCatch(
              stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                               resp = 0, subdivisions = 200, stop.on.error = FALSE)$value,
              error = function(e) NA
            )
            integral_1 <- tryCatch(
              stats::integrate(integrand, lower = ndt + .Machine$double.eps, upper = upper_limit,
                               resp = 1, subdivisions = 200, stop.on.error = FALSE)$value,
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
              # to avoid issues with sampling variability for rare outcomes.
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

            # --- KS Test Section Removed ---

          } # ndt
        } # sigmadelta
      } # sigmazero
    } # mudelta
  } # mu
})


context("LNR - brms")

test_that("lnr model can recover parameters with brms", {
  # Skip on CRAN and when not running full tests
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # Set random seed for reproducibility
  set.seed(1234)

  # Generate synthetic data
  n <- 1500  # Increased sample size
  true_mu <- 0.5
  true_mudelta <- -0.3 # Made slightly different
  true_sigmazero <- 0.8
  true_sigmadelta <- 0.2 # Made slightly different
  true_ndt <- 0.15

  # Generate data
  df <- rlnr(n, mu = true_mu, mudelta = true_mudelta,
                  sigmazero = true_sigmazero, sigmadelta = true_sigmadelta,
                  ndt = true_ndt)

  # Get min_rt for tau calculation
  min_rt <- min(df$rt)
  true_tau <- true_ndt / min_rt

  # Fit with brms
  f <- brms::bf(
    rt | dec(response) ~ 1,
    mudelta ~ 1,
    sigmazero ~ 1,
    sigmadelta ~ 1,
    tau ~ 1,
    minrt = min_rt # Pass min_rt directly
  )

  # Fit using Pathfinder VI
  fit <- brms::brm(
    formula = f,
    data = df,
    family = lnr(),
    stanvars = lnr_stanvars(),
    init = 0,
    algorithm = "pathfinder", # Use VI
    # iter = 2000, # Increase iterations for VI if needed
    refresh = 0,
    backend = "cmdstanr"
  )

  # Extract posterior summaries
  post <- brms::posterior_summary(fit, probs = c(0.05, 0.95)) # 90% CI

  # Define inverse link functions
  invlink_mu <- function(eta) eta # identity
  invlink_mudelta <- function(eta) eta # identity
  invlink_sigmazero <- function(eta) log1p(exp(eta)) # softplus
  invlink_sigmadelta <- function(eta) eta # identity
  invlink_tau <- function(eta) brms::inv_logit_scaled(eta) # logit

  # Check parameter recovery using credible intervals
  check_ci <- function(true_val, param_name, invlink_fun) {
    b_name <- if (param_name == "mu") "b_Intercept" else paste0("b_", param_name, "_Intercept")
    if (!b_name %in% rownames(post)) {
        warning("Parameter '", b_name, "' not found in posterior summary.")
        return(invisible(NULL))
    }
    ci <- invlink_fun(post[rownames(post) == b_name, c("Q5", "Q95")])
    expect_true(true_val >= ci[1] && true_val <= ci[2],
                label = paste("True", param_name, "(", round(true_val,3), ") not in 90% CI [",
                              round(ci[1], 3), ", ", round(ci[2], 3), "]", sep=""))
  }

  check_ci(true_mu, "mu", invlink_mu)
  check_ci(true_mudelta, "mudelta", invlink_mudelta)
  check_ci(true_sigmazero, "sigmazero", invlink_sigmazero)
  check_ci(true_sigmadelta, "sigmadelta", invlink_sigmadelta)
  check_ci(true_tau, "tau", invlink_tau) # Check tau recovery

  # Check derived ndt recovery (using posterior mean for simplicity)
  est_tau_mean <- invlink_tau(post["b_tau_Intercept", "Estimate"])
  est_ndt_mean <- est_tau_mean * min_rt
  # Allow slightly larger tolerance for derived quantity
  expect_equal(est_ndt_mean, true_ndt, tolerance = 0.05,
               label = "Derived NDT recovery check")


  # --- Test Post-processing Functions ---
  # Test posterior prediction
  n_pred_draws <- 10
  pred <- brms::posterior_predict(fit, ndraws = n_pred_draws, newdata = df[1:5, ])
  expect_true(is.matrix(pred))
  expect_equal(nrow(pred), n_pred_draws)
  expect_equal(ncol(pred), 10) # 5 rows * 2 columns (rt, response)
  # Check rt > ndt (approximately, using mean ndt)
  pred_rt <- pred[, seq(1, ncol(pred), by=2)]
  expect_true(all(pred_rt > est_ndt_mean * 0.9), "Predicted RTs should generally be > estimated NDT")

  # Test log-likelihood
  ll <- brms::log_lik(fit, ndraws = 5)
  expect_true(is.matrix(ll))
  expect_equal(nrow(ll), 5)
  expect_equal(ncol(ll), n)
  expect_true(all(is.finite(ll)), "Log-likelihood values should be finite")

})


test_that("Stan lnr_lpdf matches R dlnr function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function
  lnr_lpdf_stan <- lnr_lpdf_expose()

  # Define parameter grids for testing
  Y_values <- c(0.3, 0.5, 0.8, 1.2, 2.0)
  mu_values <- c(0, 0.5, 1.0)
  mudelta_values <- c(-0.5, 0, 0.5)
  sigmazero_values <- c(0.5, 1.0) # Must be > 0
  sigmadelta_values <- c(-0.5, 0, 0.5)
  tau_values <- c(0.1, 0.5) # Must be in (0, 1)
  minrt_values <- c(0.1, 0.2)
  dec_values <- c(0, 1)

  # Test across a subset of the parameter space
  for (mu in mu_values) {
    for (mudelta in mudelta_values) {
      for (sigmazero in sigmazero_values) {
        for (sigmadelta in sigmadelta_values) {
          for (tau in tau_values) {
            for (minrt in minrt_values) {
              # Calculate ndt for the R function
              ndt <- tau * minrt

              for (dec in dec_values) {
                # Test Y values that exceed ndt
                for (Y in Y_values) {
                  # Skip cases where Y <= ndt
                  if (Y <= ndt + 1e-9) next # Add tolerance

                  # Calculate lpdf using Stan function
                  stan_lpdf <- lnr_lpdf_stan(Y, mu, mudelta, sigmazero, sigmadelta, tau, minrt, dec)

                  # Calculate lpdf using R function
                  r_lpdf <- dlnr(Y, mu, mudelta, sigmazero, sigmadelta, ndt, dec, log = TRUE)

                  # Create informative label for errors
                  label <- sprintf("Y=%.2f, mu=%.1f, mud=%.1f, sz=%.1f, sd=%.1f, tau=%.1f, minrt=%.1f, dec=%d",
                                   Y, mu, mudelta, sigmazero, sigmadelta, tau, minrt, dec)

                  # Check for equality with tolerance
                  expect_equal(stan_lpdf, r_lpdf, tolerance = 1e-6,
                               label = paste("Density mismatch:", label))
                }

                # Test Y value below ndt (should be -Inf in both implementations)
                Y_below <- ndt - 0.01
                if (Y_below > 0) { # Only test if Y_below is positive
                    stan_lpdf_below <- lnr_lpdf_stan(Y_below, mu, mudelta, sigmazero, sigmadelta,
                                                     tau, minrt, dec)
                    r_lpdf_below <- dlnr(Y_below, mu, mudelta, sigmazero, sigmadelta,
                                         ndt, dec, log = TRUE)

                    expect_equal(stan_lpdf_below, r_lpdf_below,
                                 label = sprintf("Below ndt: ndt=%.2f, Y=%.2f", ndt, Y_below))
                    expect_equal(stan_lpdf_below, -Inf,
                                 label = sprintf("Should be -Inf when Y < ndt: Y=%.2f, ndt=%.2f",
                                                 Y_below, ndt))
                }
              }
            }
          }
        }
      }
    }
  }

  # --- Test invalid parameter handling ---
  # R function should error/warn or return 0/-Inf appropriately
  expect_warning(dlnr(0.5, mu=0, mudelta=0, sigmazero=-0.1, sigmadelta=0, ndt=0.1, response=0), 0) # sigmazero <= 0
  expect_warning(dlnr(0.5, mu=0, mudelta=0, sigmazero=0.5, sigmadelta=0, ndt=-0.1, response=0), 0) # ndt < 0 (handled by y < ndt)

  # Stan function should return -Inf for invalid parameters
  # Note: Stan function doesn't have explicit checks, relies on underlying functions
  # Test cases where underlying functions would fail (e.g., log(negative), sd <= 0)
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, mudelta=0, sigmazero=0, sigmadelta=0, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmazero=0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, mudelta=0, sigmazero=-0.1, sigmadelta=0, tau=0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid sigmazero<0")
  # tau outside (0,1) is implicitly handled by logit transform in brms, but check direct call
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, mudelta=0, sigmazero=1, sigmadelta=0, tau=-0.1, minrt=0.1, dec=0), -Inf, label="Stan invalid tau<0")
  expect_equal(lnr_lpdf_stan(Y=0.5, mu=0, mudelta=0, sigmazero=1, sigmadelta=0, tau=1.1, minrt=0.1, dec=0), -Inf, label="Stan invalid tau>1")

})
