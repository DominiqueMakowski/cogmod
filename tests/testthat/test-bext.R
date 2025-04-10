context("BEXT")

test_that("rbext compared to rordbeta with equivalent parameterization", {
  # Skip if ordbetareg package not available
  skip_if_not_installed("ordbetareg")

  n <- 5000
  set.seed(123)

  mu_values <- c(0.1, 0.5, 0.9)
  phi_values <- c(0.5, 1, 2)
  pex_values <-  c(0, 0.2, 0.5, 0.8, 1)
  bex_values <- c(0, 0.2, 0.5, 0.8, 1)

  for(mu in mu_values) {
    for(phi in phi_values) {
      for(pex in pex_values) {
        for(bex in bex_values) {
          left_prob <- (1 - bex) * pex
          right_prob <- bex * pex

          # In ordbeta, the cutpoints are a function of mu
          cutpoints <- c(qlogis(mu) - qlogis(1 - ((1 - bex) * pex)),
                         qlogis(mu) - qlogis(bex * pex))

          x_bext <- rbext(n, mu = mu, phi = phi, pex = pex, bex = bex)

          # Make sure it generates the right data
          pzero_bext <- sum(x_bext == 0) / n
          pones_bext <- sum(x_bext == 1) / n
          expect_equal(pzero_bext + pones_bext, pex, tolerance = 0.1,
                       label=  "pex != empirical proportion of extremes (bext)" )

          # Compare against ordbeta
          if(all(is.finite(cutpoints)) && diff(cutpoints) < -1e-16) next  # skip if cutpoints are not valid
          x_ordbeta <- ordbetareg::rordbeta(n, mu = mu, phi = phi * 2, cutpoints = cutpoints)

          l_bext <- paste0("rbext(mu=", mu, ", phi=", ", pex=", pex, ", bex=", bex, ")")
          l_ordbeta <- paste0("rordbeta(mu=", mu, ", phi=", phi * 2, ", cutpoints=c(", round(cutpoints[1], 3), ", ", round(cutpoints[2], 3), "))")
          l <- paste("Comparing", l_bext, "vs.", l_ordbeta)

          # Visualize two histograms
          # hist(x_bext, breaks = seq(-0.05, 1.05, by = 0.1))
          # hist(x_ordbeta, breaks = seq(-0.05, 1.05, by = 0.1), col = "red", add = TRUE)

          # Compare the histograms
          h_bext <- hist(x_bext, breaks = seq(-0.05, 1.05, by = 0.1), plot = FALSE)
          h_ordbeta <- hist(x_ordbeta, breaks = seq(-0.05, 1.05, by = 0.1), plot = FALSE)
          expect_equal(mean(abs(h_bext$density - h_ordbeta$density)), 0, tolerance = 0.1,
                       label = paste0(l, ": Histogram comparison"))

          # Compare the proportion of extremes
          pzero_ordbeta <- sum(x_ordbeta == 0) / n
          expect_equal(pzero_bext, pzero_ordbeta, tolerance = 0.05,
                       label = paste0(l, ": Proportion of zeros (bext: ", round(pzero_bext, 2),
                                      ", ordbeta: ", round(pzero_ordbeta, 2), ")"))

          pones_ordbeta <- sum(x_ordbeta == 1) / n
          expect_equal(pones_bext, pones_ordbeta, tolerance = 0.05,
                       label = paste0(l, ": Proportion of ones (bext: ", round(pones_bext, 2),
                                      ", ordbeta: ", round(pones_ordbeta, 2), ")"))

          # Compare mean of middle values
          meanmid_bext <- mean(x_bext[x_bext > 0 & x_bext < 1])
          meanmid_ordbeta <- mean(x_ordbeta[x_ordbeta > 0 & x_ordbeta < 1])
          if(pex < 1) {
            expect_equal(meanmid_bext, meanmid_ordbeta, tolerance = 0.05,
                         label = paste0(l, ": Mean of middle values (bext: ", round(meanmid_bext, 2),
                                        ", ordbeta: ", round(meanmid_ordbeta, 2), ")"))
          }
        }
      }
    }
  }
})


test_that("dbext matches rbext empirical distribution and integrates correctly", {
  n <- 10000 # Use a larger sample size for better empirical estimates
  set.seed(456)

  # Tolerance based roughly on sqrt(p*(1-p)/n), max variance at p=0.5 -> sqrt(0.25/10000) = 0.005
  # Use a slightly larger tolerance to account for multiple sources of variation.
  tol_prob <- 0.05
  tol_integration <- 0.05 # Tolerance for numerical integration check
  tol_ks_p <- 0.00001 # Minimum p-value for KS test to pass

  mu_values <- c(0.2, 0.5, 0.8)
  phi_values <- c(1, 3, 10) # phi=1 -> precision=2 (uniform for mu=0.5)
  pex_values <-  c(0, 0.1, 0.5, 0.9, 1)
  bex_values <- c(0, 0.3, 0.5, 0.7, 1)

  for(mu in mu_values) {
    for(phi in phi_values) {
      for(pex in pex_values) {
        for(bex in bex_values) {

          l <- paste0("params(mu=", mu, ", phi=", phi, ", pex=", pex, ", bex=", bex, ")")

          # Generate sample
          x_sample <- rbext(n, mu = mu, phi = phi, pex = pex, bex = bex)

          # Calculate empirical probabilities
          emp_p0 <- mean(x_sample == 0)
          emp_p1 <- mean(x_sample == 1)
          emp_pmid <- mean(x_sample > 0 & x_sample < 1)

          # Calculate theoretical probabilities/densities at extremes
          theo_p0 <- dbext(0, mu = mu, phi = phi, pex = pex, bex = bex)
          theo_p1 <- dbext(1, mu = mu, phi = phi, pex = pex, bex = bex)
          expected_pmid <- 1 - pex # Theoretical total probability for (0, 1)

          # --- Check Point Masses ---
          expect_equal(emp_p0, theo_p0, tolerance = tol_prob,
                       label = paste(l, "- P(0): Empirical vs Theoretical"))
          expect_equal(emp_p1, theo_p1, tolerance = tol_prob,
                       label = paste(l, "- P(1): Empirical vs Theoretical"))
          expect_equal(emp_pmid, expected_pmid, tolerance = tol_prob,
                       label = paste(l, "- P(0<x<1): Empirical vs Theoretical (1-pex)"))

          # --- Integration Check ---
          # Define the density function for integration (only the continuous part)
          # dbext already includes the (1-pex) scaling for the continuous part
          integrand <- function(x_int) {
            dbext(x_int, mu = mu, phi = phi, pex = pex, bex = bex)
          }

          # Integrate from slightly above 0 to slightly below 1
          # Only integrate if there's probability mass in the middle (pex < 1)
          integral_val <- 0
          if (pex < 1) {
             integral_result <- tryCatch(
                stats::integrate(integrand, lower = .Machine$double.eps, upper = 1 - .Machine$double.eps,
                                 subdivisions = 500, rel.tol = tol_integration)$value,
                error = function(e) {
                    warning("Integration failed for ", l, ": ", e$message)
                    NA # Handle potential integration errors
                }
             )
             if (!is.na(integral_result)) {
                 integral_val <- integral_result
                 # Check if the integral of the continuous part matches 1-pex
                 expect_equal(integral_val, expected_pmid, tolerance = tol_integration,
                              label = paste(l, "- Integral of continuous density matches 1-pex"))
             }
          }

          # Check if total probability (point masses + integral) is 1
          total_prob <- theo_p0 + theo_p1 + integral_val
          expect_equal(total_prob, 1, tolerance = tol_integration,
                       label = paste(l, "- Theoretical density integrates to 1"))


          # --- KS Test for Continuous Part ---
          if (pex < 1 && emp_pmid > 0.01) { # Only if pex < 1 and enough middle values empirically
            x_middle <- x_sample[x_sample > 0 & x_sample < 1]
            if (length(x_middle) > 5) { # KS test needs a reasonable number of points
              # Theoretical parameters for the underlying Beta (unscaled density)
              shape1 <- mu * phi * 2
              shape2 <- (1 - mu) * phi * 2

              # Add jitter to break ties for KS test
              # The amount of jitter should be very small relative to the data range
              jitter_amount <- min(diff(sort(unique(x_middle)))) / 1000 # Heuristic for small jitter
              if (!is.finite(jitter_amount) || jitter_amount == 0) jitter_amount <- 1e-9 # Fallback jitter
              x_middle_jittered <- x_middle + runif(length(x_middle), -jitter_amount, jitter_amount)
              # Ensure jittered values stay within (0, 1) - might be overly cautious but safe
              x_middle_jittered <- pmax(1e-12, pmin(1 - 1e-12, x_middle_jittered))


              # Perform KS test comparing jittered empirical middle data to theoretical Beta CDF
              # The comparison distribution is the *unscaled* Beta
              ks_result <- stats::ks.test(x_middle_jittered, "pbeta", shape1 = shape1, shape2 = shape2)

              # Expect p-value to be reasonably large (distributions are similar)
              expect_gt(ks_result$p.value, tol_ks_p,
                        label = paste0(l, "- KS test p-value for middle distribution (p=",
                                      format(ks_result$p.value, digits=3), ")"))
            }
          }

          # --- Mean Check for Middle Values ---
          if (pex < 1 && emp_pmid > 0.01) {
            x_middle <- x_sample[x_sample > 0 & x_sample < 1]
             if (length(x_middle) > 0) {
                theo_mean_beta <- mu # Mean of the underlying Beta
                emp_mean_middle <- mean(x_middle)
                # Allow larger tolerance for mean comparison, depends heavily on phi
                mean_tolerance <- max(0.1, 1 / sqrt(phi * length(x_middle)/n)) # Heuristic tolerance, adjusted for sample size
                expect_equal(emp_mean_middle, theo_mean_beta, tolerance = mean_tolerance,
                             label = paste(l, "- Mean of middle values"))
             }
          }
        }
      }
    }
  }
})




context("BEXT - brms")

test_that("BEXT model can recover parameters with brms using variational inference", {
  skip_on_cran()
  skip_if_not_installed("brms")
  skip_if_not_installed("cmdstanr")

  # --- 1. Set up Simulation ---
  true_mu <- 0.6
  true_phi <- 4
  true_pex <- 0.15
  true_bex <- 0.4
  n_obs <- 1500 # Increased sample size for better stability

  # Generate synthetic data with known parameters
  set.seed(123)
  df <- data.frame(
    y = rbext(n = n_obs, mu = true_mu, phi = true_phi,
              pex = true_pex, bex = true_bex)
  )

  # --- 2. Define and Fit brms Model ---
  # Formula for intercept-only model (recovering the global parameters)
  f <- brms::bf(
    y ~ 1,
    phi ~ 1,
    pex ~ 1,
    bex ~ 1,
    family = bext() # Use the custom family
  )

  # Fit model using Pathfinder (VI algorithm)
  # Note: VI provides an approximation to the posterior. MCMC might be more accurate but slower.
  # Using default iterations for pathfinder.
  m <- brms::brm(
    formula = f,
    data = df,
    stanvars = bext_stanvars(), # Include Stan functions
    backend = "cmdstanr",
    seed = 123,
    refresh = 0,
    algorithm = "pathfinder" # Use Pathfinder VI
  )

  # --- 3. Check Parameter Recovery ---
  # Extract posterior summary (includes mean, sd, credible intervals)
  post_summary <- brms::posterior_summary(m, probs = c(0.05, 0.95)) # 90% CI

  # Define inverse link functions for clarity
  invlink_mu <- function(eta) brms::inv_logit_scaled(eta)
  invlink_phi <- function(eta) log(1 + exp(eta))  # solftplus link
  invlink_pex <- function(eta) brms::inv_logit_scaled(eta)
  invlink_bex <- function(eta) brms::inv_logit_scaled(eta)

  # Check if true parameters fall within the 90% credible interval
  # Mu
  mu_ci <- invlink_mu(post_summary[rownames(post_summary) == "b_Intercept", c("Q5", "Q95")])
  expect_true(true_mu >= mu_ci[1] && true_mu <= mu_ci[2],
              label = paste("True mu (", true_mu, ") not in 90% CI [", round(mu_ci[1], 3), ", ", round(mu_ci[2], 3), "]", sep=""))

  # Phi
  phi_ci <- invlink_phi(post_summary[rownames(post_summary) == "b_phi_Intercept", c("Q5", "Q95")])
  expect_true(true_phi >= phi_ci[1] && true_phi <= phi_ci[2],
              label = paste("True phi (", true_phi, ") not in 90% CI [", round(phi_ci[1], 3), ", ", round(phi_ci[2], 3), "]", sep=""))

  # Pex
  pex_ci <- invlink_pex(post_summary[rownames(post_summary) == "b_pex_Intercept", c("Q5", "Q95")])
  expect_true(true_pex >= pex_ci[1] && true_pex <= pex_ci[2],
              label = paste("True pex (", true_pex, ") not in 90% CI [", round(pex_ci[1], 3), ", ", round(pex_ci[2], 3), "]", sep=""))

  # Bex
  bex_ci <- invlink_bex(post_summary[rownames(post_summary) == "b_bex_Intercept", c("Q5", "Q95")])
  expect_true(true_bex >= bex_ci[1] && true_bex <= bex_ci[2],
              label = paste("True bex (", true_bex, ") not in 90% CI [", round(bex_ci[1], 3), ", ", round(bex_ci[2], 3), "]", sep=""))


  # --- 4. Test Post-processing Functions ---
  # Test posterior prediction (basic checks)
  n_pred_draws <- 10
  pred <- brms::posterior_predict(m, ndraws = n_pred_draws)
  expect_equal(nrow(pred), n_pred_draws)
  expect_equal(ncol(pred), n_obs)
  expect_true(all(pred >= 0 & pred <= 1), "Posterior predictions outside [0, 1]")
  expect_false(any(is.na(pred)), "NA values in posterior predictions")

  # Test with newdata
  n_newdata <- 5
  pred_new <- brms::posterior_predict(m, ndraws = n_pred_draws, newdata = df[1:n_newdata, ])
  expect_equal(nrow(pred_new), n_pred_draws)
  expect_equal(ncol(pred_new), n_newdata)
  expect_true(all(pred_new >= 0 & pred_new <= 1), "Posterior predictions (newdata) outside [0, 1]")
  expect_false(any(is.na(pred_new)), "NA values in posterior predictions (newdata)")

  # Test log-likelihood calculation (basic checks)
  n_ll_draws <- 5
  ll <- brms::log_lik(m, ndraws = n_ll_draws)
  expect_equal(nrow(ll), n_ll_draws)
  expect_equal(ncol(ll), n_obs)
  expect_true(all(is.finite(ll)), "Non-finite values found in log-likelihood")
})



test_that("Stan bext_lpdf matches R dbext function", {
  skip_on_cran()
  skip_if_not_installed("cmdstanr")

  # Expose the Stan function if possible
  bext_lpdf <- bext_lpdf_expose()

  # --- Define parameter grids for testing ---
  y_values <- c(0, 0.01, 0.5, 0.99, 1) # Test boundaries and mid-points
  mu_values <- c(0.2, 0.5, 0.8)
  phi_values <- c(1, 5, 10) # Include phi=1 (uniform beta case)
  pex_values <- c(0, 0.3, 0.7, 1) # Include boundary cases for pex
  bex_values <- c(0, 0.5, 1)   # Include boundary cases for bex

  # --- Loop through parameter combinations ---
  for (mu in mu_values) {
    for (phi in phi_values) {
      for (pex in pex_values) {
        for (bex in bex_values) {
          # Test across different y values for the current parameter set
          for (y in y_values) {
            label <- sprintf("y=%.2f, mu=%.1f, phi=%.1f, pex=%.1f, bex=%.1f",
                             y, mu, phi, pex, bex)

            # Calculate log-density using Stan function
            stan_log_lik <- bext_lpdf(y, mu, phi, pex, bex)

            # Calculate log-density using R function
            r_log_lik <- dbext(y, mu, phi, pex, bex, log = TRUE)

            # Compare log-likelihoods (handles -Inf correctly)
            expect_equal(stan_log_lik, r_log_lik, tolerance = 1e-6,
                         label = paste("Log-likelihood comparison:", label))
          }
        }
      }
    }
  }

  # --- Test invalid parameter handling ---
  # R function should error
  expect_warning(dbext(0.5, mu = -0.1, phi = 5, pex = 0.1, bex = 0.5), "mu must be between 0 and 1")
  expect_warning(dbext(0.5, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), "phi must be positive")
  expect_warning(dbext(0.5, mu = 0.5, phi = 5, pex = 1.1, bex = 0.5), "pex must be between 0 and 1")
  expect_warning(dbext(0.5, mu = 0.5, phi = 5, pex = 0.1, bex = -0.2), "bex must be between 0 and 1")

  # Stan function should return -Inf for invalid parameters
  expect_equal(bext_lpdf(0.5, mu = -0.1, phi = 5, pex = 0.1, bex = 0.5), -Inf, label = "Stan invalid mu")
  expect_equal(bext_lpdf(0.5, mu = 0.5, phi = -1, pex = 0.1, bex = 0.5), -Inf, label = "Stan invalid phi")
  expect_equal(bext_lpdf(0.5, mu = 0.5, phi = 5, pex = 1.1, bex = 0.5), -Inf, label = "Stan invalid pex")
  expect_equal(bext_lpdf(0.5, mu = 0.5, phi = 5, pex = 0.1, bex = -0.2), -Inf, label = "Stan invalid bex")
  # Also check invalid y
  expect_equal(bext_lpdf(-0.1, mu = 0.5, phi = 5, pex = 0.1, bex = 0.5), -Inf, label = "Stan invalid y (<0)")
  expect_equal(bext_lpdf(1.1, mu = 0.5, phi = 5, pex = 0.1, bex = 0.5), -Inf, label = "Stan invalid y (>1)")
})
