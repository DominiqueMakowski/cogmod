# Compile the base and patched programs emitted by emit.R and compare log
# density and gradient at 50 random unconstrained points around the init.
#
#   Rscript benchmarks/stanli/brms_scalar_dpar/compile_check.R <model> [<model> ...]
#
# 2026-10-04, CmdStan 2.38, brms 2.23.1 (base) vs f131ef1 + patch (pr), all
# eleven models of emit.R: lp identical at every finite point or one ULP
# apart (gauss_threads_cens 1/50 at 1.8e-16, gamma_shape_cens at 3.3e-16),
# gradient to 1e-15 relative, finite at the same points.
root <- "benchmarks/results/stanli/brms_scalar_dpar"
# models whose init must respect a data-dependent support (fixed_param
# rejects a bad init before any comparison); the rest of the parameters
# take CmdStan's random inits
inits <- list(
  shifted_lognormal = list(Intercept_ndt = -3),  # ndt = exp(-3) < min(Y) = 0.2
  gamma_shape_cens  = list(Intercept = 1)        # inverse link: mu = inv(1)
)
gp_compile <- function(stan, dir = dirname(stan)) {
  cmdstanr::cmdstan_model(stan, compile_model_methods = TRUE, dir = dir,
                          force_recompile = TRUE, quiet = TRUE)
}
gp_methods <- function(mod, data, init) {
  fit <- mod$sample(data = data, init = init, chains = 1, iter_warmup = 1,
                    iter_sampling = 1, fixed_param = TRUE, refresh = 0, seed = 11,
                    show_messages = FALSE, show_exceptions = FALSE)
  fit$metadata()
  fit$init_model_methods(verbose = FALSE)
  fit
}
for (nm in commandArgs(TRUE)) {
  res <- tryCatch({
    init <- if (is.null(inits[[nm]])) 0 else list(inits[[nm]])
    fits <- lapply(c(base = "base", pr = "pr"), function(w) {
      dir <- file.path(root, w)
      gp_methods(gp_compile(file.path(dir, paste0(nm, ".stan"))), file.path(dir, paste0(nm, ".json")), init)
    })
    up0 <- as.numeric(fits$pr$unconstrain_draws(format = "draws_matrix")[1, ])
    K <- length(up0)
    stopifnot(K == length(fits$base$unconstrain_draws(format = "draws_matrix")[1, ]))
    set.seed(3)
    P <- lapply(1:50, function(k) up0 + rnorm(K, 0, 0.3))
    g <- lapply(fits, function(f) lapply(P, function(q)
      tryCatch(f$grad_log_prob(q), error = function(e) structure(rep(NA_real_, K), log_prob = NA_real_))))
    lp <- sapply(g, function(x) vapply(x, attr, 0, "log_prob"))
    ok <- is.finite(lp[, 1]) & is.finite(lp[, 2])
    dlp <- max(abs(lp[ok, 1] - lp[ok, 2]) / pmax(1, abs(lp[ok, 1])))
    dg <- max(mapply(function(a, b) max(abs(a - b) / pmax(1, abs(a))), g$base[ok], g$pr[ok]))
    sprintf("%-20s compiled both; K = %d; finite at %d/50 (same points: %s); lp max rel diff %.1e (identical at %d); grad max rel diff %.1e",
            nm, K, sum(ok), all(is.finite(lp[, 1]) == is.finite(lp[, 2])), dlp, sum(lp[ok, 1] == lp[ok, 2]), dg)
  }, error = function(e) sprintf("%-20s ERROR: %s", nm, conditionMessage(e)))
  cat(res, "\n")
}
