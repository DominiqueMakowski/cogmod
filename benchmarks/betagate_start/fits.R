# Step 2: does the start reach the sampler? Cold single-chain fits of
# cogmod_betagate() from cogmod_inits() (default jitter), init = 0 and Stan's
# default U(-2, 2), on three of step 1's styles x 4 seeds. Needs cmdstanr; one
# compile, about 15 minutes on a laptop.
#
#   Rscript benchmarks/betagate_start/fits.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
suppressMessages(library(cmdstanr))
work <- "benchmarks/betagate_start/work"
dir.create(work, showWarnings = FALSE, recursive = TRUE)

styles <- list(
  default = c(mu = 0.50, phi = 3,   pex = 0.10, bex = 0.50),
  ushape  = c(mu = 0.50, phi = 0.5, pex = 0.10, bex = 0.50),
  ceiling = c(mu = 0.85, phi = 5,   pex = 0.30, bex = 0.85)
)
f <- brms::bf(y ~ 1, phi ~ 1, pex ~ 1, bex ~ 1, family = cogmod_betagate())
rows <- list()
for (st in names(styles)) {
  s <- styles[[st]]
  set.seed(7)
  d <- data.frame(y = rcogmod_betagate(2000, mu = s[["mu"]], phi = s[["phi"]],
                                       pex = s[["pex"]], bex = s[["bex"]]))
  pr <- cogmod_priors(f, d)
  code <- brms::make_stancode(f, data = d, prior = pr, stanvars = cogmod_stanvars(f),
                              backend = "cmdstanr")
  sdata <- brms::make_standata(f, data = d, prior = pr, stanvars = cogmod_stanvars(f))
  sdata <- lapply(unclass(sdata), identity)
  mod <- cmdstan_model(write_stan_file(code, dir = work), quiet = TRUE)
  arms <- list(cogmod = cogmod_inits(f, d), zero = 0, random = NULL)
  for (arm in names(arms)) for (seed in 1:4) {
    init <- if (arm == "cogmod") list(arms$cogmod(1)) else arms[[arm]]
    fit <- mod$sample(data = sdata, init = init, seed = 100 + seed, chains = 1,
                      iter_warmup = 1000, iter_sampling = 1000, save_warmup = TRUE,
                      refresh = 0, show_messages = FALSE, show_exceptions = FALSE)
    sd <- fit$sampler_diagnostics(inc_warmup = TRUE, format = "draws_df")
    w <- sd$.iteration <= 1000
    dr <- fit$draws(c("Intercept", "Intercept_phi", "Intercept_pex", "Intercept_bex"),
                    format = "draws_df")
    summ <- posterior::summarise_draws(dr, "median", "ess_bulk")
    rows[[length(rows) + 1]] <- data.frame(
      style = st, arm = arm, seed = seed,
      warmup_leapfrog = sum(sd$n_leapfrog__[w]),
      sampling_leapfrog = sum(sd$n_leapfrog__[!w]),
      div = sum(sd$divergent__[!w]), div_warmup = sum(sd$divergent__[w]),
      min_ess = round(min(summ$ess_bulk)),
      mu = round(summ$median[1], 3), phi = round(summ$median[2], 3),
      pex = round(summ$median[3], 3), bex = round(summ$median[4], 3)
    )
  }
}
res <- do.call(rbind, rows)
dir.create("benchmarks/betagate_start/results", showWarnings = FALSE)
utils::write.csv(res, "benchmarks/betagate_start/results/fits.csv", row.names = FALSE)
options(width = 200)
print(aggregate(cbind(warmup_leapfrog, sampling_leapfrog, div, min_ess) ~ style + arm,
                data = res, FUN = stats::median))
