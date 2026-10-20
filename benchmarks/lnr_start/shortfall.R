# Step 1 of the LNR start check (README.md): how much log-likelihood the
# whole cold start leaves on the table, as a function of where `nuone` starts.
#
# Every other parameter sits at its 0.3.5 start (mu = 0.7, both sigma 0.5,
# sigmabias 0 or 0.5, ndt at half the data's first percentile, poutlier 0.02)
# and `nuone` is moved through 0.7 (the current start, a 50/50 race), 0.35,
# 0, -0.35 and -0.7. The shortfall is the log-likelihood at the generating
# values minus the log-likelihood at the start, over 2,000 trials, so that
# the figures are comparable with NEWS 0.3.4's for cogmod_lba2().
#
# Two generating sigma pairs: equal (0.5, 0.5), the start's own, which
# isolates `nuone`; and wide (0.6, 1.2), the ratio the decision_making
# article's fit and the Illusion Game fits both show (sigmaone about twice
# sigmazero). Two model variants: sigmabias fixed at 0 (generated at 0), and
# sigmabias free, started at 0.5 and generated at 0.2 - the Illusion Game
# posteriors sit 1.2-3 link units below the 0.5 start, which is 0.03-0.18.
#
# Part B does the same on the decision_making article's data (4,620 trials,
# 7.25% errors), with the fitted model's posterior medians as the reference.
#
#   Rscript benchmarks/lnr_start/shortfall.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))
here <- "benchmarks/lnr_start"
dir.create(file.path(here, "results"), showWarnings = FALSE)
softplus <- function(x) log1p(exp(x))

start <- list(mu = 0.7, sigmazero = 0.5, sigmaone = 0.5, poutlier = 0.02)
cands <- c(0.7, 0.35, 0, -0.35, -0.7)
targets <- c(0.03, 0.12, 0.26, 0.43)
gens <- list(
  equal = list(nuzero = 1.0, sigmazero = 0.5, sigmaone = 0.5),
  wide = list(nuzero = 1.0, sigmazero = 0.6, sigmaone = 1.2)
)
variants <- list(fixed = c(gen = 0, start = 0), free = c(gen = 0.2, start = 0.5))
ndt_true <- 0.25
pout_true <- 0.02
n <- 2000
seeds <- 1:5

# P(accumulator 1 wins). Closed form at sigmabias = 0: log T1 - log T0 is
# Normal(nu0 - nu1, sqrt(s0^2 + s1^2)). Otherwise by simulation.
err_rate <- function(nuone, g, sigmabias) {
  if (sigmabias == 0) {
    return(stats::pnorm((nuone - g$nuzero) / sqrt(g$sigmazero^2 + g$sigmaone^2)))
  }
  set.seed(1)
  sim <- rcogmod_lnr(2e5, nuzero = g$nuzero, nuone = nuone,
                     sigmazero = g$sigmazero, sigmaone = g$sigmaone,
                     ndt = 0, sigmabias = sigmabias, poutlier = 0)
  mean(sim$response == 1)
}
solve_nuone <- function(p, g, sigmabias) {
  stats::uniroot(function(x) err_rate(x, g, sigmabias) - p, c(-6, 6),
                 tol = 1e-4)$root
}
ll <- function(d, nuzero, nuone, sigmazero, sigmaone, ndt, sigmabias,
               poutlier) {
  sum(dcogmod_lnr(d$rt, nuzero = nuzero, nuone = nuone,
                  sigmazero = sigmazero, sigmaone = sigmaone, ndt = ndt,
                  response = d$response, sigmabias = sigmabias,
                  poutlier = poutlier, log = TRUE))
}

# Part A: simulated data -----------------------------------------------------
rows <- list()
for (gn in names(gens)) {
  g <- gens[[gn]]
  for (vn in names(variants)) {
    v <- variants[[vn]]
    for (p in targets) {
      nu1 <- solve_nuone(p, g, v[["gen"]])
      for (s in seeds) {
        set.seed(1000 * s + round(100 * p))
        sim <- rcogmod_lnr(n, nuzero = g$nuzero, nuone = nu1,
                           sigmazero = g$sigmazero, sigmaone = g$sigmaone,
                           ndt = ndt_true, sigmabias = v[["gen"]],
                           poutlier = pout_true)
        d <- data.frame(rt = sim$rt, response = sim$response)
        ll_true <- ll(d, g$nuzero, nu1, g$sigmazero, g$sigmaone, ndt_true,
                      v[["gen"]], pout_true)
        ndt0 <- cogmod:::.ndt_start(d$rt)
        for (cand in cands) {
          ll0 <- ll(d, start$mu, cand, start$sigmazero, start$sigmaone, ndt0,
                    v[["start"]], start$poutlier)
          rows[[length(rows) + 1]] <- data.frame(
            gen = gn, variant = vn, target = p, nuone_true = nu1,
            err_obs = mean(d$response == 1), seed = s, nuone_start = cand,
            ll_true = ll_true, ll_start = ll0, shortfall = ll_true - ll0
          )
        }
      }
    }
  }
}
res <- do.call(rbind, rows)
utils::write.csv(res, file.path(here, "results", "shortfall_sim.csv"),
                 row.names = FALSE)

# Mean over seeds: the shortfall per candidate, and its reduction from the
# current start.
agg <- stats::aggregate(shortfall ~ gen + variant + target + nuone_true + nuone_start,
                        res, mean)
wide <- stats::reshape(agg, idvar = c("gen", "variant", "target", "nuone_true"),
                       timevar = "nuone_start", direction = "wide")
names(wide) <- sub("shortfall.", "nuone=", names(wide), fixed = TRUE)
wide <- wide[order(wide$gen, wide$variant, wide$target), ]
wide$nuone_true <- round(wide$nuone_true, 2)
cols <- grep("^nuone=", names(wide), value = TRUE)
wide[cols] <- round(wide[cols])
wide$gain_0.7_to_0 <- wide[["nuone=0.7"]] - wide[["nuone=0"]]
wide$best <- as.numeric(sub("nuone=", "", cols))[apply(wide[cols], 1, which.min)]
rownames(wide) <- NULL
cat("\nPart A - shortfall of the whole start (log-lik units per 2,000 trials, mean of 5 seeds)\n")
print(wide, row.names = FALSE)
utils::write.csv(wide, file.path(here, "results", "shortfall_sim_summary.csv"),
                 row.names = FALSE)

# Part B: the decision_making article's data ---------------------------------
suppressMessages(library(brms))
m <- readRDS("vignettes/models/m_lnr.rds")
d <- m$data
fe <- brms::fixef(m)[, "Estimate"]
pout <- brms::posterior_summary(m, variable = "poutlier")[, "Estimate"]
sp <- as.numeric(d$Condition == "Speed")
ll_fit <- sum(dcogmod_lnr(
  d$RT, nuzero = fe[["Intercept"]] + sp * fe[["ConditionSpeed"]],
  nuone = fe[["nuone_Intercept"]] + sp * fe[["nuone_ConditionSpeed"]],
  sigmazero = softplus(fe[["sigmazero_Intercept"]]),
  sigmaone = softplus(fe[["sigmaone_Intercept"]]),
  ndt = exp(fe[["ndt_Intercept"]] + sp * fe[["ndt_ConditionSpeed"]]),
  response = d$Error, sigmabias = 0, poutlier = pout, log = TRUE))
ndt0 <- cogmod:::.ndt_start(d$RT)
grid <- c(0.7, 0.35, 0, -0.35, -0.7, -1.0)
b <- do.call(rbind, lapply(grid, function(cand) {
  ll0 <- ll(data.frame(rt = d$RT, response = d$Error), start$mu, cand,
            start$sigmazero, start$sigmaone, ndt0, 0, start$poutlier)
  data.frame(nuone_start = cand,
             implied_err = stats::pnorm((cand - start$mu) / sqrt(2 * 0.5^2)),
             ll_start = ll0, shortfall = ll_fit - ll0)
}))
b$vs_current <- b$shortfall[1] - b$shortfall
b$implied_err <- round(b$implied_err, 3)
b[c("ll_start", "shortfall", "vs_current")] <- round(b[c("ll_start", "shortfall", "vs_current")])
cat(sprintf("\nPart B - decision_making data: %d trials, %.1f%% errors, posterior-median log-lik %.0f\n",
            nrow(d), 100 * mean(d$Error), ll_fit))
print(b, row.names = FALSE)
utils::write.csv(b, file.path(here, "results", "shortfall_vignette.csv"),
                 row.names = FALSE)
