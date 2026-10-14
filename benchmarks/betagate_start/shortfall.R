# Step 1: where should cogmod_betagate() start? Log-likelihood shortfall per
# observation of a start against each simulated data set's own generating
# values, over eight response styles - the measure the choco and betadiscrete
# starts in cogmod_inits.R were chosen on. R only, no Stan.
#
#   Rscript benchmarks/betagate_start/shortfall.R
suppressMessages(pkgload::load_all(".", quiet = TRUE))

styles <- list(
  default   = c(mu = 0.50, phi = 3,   pex = 0.10, bex = 0.50),
  flat      = c(mu = 0.50, phi = 1,   pex = 0.05, bex = 0.50),
  confident = c(mu = 0.50, phi = 8,   pex = 0.05, bex = 0.50),
  high      = c(mu = 0.75, phi = 4,   pex = 0.10, bex = 0.70),
  low       = c(mu = 0.25, phi = 4,   pex = 0.10, bex = 0.30),
  extreme   = c(mu = 0.50, phi = 2,   pex = 0.40, bex = 0.50),
  ushape    = c(mu = 0.50, phi = 0.5, pex = 0.10, bex = 0.50),
  ceiling   = c(mu = 0.85, phi = 5,   pex = 0.30, bex = 0.85)
)

n <- 4000
set.seed(2026)
data <- lapply(styles, function(s) {
  rcogmod_betagate(n, mu = s[["mu"]], phi = s[["phi"]], pex = s[["pex"]], bex = s[["bex"]])
})
ll <- function(y, mu, phi, pex, bex) {
  mean(dcogmod_betagate(y, mu = mu, phi = phi, pex = pex, bex = bex, log = TRUE))
}
truth <- mapply(function(y, s) ll(y, s[["mu"]], s[["phi"]], s[["pex"]], s[["bex"]]),
                data, styles)

# `mu` and `bex` start at 0.5, where nothing favours a side; `phi` and `pex`
# move, log(2) and 0.5 being where init = 0 puts them.
grid <- expand.grid(phi = c(log(2), 1, 1.5, 2, 3, 4, 6), pex = c(0.5, 0.3, 0.2, 0.1, 0.05))
short <- t(apply(grid, 1, function(g) {
  truth - vapply(data, ll, numeric(1), mu = 0.5, phi = g[["phi"]], pex = g[["pex"]], bex = 0.5)
}))
res <- cbind(grid, round(short, 3), worst = round(apply(short, 1, max), 3),
             mean = round(rowMeans(short), 3))
res$phi <- round(res$phi, 2)
res <- res[order(res$worst), ]

# What brms does with no init: U(-2, 2) on every unconstrained coordinate.
set.seed(1)
random_worst <- replicate(2000, {
  u <- stats::runif(4, -2, 2)
  max(truth - vapply(data, ll, numeric(1), mu = stats::plogis(u[1]),
                     phi = log1p(exp(u[2])), pex = stats::plogis(u[3]),
                     bex = stats::plogis(u[4])))
})

dir.create("benchmarks/betagate_start/results", showWarnings = FALSE)
utils::write.csv(res, "benchmarks/betagate_start/results/shortfall.csv", row.names = FALSE)
options(width = 200)
print(res, row.names = FALSE)
cat("\ninit = 'random', worst shortfall over the styles, quantiles over 2000 draws:\n")
print(round(stats::quantile(random_worst, c(0.1, 0.5, 0.9)), 3))
