# Example 2 of the paper: bounded ratings with betagate().
# Run from the `paper/` directory. Reproduces @tbl-betagate.
#
# The manipulation raises both the underlying tendency (`mu`) and the
# inclination to answer at the top of the scale (`pex`, `bex`), so the mean
# rating moves for two distinct reasons. A mean comparison reports one number;
# Beta-Gate separates them. Consistency (`phi`) is left unchanged as a
# negative control.

library(cogmod)
library(brms)

set.seed(7)

n <- 500
true <- list(
  A = list(mu = 0.55, phi = 3, pex = 0.10, bex = 0.50),
  B = list(mu = 0.70, phi = 3, pex = 0.35, bex = 0.80)
)

df <- do.call(rbind, lapply(names(true), function(k) {
  p <- true[[k]]
  data.frame(
    Condition = k,
    score = rbetagate(n, mu = p$mu, phi = p$phi, pex = p$pex, bex = p$bex)
  )
}))
df$Condition <- factor(df$Condition)

# What a mean comparison sees.
print(aggregate(score ~ Condition, df, function(v) {
  round(c(mean = mean(v), prop_0 = mean(v == 0), prop_1 = mean(v == 1)), 3)
}))
print(t.test(score ~ Condition, data = df))

f <- bf(
  score ~ Condition,
  phi ~ Condition,
  pex ~ Condition,
  bex ~ Condition,
  family = betagate()
)

m <- brm(
  f,
  data = df,
  stanvars = betagate_stanvars(),
  init = 0,
  chains = 4, cores = 4, iter = 2000, seed = 1,
  control = list(adapt_delta = 0.95),
  backend = "cmdstanr",
  file = "m_betagate_ratings", file_refit = "on_change"
)
print(m)

# Response scale: logit link on mu/pex/bex, softplus on phi.
d <- brms::as_draws_df(m)
softplus <- function(x) log1p(exp(x))
cell <- function(v) {
  sprintf("%.2f [%.2f, %.2f]", mean(v), quantile(v, .025), quantile(v, .975))
}

pars <- list(
  mu  = list(link = plogis,   int = "b_Intercept",     slope = "b_ConditionB"),
  phi = list(link = softplus, int = "b_phi_Intercept", slope = "b_phi_ConditionB"),
  pex = list(link = plogis,   int = "b_pex_Intercept", slope = "b_pex_ConditionB"),
  bex = list(link = plogis,   int = "b_bex_Intercept", slope = "b_bex_ConditionB")
)

for (nm in names(pars)) {
  p <- pars[[nm]]
  a <- p$link(d[[p$int]])
  b <- p$link(d[[p$int]] + d[[p$slope]])
  cat(sprintf(
    "%-4s A = %-20s B = %-20s  (generated: %.2f -> %.2f)\n",
    nm, cell(a), cell(b), true$A[[nm]], true$B[[nm]]
  ))
}

cat("\nDone.\n")
