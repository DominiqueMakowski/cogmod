# Shared set-up for the ex-Wald mode investigation (README.md). Sourced by the
# other scripts, from the package root:
#
#   EXW_TREE  package tree to load_all() (default "."). The likelihood is the
#             same in 0.3.4 and 0.3.5; only the priors differ, and those are
#             written out below rather than read off the tree.
#   EXW_OUT   where results go (default benchmarks/exwald_modes/results)
#
# What it provides: the RT article's data and formula, an R twin of the
# brms program's lp__ on Stan's unconstrained scale for each prior set, and
# helpers to read the four modes off a parameter vector.
suppressMessages(pkgload::load_all(Sys.getenv("EXW_TREE", "."), quiet = TRUE))
OUT <- Sys.getenv("EXW_OUT", "benchmarks/exwald_modes/results")
dir.create(OUT, recursive = TRUE, showWarnings = FALSE)

# Three participants of rtdists' speed_acc, correct responses, RT <= 2 s:
# exactly what vignettes/models/m_wald4.rds and m_exwald.rds were fitted to.
D <- as.data.frame(readRDS("vignettes/models/m_wald4.rds")$data)[, c("RT", "Condition")]
FORMULA <- brms::bf(RT ~ Condition, boundary ~ Condition, ndt ~ Condition,
                    tau ~ Condition, family = cogmod_exwald())
Y <- D$RT
XC <- as.numeric(D$Condition == "Speed")
MX <- mean(XC)        # brms centres the Speed dummy: Accuracy -MX, Speed 1 - MX
XC <- XC - MX

# Stan's declaration order, which is also the order of every parameter vector
# here. poutlier enters on the logit scale, as Stan samples it.
PN <- c("b", "Intercept", "b_boundary", "Intercept_boundary", "b_tau",
        "Intercept_tau", "b_ndt", "Intercept_ndt", "poutlier_u")
softplus <- function(x) ifelse(x > 30, x, log1p(exp(x)))
softplus_inv <- function(x) log(expm1(x))

# The priors of each set, as the generated programs state them. 0.3.4 is
# 5382c79 (arm A of benchmarks/driftone_priors/); 0.3.5 adds the boundary rows
# and moves poutlier's rate from 100 to 280; P1 is 0.3.5 plus the drift prior
# this investigation proposes (README.md, "A fix").
lprior_fun <- function(u, prior = "0.3.5") {
  p <- stats::plogis(u[9])
  lp <- if (prior == "P1") {
    stats::dnorm(u[2], 3, 1, log = TRUE) + stats::dnorm(u[1], 0, 1, log = TRUE)
  } else {
    stats::dt(u[2] / 2.5, 3, log = TRUE) - log(2.5)    # brms's default, slope flat
  }
  lp <- lp + stats::dnorm(u[5], 0, 0.5, log = TRUE) + stats::dnorm(u[6], -1.5, 0.7, log = TRUE) +
    stats::dnorm(u[7], 0, 0.2, log = TRUE) + stats::dnorm(u[8], -1.2, 0.5, log = TRUE)
  if (prior == "0.3.4") {
    lp + stats::dexp(p, 100, log = TRUE) - stats::pexp(1, 100, log.p = TRUE)
  } else {
    lp + stats::dnorm(u[3], 0, 0.5, log = TRUE) + stats::dnorm(u[4], 0.5, 1, log = TRUE) +
      stats::dexp(p, 280, log = TRUE) - stats::pexp(1, 280, log.p = TRUE)
  }
}

loglik_fun <- function(u) {
  sum(dcogmod_exwald(Y, mu = softplus(u[2] + XC * u[1]), boundary = softplus(u[4] + XC * u[3]),
                     tau = softplus(u[6] + XC * u[5]), ndt = exp(u[8] + XC * u[7]),
                     poutlier = stats::plogis(u[9]), log = TRUE))
}

# lp__ itself: brms writes every term with target +=, so constants included,
# plus the log Jacobian of poutlier's logit transform. Matches Stan's lp__ to
# 2e-11 (numerics.R).
lp_fun <- function(u, prior = "0.3.5") {
  p <- stats::plogis(u[9])
  loglik_fun(u) + lprior_fun(u, prior) + log(p) + log1p(-p)
}

# Natural parameters of each condition
natpar <- function(u) {
  x <- c(Accuracy = -MX, Speed = 1 - MX)
  data.frame(cond = names(x), mu = softplus(u[2] + x * u[1]),
             boundary = softplus(u[4] + x * u[3]), tau = softplus(u[6] + x * u[5]),
             ndt = exp(u[8] + x * u[7]), poutlier = stats::plogis(u[9]), row.names = NULL)
}

# W where the Wald's drift sets the right tail (mu^2 / 2 < 1 / tau, the
# Faddeeva branch), E where the exponential does (the closed-form branch);
# Accuracy first. The seam between the branches, mu^2 = 2 / tau, is exactly
# where the two tail rates are equal.
mode_label <- function(u) {
  n <- natpar(u)
  paste0(ifelse(n$mu^2 - 2 / n$tau > 0, "E", "W"), collapse = "")
}

# Parameter vectors <-> the init lists cmdstanr takes
as_init <- function(u) {
  list(b = u[1], Intercept = u[2], b_boundary = u[3], Intercept_boundary = u[4],
       b_tau = u[5], Intercept_tau = u[6], b_ndt = u[7], Intercept_ndt = u[8],
       poutlier = stats::plogis(u[9]))
}
draws_to_u <- function(x) {
  u <- cbind(x$`b[1]`, x$Intercept, x$`b_boundary[1]`, x$Intercept_boundary, x$`b_tau[1]`,
             x$Intercept_tau, x$`b_ndt[1]`, x$Intercept_ndt, stats::qlogis(x$poutlier))
  colnames(u) <- PN
  u
}

# The program for a prior set, built the way brm() would build it, from the
# 0.3.5 tree. 0.3.4 is that with the boundary rows emptied (brms then leaves
# them flat, as 0.3.4 did) and poutlier's old rate; P1 replaces the
# response's intercept row and fills its slope row, which is all the proposal
# changes. numerics.R checks the 0.3.4 program against 5382c79's own.
make_program <- function(prior = "0.3.5") {
  pri <- suppressMessages(suppressWarnings(cogmod_priors(FORMULA, D)))
  if (prior == "0.3.4") {
    pri$prior[pri$dpar == "boundary" & pri$class %in% c("Intercept", "b")] <- ""
    pri$prior[pri$class == "poutlier"] <- "exponential(100)"
  }
  if (prior == "P1") {
    pri$prior[pri$class == "Intercept" & !nzchar(pri$dpar)] <- "normal(3, 1)"
    row <- pri$class == "b" & !nzchar(pri$dpar) & !nzchar(pri$coef)
    if (any(row)) pri$prior[row] <- "normal(0, 1)" else
      pri <- pri + brms::set_prior("normal(0, 1)", class = "b")
  }
  sv <- cogmod_stanvars(FORMULA)
  list(prior = pri,
       code = suppressWarnings(brms::make_stancode(FORMULA, data = D, prior = pri, stanvars = sv)),
       sdata = unclass(suppressWarnings(brms::make_standata(FORMULA, data = D, prior = pri, stanvars = sv))))
}
