# Reading the family off what the user passed
# ===========================================
#
# cogmod_priors(), cogmod_inits() and cogmod_stanvars() all take the model
# rather than the family, so the family is named once - in bf() - and cannot
# drift out of step with the rest of the call. These pull it back out.

# A brmsformula, a brmsfit and a family object are all lists carrying the
# family under `$family`; a family object is its own answer. A plain formula
# carries nothing, and gives NULL so the caller can say so.
#' @keywords internal
.cogmod_family <- function(x) {
  if (inherits(x, c("customfamily", "brmsfamily", "family"))) return(x)
  if (is.list(x)) return(x[["family"]])
  NULL
}

# brms families carry $name; stats families only $family.
#' @keywords internal
.family_name <- function(family) {
  if (is.null(family)) return(NULL)
  if (is.character(family)) return(family)
  nm <- family[["name"]]
  if (is.null(nm)) nm <- family[["family"]]
  nm
}

# The link of every dpar, named. brms stores the first dpar's link as `$link`
# and the rest as `$link_<dpar>`, so the two forms are folded back together
# here. Reading them off the family rather than the registry is what makes
# cogmod_gamma(link_mu = "log") behave.
#' @keywords internal
.family_links <- function(family) {
  dpars <- family[["dpars"]]
  if (is.null(dpars)) return(stats::setNames(character(0), character(0)))
  out <- vapply(seq_along(dpars), function(i) {
    v <- family[[paste0("link_", dpars[i])]]
    if (is.null(v) && i == 1L) v <- family[["link"]]
    if (is.null(v)) NA_character_ else as.character(v)[1]
  }, character(1))
  stats::setNames(out, dpars)
}


# How far one unit of a smooth's SD moves the linear predictor
# ============================================================
#
# brms writes the penalised part of a smooth as random effects: it adds
# `sds * Zs %*% zs` to the linear predictor with zs ~ N(0, 1), which has SD
# `sds * ||Zs[i, ]||` at row i. The size of ||Zs[i, ]|| is set by the basis mgcv
# builds, and it varies a lot. Measured with covariates uniform on [-1, 1]
# (2026-10-01), the root mean square over rows is about 0.45 for s(x), 0.26 for
# s(x, k = 5), 0.6 for s(x, y), 1.75 for s(x, bs = "cr", k = 5), 6 for
# s(x, bs = "cr") and 0.08-0.12 for t2(x, z, k = c(5, 5), bs = c("cr", "cr")).
# Covariates on a coarse grid or skewed move each of these by a factor of up
# to two. The number of rows barely matters (within 15% from 1,000 to 100,000)
# and the units of the covariates not at all. So a fixed number for `sds`, as
# a prior or as a start, is a different amount of wiggle on every basis: 70
# times more on that s(x, bs = "cr") than on that t2. Dividing by this scale
# converts it to link units, which is what cogmod_priors() and cogmod_inits()
# do.
#
# One value per Zs block, i.e. per penalty, and per level of a `by` factor,
# which brms numbers as a term of its own. `stem` picks out the blocks behind
# one `sds` parameter ("1" for mu's first smooth, "ndt_1" for ndt's) in penalty
# order. NULL takes every block. A block that gives no positive finite scale
# comes back NA.
#' @keywords internal
.zs_scales <- function(sdata, stem = NULL) {
  nm <- names(sdata)
  if (is.null(stem)) {
    blocks <- nm[startsWith(nm, "Zs_")]
  } else if (paste0("Zs_", stem) %in% nm) {
    # Older brms declared one scalar `sds_<k>_<j>` per penalty, which names its
    # block exactly.
    blocks <- paste0("Zs_", stem)
  } else {
    prefix <- paste0("Zs_", stem, "_")
    rest <- substring(nm, nchar(prefix) + 1L)
    hit <- startsWith(nm, prefix) & grepl("^[0-9]+$", rest)
    blocks <- nm[hit][order(as.integer(rest[hit]))]
  }
  vapply(blocks, function(b) {
    z <- as.matrix(sdata[[b]])
    s <- if (length(z)) sqrt(mean(rowSums(z^2))) else NA_real_
    if (is.finite(s) && s > 0) s else NA_real_
  }, numeric(1))
}


# The scales of one smooth term, as get_prior() labels it in `coef`, from a
# model holding that term alone. brms builds each smooth's basis on its own,
# so the blocks come out as they do in the full model, give or take rows the
# full model drops for missing values in other variables. Building the term
# alone avoids reconstructing which `Zs_<k>` blocks belong to it: brms numbers
# them in term order with a `by` factor expanded, while get_prior() sorts its
# rows by label. On the IGC Muller-Lyer data (324,000 rows) a t2() takes 1.8 s.
#
# Two scales come back, or NULL when the term cannot be built:
#
# `sds`: the geometric mean over the term's Zs blocks of .zs_scales(), because
# brms takes one sds prior per term for all of them (the three penalties of
# the t2() above differ by about 20%). NA if any block has no scale.
#
# `bs`: the root mean square over rows of each column of the term's
# UNPENALISED part, `Xs` - the linear trend of an s(x) or the marginal trends
# of a t2() - named as brms names the columns and as get_prior() labels the
# matching `b` rows (`sx_1`, `t2xw_1`, `sw:ga_1`). mgcv scales these columns
# without regard to the covariate's units and differently on every basis:
# measured with x on [-1, 1] (2026-10-02), 0.16 for s(x) and for s(age) alike,
# 0.11 per column for s(x, w), 0.8-1.0 for a t2(..., bs = c("cr", "cr")),
# 1.0 for s(x, bs = "cr", k = 5) and 2.0 for s(x, bs = "cr"). So a prior
# stated on the coefficient means a different amount of trend on every basis,
# twentyfold between the ends, and for a tp or cr basis the penalised part has
# no linear component to make up the difference. A column whose scale is not
# a positive finite number is NA.
#' @keywords internal
.smooth_term_scales <- function(term, formula, data, ...) {
  f <- if (inherits(formula, "brmsformula")) formula$formula else formula
  env <- if (inherits(f, "formula")) environment(f) else NULL
  if (is.null(env)) env <- globalenv()
  y <- ".cogmod_scale_y"
  one <- tryCatch(stats::as.formula(paste(y, "~", term), env = env),
                  error = function(e) NULL)
  if (is.null(one)) return(NULL)
  data[[y]] <- 0
  sdata <- tryCatch(
    suppressWarnings(brms::make_standata(one, data = data, ...)),
    error = function(e) NULL
  )
  if (is.null(sdata)) return(NULL)
  s <- .zs_scales(sdata)
  sds <- if (!length(s) || anyNA(s)) NA_real_ else exp(mean(log(s)))
  Xs <- sdata$Xs
  bs <- if (is.null(Xs) || !length(Xs)) {
    numeric(0)
  } else {
    Xs <- as.matrix(Xs)
    r <- sqrt(colMeans(Xs^2))
    r[!is.finite(r) | r <= 0] <- NA_real_
    r
  }
  list(sds = sds, bs = bs)
}


# Warning when the evidence scale is left free
# ===========================================
#
# The LBA families have a likelihood that is EXACTLY constant along a ray:
# multiply every parameter named in the registry's `scale_ray` slot by any
# c > 0 and every finishing time `(b - z) / v` is unchanged. Nothing in the
# data can distinguish a point on that ray from any other, so unless one of
# them is pinned in `bf()`, what comes back for the individual parameters is
# whatever the priors happen to say about that direction. The RT distribution
# they jointly generate is perfectly well identified either way, which is
# exactly what makes this quiet: the fit converges, `pp_check()` looks right,
# and only the parameter estimates are meaningless.
#
# Fixing any ONE member of the ray at a NON-ZERO value pins it - `sigmazero = 1`
# conventionally, but `boundary = 1` or `sigmabias = 0.5` do the job just as
# well. Two things that look like they should pin it do not:
#
#   - Omitting a dpar from `bf()`. brms then declares it as a free auxiliary
#     parameter, which leaves the ray exactly as free as before. That
#     distinction is why this looks at `$pfix` rather than at `$pforms`.
#   - Fixing one at ZERO. Scaling multiplies every member by c and c * 0 = 0, so
#     zero is the one value the ray runs through unchanged. `sigmabias = 0` is
#     the recinormal (LATER), and it is a genuinely useful thing to write - but
#     it removes a parameter from the ray rather than pinning the ray, and the
#     remaining members still need one of their own.
#
# `.warn_scale_ray()` is called from cogmod_stanvars(), which is the one of the
# three generics a fit cannot skip.
#' @keywords internal
.warn_scale_ray <- function(formula) {
  if (!inherits(formula, "brmsformula")) return(invisible(NULL))
  fam <- .family_name(.cogmod_family(formula))
  if (is.null(fam)) return(invisible(NULL))
  spec <- if (is.null(.SHIFTED[[fam]])) .CHOICE[[fam]] else .SHIFTED[[fam]]
  ray <- spec[["scale_ray"]]
  if (is.null(ray)) return(invisible(NULL))
  # `mu` is the response's own formula, so it is never in $pfix and cannot be
  # the one the user pinned. It stays in `ray` because it is genuinely on the
  # ray - it just is not a way out of it.
  pinned <- formula$pfix[intersect(names(formula$pfix), ray)]
  # Anything that is not a single number counts as a pin, which is what this
  # did before the zero case was carved out: only an unambiguous zero is
  # treated as failing to pin.
  pins <- vapply(pinned, function(v) {
    num <- tryCatch(suppressWarnings(as.numeric(v)),
                    error = function(e) NA_real_)
    length(num) != 1 || is.na(num) || num != 0
  }, logical(1))
  if (any(pins)) return(invisible(NULL))

  # The message spells the invariance out rather than naming it. "The evidence
  # scale is arbitrary" means nothing to someone who has not already worked out
  # what is wrong; "multiplying these four by a common constant leaves the
  # likelihood unchanged" is the same fact in a form the reader can check.
  sd <- if (identical(fam, "cogmod_lba1")) "sigma" else "sigmazero"
  zero <- names(pins)[!pins]
  warning(
    fam, "(): the evidence scale is arbitrary - multiplying ",
    paste0("`", ray, "`", collapse = ", "),
    " by a common constant leaves the likelihood exactly unchanged - so one of",
    " them has to be fixed to a non-zero value. Add `", sd,
    " = 1` to the formula.",
    if (length(zero)) {
      paste0(" `", zero[1], " = 0` does not count: zero times any constant is",
             " still zero, so fixing it constrains nothing and the others stay",
             " free to scale.")
    } else "",
    " See ?", sub("^cogmod_", "rcogmod_", fam), ".",
    call. = FALSE
  )
  invisible(NULL)
}


# Stable computation of log(1 - exp(x))
# A stable version of log1m_exp (assumes x is scalar or vector)
#' @keywords internal
.log1m_exp <- function(x) {
  # Computes log(1 - exp(x)) in a numerically stable way.
  # x should be negative; for vectorized input, use ifelse.
  ifelse(x < log(0.5), log1p(-exp(x)), log(-expm1(x)))
}






# Stable log of a two-component mixture, mirroring Stan's log_mix():
#   log(theta * exp(lp1) + (1 - theta) * exp(lp2))
# Handles lp2 = -Inf (a component with zero density), which is the normal case
# for shifted distributions evaluated below their shift.
#' @keywords internal
.log_mix <- function(theta, lp1, lp2) {
  a <- log(theta) + lp1
  b <- log1p(-theta) + lp2
  m <- pmax(a, b)
  out <- m + log(exp(a - m) + exp(b - m))
  # Both components vanish; pmax() is -Inf and the shift would give NaN.
  out[!is.finite(m)] <- -Inf
  out
}


# Replace every input the density cannot be evaluated at with one it can, and
# report which those were, so the caller can overwrite their result with -Inf.
#
# This exists because "cannot be evaluated at" and "returns zero density at" are
# not the same thing here. The density cores branch on comparisons - the Wald's
# `if (any(A / sqrt(t) < .RDM_EPS_A))`, the LBA's and the DDM's equivalents -
# and a comparison against NA is NA, which makes a scalar `if()` throw rather
# than return. So a missing reaction time, or a missing parameter, aborted the
# whole call instead of producing the zero that `.ldec()` and `.ldec_choice()`
# both go on to write for it. Measured against the sixteen mixture families:
# four threw on a missing time (cogmod_lba1, cogmod_lba2, cogmod_rdm,
# cogmod_ddm) and two on an infinite one (cogmod_loggamma, cogmod_lba2); the
# rest returned the zero already, their densities being arithmetic all the way
# down.
#
# The substitute for a bad parameter is `spec$init`, the registry's own starting
# point for that family - by construction a set the density accepts. The
# substitute for a bad time is 1, for the same reason and no other. Both results
# are discarded.
#' @keywords internal
.dens_mask <- function(spec, t, p) {
  ok <- is.finite(t) & t > 0
  for (d in spec$dpars) {
    v <- p[[d]]
    nf <- !is.finite(v)
    if (any(nf)) {
      # `ok[] <-` rather than `ok <-`: `t` is a draws x observations matrix when
      # p_outlier() is the caller, and the dim has to survive.
      ok[] <- ok & rep_len(!nf, length(ok))
      v[nf] <- spec$init[[d]]
      p[[d]] <- v
    }
  }
  if (!all(ok)) t[!ok] <- 1
  list(t = t, p = p, ok = ok)
}


# Rejection sampling algorithm by Robert (Stat. Comp (1995), 5, 121-5)
# for simulating from the truncated normal distribution.
# Copied from the msm package
#' @keywords internal
.rnorm_truncated <- function(n, mean = 0, sd = 1, lower = -Inf, upper = Inf) {
  if (length(n) > 1) {
    n <- length(n)
  }
  mean <- rep(mean, length = n)
  sd <- rep(sd, length = n)
  lower <- rep(lower, length = n)
  upper <- rep(upper, length = n)
  ret <- numeric(n)
  ind <- seq(length.out = n)

  sdzero <- sd < .Machine$double.eps
  ## return the mean, unless mean is outside the range, then return nan
  sdna <- sdzero & ((mean < lower) | (mean > upper))

  lower <- (lower - mean) / sd ## Algorithm works on mean 0, sd 1 scale
  upper <- (upper - mean) / sd
  nas <- is.na(mean) | is.na(sd) | is.na(lower) | is.na(upper) | sdna
  if (any(nas)) warning("NAs produced")
  ## Different algorithms depending on where upper/lower limits lie.
  alg <- ifelse(
    ((lower > upper) | nas),
    -1, # return NaN
    ifelse(
      sdzero,
      4, # SD zero, so set the sampled value to the mean.
      ifelse(
        ((lower < 0 & upper == Inf) |
          (lower == -Inf & upper > 0) |
          (is.finite(lower) & is.finite(upper) & (lower < 0) & (upper > 0) & (upper - lower > sqrt(2 * pi)))
        ),
        0, # standard "simulate from normal and reject if outside limits" method. Use if bounds are wide.
        ifelse(
          (lower >= 0 & (upper > lower + 2 * sqrt(exp(1)) /
            (lower + sqrt(lower^2 + 4)) * exp((lower * 2 - lower * sqrt(lower^2 + 4)) / 4))),
          1, # rejection sampling with exponential proposal. Use if lower >> mean
          ifelse(upper <= 0 & (-lower > -upper + 2 * sqrt(exp(1)) /
            (-upper + sqrt(upper^2 + 4)) * exp((upper * 2 - -upper * sqrt(upper^2 + 4)) / 4)),
          2, # rejection sampling with exponential proposal. Use if upper << mean.
          3
          )
        )
      )
    )
  ) # rejection sampling with uniform proposal. Use if bounds are narrow and central.

  ind.nan <- ind[alg == -1]
  ind.no <- ind[alg == 0]
  ind.expl <- ind[alg == 1]
  ind.expu <- ind[alg == 2]
  ind.u <- ind[alg == 3]
  ind.sd0 <- ind[alg == 4]
  ret[ind.nan] <- NaN
  ret[ind.sd0] <- 0 # SD zero, so set the sampled value to the mean.
  while (length(ind.no) > 0) {
    y <- stats::rnorm(length(ind.no))
    done <- which(y >= lower[ind.no] & y <= upper[ind.no])
    ret[ind.no[done]] <- y[done]
    ind.no <- setdiff(ind.no, ind.no[done])
  }
  stopifnot(length(ind.no) == 0)
  while (length(ind.expl) > 0) {
    a <- (lower[ind.expl] + sqrt(lower[ind.expl]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expl), a) + lower[ind.expl]
    u <- stats::runif(length(ind.expl))
    done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= upper[ind.expl]))
    ret[ind.expl[done]] <- z[done]
    ind.expl <- setdiff(ind.expl, ind.expl[done])
  }
  stopifnot(length(ind.expl) == 0)
  while (length(ind.expu) > 0) {
    a <- (-upper[ind.expu] + sqrt(upper[ind.expu]^2 + 4)) / 2
    z <- stats::rexp(length(ind.expu), a) - upper[ind.expu]
    u <- stats::runif(length(ind.expu))
    done <- which((u <= exp(-(z - a)^2 / 2)) & (z <= -lower[ind.expu]))
    ret[ind.expu[done]] <- -z[done]
    ind.expu <- setdiff(ind.expu, ind.expu[done])
  }
  stopifnot(length(ind.expu) == 0)
  while (length(ind.u) > 0) {
    z <- stats::runif(length(ind.u), lower[ind.u], upper[ind.u])
    rho <- ifelse(lower[ind.u] > 0,
      exp((lower[ind.u]^2 - z^2) / 2), ifelse(upper[ind.u] < 0,
        exp((upper[ind.u]^2 - z^2) / 2),
        exp(-z^2 / 2)
      )
    )
    u <- stats::runif(length(ind.u))
    done <- which(u <= rho)
    ret[ind.u[done]] <- z[done]
    ind.u <- setdiff(ind.u, ind.u[done])
  }
  stopifnot(length(ind.u) == 0)
  ret * sd + mean
}
