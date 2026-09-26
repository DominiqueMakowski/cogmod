# Prototype likelihoods for cogmod_lnr(), for bench.R. Not package code: it
# reads the package's internals (.choice_lpdf(), .LNR_STAN_PRELUDE,
# .POUTLIER_SCALE) and is sourced after pkgload::load_all(). See README.md.
#
# brms, with custom_family(loop = FALSE), writes ONE call
#   target += cogmod_lnr_lpdf(Y | mu, nuone, sigmazero, ..., dec);
# passing `vector` for every predicted dpar and `real` for aux / fixed ones,
# and the whole `dec` array (vars = "dec", not "dec[n]"). So the overload has
# to be generated for the formula: 2^6 signatures, `mu` always a vector.
# V0-V2 are that overload in three styles, each defined NEXT TO the package's
# scalar cogmod_lnr_lpdf() (Stan overloads user functions by signature), which
# they call for the sigmabias > 0 path. V3 is not vectorised at all.
#
# None of these guards against brms threading: with threads and loop = FALSE,
# brms slices Y and the dpars but passes the full `dec`, so a real
# implementation would need `if (size(dec) != rows(Y)) reject(...)`.

.LNR_ARGS <- c("mu", "nuone", "sigmazero", "sigmaone", "sigmabias", "ndt", "poutlier")

# TRUE = vector, per dpar, read off the formula the way brms decides it:
# predicted dpars are vectors, fixed (`sigmabias = 0`) and omitted ones reals.
lnr_arg_types <- function(f) {
  bt <- brms::brmsterms(f)
  pred <- names(bt$dpars)
  stats::setNames(.LNR_ARGS %in% pred, .LNR_ARGS)
}

lnr_vec_family <- function() {
  fam <- cogmod_lnr()
  fam$loop <- FALSE
  fam$vars <- "dec"
  fam
}

# ---- pieces ------------------------------------------------------------------

.sig <- function(is_vec) {
  paste(c("vector Y",
          sprintf("%s %s", ifelse(is_vec, "vector", "real"), names(is_vec)),
          "array[] int dec"), collapse = ", ")
}
.el <- function(is_vec, nm, idx) if (is_vec[[nm]]) sprintf("%s[%s]", nm, idx) else nm
.lo <- function(is_vec, nm) if (is_vec[[nm]]) sprintf("min(%s)", nm) else nm
.hi <- function(is_vec, nm) if (is_vec[[nm]]) sprintf("max(%s)", nm) else nm

.checks <- function(v) {
  sprintf(paste0(
    "  if (%s <= 0 || %s <= 0 || %s < 0 || %s < 0 || %s < 0 || %s > 1)\n",
    "    return negative_infinity();\n",
    "  if (min(dec) < 0 || max(dec) > 1) return negative_infinity();\n",
    "  if (min(Y) <= 0) return negative_infinity();\n"),
    .lo(v, "sigmazero"), .lo(v, "sigmaone"), .lo(v, "sigmabias"),
    .lo(v, "ndt"), .lo(v, "poutlier"), .hi(v, "poutlier"))
}

# Outlier constants, exactly as .choice_lpdf() writes them.
.lc <- function() formatC(log(2) - 0.5 * log(2 * pi * .POUTLIER_SCALE^2) - log(2),
                          format = "g", digits = 17, width = 1)
.k2 <- function() formatC(1 / (2 * .POUTLIER_SCALE^2), format = "g", digits = 15, width = 1)

# The scalar call, element n.
.scalar_call <- function(v, idx = "n") {
  sprintf("cogmod_lnr_lpdf(Y[%s] | %s, dec[%s])", idx,
          paste(vapply(names(v), function(nm) .el(v, nm, idx), ""), collapse = ", "), idx)
}

# ---- V0: the loop moved inside the function ------------------------------------

lnr_v0 <- function(v) {
  sprintf("
real cogmod_lnr_lpdf(%s) {
  int N = rows(Y);
  real lp = 0;
  for (n in 1:N) lp += %s;
  return lp;
}
", .sig(v), .scalar_call(v))
}

# ---- V1: one fused scalar loop, hoisted constants ---------------------------------
# Same per-element arithmetic as the scalar function at sigmabias = 0, in one
# loop with no nested user-function calls on the common path: log(t) is
# computed once for winner and loser, and log(poutlier) once per call when
# poutlier is a real.

lnr_v1 <- function(v) {
  pv <- v[["poutlier"]]
  lpo <- if (pv) "log(poutlier[n])" else "lpo"
  pre <- if (pv) "" else "  real lpo = log(poutlier);\n"
  sprintf("
real cogmod_lnr_lpdf(%s) {
  int N = rows(Y);
%s  if (%s > 0) {
    real lp = 0;
    for (n in 1:N) lp += %s;
    return lp;
  }
%s  real lp = 0;
  for (n in 1:N) {
    real t = Y[n] - %s;
    if (t <= 0) {
      lp += %s + (%s - %s * square(Y[n]));
    } else {
      real lt = log(t);
      real ld;
      if (dec[n] == 0) {
        ld = lognormal_lpdf(t | -%s, %s)
             + cogmod_log_Phi(-(lt + %s) / %s);
      } else {
        ld = lognormal_lpdf(t | -%s, %s)
             + cogmod_log_Phi(-(lt + %s) / %s);
      }
      lp += log_mix(%s, %s - %s * square(Y[n]), ld);
    }
  }
  return lp;
}
", .sig(v), .checks(v), .hi(v, "sigmabias"), .scalar_call(v), pre,
    .el(v, "ndt", "n"), lpo, .lc(), .k2(),
    .el(v, "mu", "n"), .el(v, "sigmazero", "n"), .el(v, "nuone", "n"), .el(v, "sigmaone", "n"),
    .el(v, "nuone", "n"), .el(v, "sigmaone", "n"), .el(v, "mu", "n"), .el(v, "sigmazero", "n"),
    .el(v, "poutlier", "n"), .lc(), .k2())
}

# ---- V2: vector expressions over index sets ----------------------------------------
# Elements are split once per call into three index sets from the VALUES (t <= 0,
# winner 0, winner 1); each set is then one block of vector arithmetic. The
# normal tail is cogmod_log_Phi() over a further split by branch, so the
# vector code carries no per-element branching at all.

.LOG_PHI_VEC <- "
// cogmod_log_Phi() over a vector, one branch per index set so that every
// element takes exactly the arithmetic the scalar function gives it.
vector cogmod_log_Phi_vec(vector x) {
  int N = rows(x);
  int na = 0;
  int nb = 0;
  for (n in 1:N) {
    if (x[n] < -25) na += 1;
    else if (x[n] > 0) nb += 1;
  }
  array[na] int ia;
  array[nb] int ib;
  array[N - na - nb] int ic;
  {
    int a = 1;
    int b = 1;
    int c = 1;
    for (n in 1:N) {
      if (x[n] < -25) { ia[a] = n; a += 1; }
      else if (x[n] > 0) { ib[b] = n; b += 1; }
      else { ic[c] = n; c += 1; }
    }
  }
  vector[N] out;
  if (na > 0) {
    vector[na] xa = x[ia];
    vector[na] z = inv_square(xa);
    out[ia] = -0.5 * square(xa) - log(-xa) - 0.91893853320467274
              + log(1 + z .* (-1 + z .* (3 + z .* (-15 + z .* (105 - 945 * z)))));
  }
  if (nb > 0) out[ib] = log1p(-0.5 * erfc(x[ib] * 0.7071067811865476));
  if (N - na - nb > 0) out[ic] = log(0.5 * erfc(-x[ic] * 0.7071067811865476));
  return out;
}
"

lnr_v2 <- function(v) {
  el <- function(nm, idx) .el(v, nm, idx)
  # elementwise divide for a vector dpar, scalar divide otherwise
  dv <- function(nm, idx) if (v[[nm]]) sprintf("./ %s[%s]", nm, idx) else sprintf("/ %s", nm)
  lg <- function(nm, idx) if (v[[nm]]) sprintf("log(%s[%s])", nm, idx) else sprintf("log(%s)", nm)
  block <- function(idx, n, w_nu, w_s, l_nu, l_s) {
    sprintf("
  if (%s > 0) {
    vector[%s] lt = log(t[%s]);
    vector[%s] ld = -lt - %s - 0.91893853320467274
                    - 0.5 * square((lt + %s) %s)
                    + cogmod_log_Phi_vec(-(lt + %s) %s);
    lp += sum(log_sum_exp(%s + (%s - %s * square(Y[%s])), %s + ld));
  }", n, n, idx, n, lg(w_s, idx), el(w_nu, idx), dv(w_s, idx),
      el(l_nu, idx), dv(l_s, idx),
      if (v[["poutlier"]]) sprintf("log(poutlier[%s])", idx) else "log(poutlier)",
      .lc(), .k2(), idx,
      if (v[["poutlier"]]) sprintf("log1m(poutlier[%s])", idx) else "log1m(poutlier)")
  }
  sprintf("
real cogmod_lnr_lpdf(%s) {
  int N = rows(Y);
%s  if (%s > 0) {
    real lp = 0;
    for (n in 1:N) lp += %s;
    return lp;
  }
  vector[N] t = Y - ndt;
  int nn = 0;
  int n0 = 0;
  for (n in 1:N) {
    if (t[n] <= 0) nn += 1;
    else if (dec[n] == 0) n0 += 1;
  }
  int n1 = N - nn - n0;
  array[nn] int ineg;
  array[n0] int i0;
  array[n1] int i1;
  {
    int a = 1;
    int b = 1;
    int c = 1;
    for (n in 1:N) {
      if (t[n] <= 0) { ineg[a] = n; a += 1; }
      else if (dec[n] == 0) { i0[b] = n; b += 1; }
      else { i1[c] = n; c += 1; }
    }
  }
  real lp = 0;
  if (nn > 0) lp += %s + sum(%s - %s * square(Y[ineg]));%s%s
  return lp;
}
", .sig(v), .checks(v), .hi(v, "sigmabias"), .scalar_call(v),
    if (v[["poutlier"]]) "sum(log(poutlier[ineg]))" else "nn * log(poutlier)",
    .lc(), .k2(),
    block("i0", "n0", "mu", "sigmazero", "nuone", "sigmaone"),
    block("i1", "n1", "nuone", "sigmaone", "mu", "sigmazero"))
}

# ---- stanvars -----------------------------------------------------------------

# The package's scalar function (prelude + cogmod_lnr_lpdf for one observation)
# followed by the vectorised overload of `style`.
lnr_vec_stanvars <- function(f, style = c("v0", "v1", "v2")) {
  style <- match.arg(style)
  v <- lnr_arg_types(f)
  vec <- switch(style, v0 = lnr_v0(v), v1 = lnr_v1(v),
                v2 = paste0(.LOG_PHI_VEC, lnr_v2(v)))
  brms::stanvar(scode = paste0(.choice_lpdf("cogmod_lnr"), vec), block = "functions")
}

# ---- V3: NOT vectorised - the looped family, with a leaner per-trial density ----
# Same brms code as the base (loop = TRUE); only the scalar lpdf changes. At
# sigmabias = 0 the winner is one lognormal_lpdf() and the loser one
# lognormal_lccdf() - fused built-ins, one tape entry each with analytic
# partials - instead of ~9 separate operations through cogmod_log_Phi().
# lognormal_lccdf() IS log(0.5 * erfc(z / sqrt(2))), cogmod_log_Phi()'s body,
# and fails only where erfc underflows (z ~ 37.5), so it is guarded by the
# same switch to the asymptotic series at z = 25 that cogmod_log_Phi() uses.
lnr_v3_stanvars <- function() {
  lean <- sprintf("
real cogmod_lnr_lpdf(real Y, real mu, real nuone, real sigmazero, real sigmaone, real sigmabias, real ndt, real poutlier, int dec) {
    if (sigmazero <= 0 || sigmaone <= 0 || sigmabias < 0 || ndt < 0 || poutlier < 0 || poutlier > 1) {
      return negative_infinity();
    }
    if (dec < 0 || dec > 1) return negative_infinity();
    if (Y <= 0) return negative_infinity();
    real t_adj = Y - ndt;
    if (t_adj <= 0) return log(poutlier) + (%s - %s * square(Y));
    real lp_dec;
    if (sigmabias != 0) {
      lp_dec = dec == 0
        ? cogmod_lnr_decision_lpdf(t_adj | mu, sigmazero, nuone, sigmaone, sigmabias)
        : cogmod_lnr_decision_lpdf(t_adj | nuone, sigmaone, mu, sigmazero, sigmabias);
    } else {
      real lt = log(t_adj);
      if (dec == 0) {
        real z = (lt + nuone) / sigmaone;
        lp_dec = lognormal_lpdf(t_adj | -mu, sigmazero)
                 + (z < 25 ? lognormal_lccdf(t_adj | -nuone, sigmaone) : cogmod_log_Phi(-z));
      } else {
        real z = (lt + mu) / sigmazero;
        lp_dec = lognormal_lpdf(t_adj | -nuone, sigmaone)
                 + (z < 25 ? lognormal_lccdf(t_adj | -mu, sigmazero) : cogmod_log_Phi(-z));
      }
    }
    return log_mix(poutlier, %s - %s * square(Y), lp_dec);
}
", .lc(), .k2(), .lc(), .k2())
  brms::stanvar(scode = paste0(get(".LNR_STAN_PRELUDE"), lean), block = "functions")
}
