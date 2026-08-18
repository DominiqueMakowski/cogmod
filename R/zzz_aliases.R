# Deprecated names
# ===============
#
# Every export the package had before the families were renamed to the
# `cogmod_*` scheme, bound to whatever it is called now. Plain bindings rather
# than wrappers, so the signature, the documentation and `args()` are the new
# function's - `rt_lognormal` and `cogmod_lognormal` are the same object.
#
# The brms hooks matter as much as the constructors: a model fitted before the
# rename carries `family$name == "rt_lognormal"`, and brms looks up
# `log_lik_rt_lognormal()` by that name when post-processing it. Without these,
# an old fit could still be loaded but not summarised.
#
# The file is named to collate last: these are assignments, not function
# bodies, so the right-hand sides have to exist by the time it is sourced.

#' Deprecated family and function names
#'
#' @description
#' The model families were renamed to a single `cogmod_*` scheme, so that every
#' family reads as one of this package's and none of them can collide with a
#' `brms` or `stats` name. The old names still work and are exact synonyms, but
#' they are no longer documented and will be removed in a future release.
#'
#' | old | new |
#' | --- | --- |
#' | `rt_lognormal()` | `cogmod_lognormal()` |
#' | `rt_loggamma()` | `cogmod_loggamma()` |
#' | `rt_invgaussian()` | `cogmod_invgaussian()` |
#' | `rt_gamma()` | `cogmod_gamma()` |
#' | `rt_invgamma()` | `cogmod_invgamma()` |
#' | `rt_weibull()` | `cogmod_weibull()` |
#' | `rt_invweibull()` | `cogmod_invweibull()` |
#' | `rt_logweibull()` | `cogmod_logweibull()` |
#' | `rt_exgaussian()` | `cogmod_exgaussian()` |
#' | `rt_lba()` | `cogmod_lba1()` |
#' | `lba()` | `cogmod_lba2()` |
#' | `lnr()` | `cogmod_lnr()` |
#' | `rdm()` | `cogmod_rdm()` |
#' | `ddm()` | `cogmod_ddm()` |
#' | `choco()` | `cogmod_choco()` |
#' | `betagate()` | `cogmod_betagate()` |
#' | `betadiscrete()` | `cogmod_betadiscrete()` |
#'
#' Every derived function follows its family: `rrt_lognormal()` is
#' `rcogmod_lognormal()`, `rt_lognormal_stanvars()` is
#' `cogmod_lognormal_stanvars()`, and so on for the density, the `*_stanvars()`,
#' the `*_lpdf_expose()` and the `brms` post-processing hooks.
#'
#' @name cogmod-deprecated
#' @keywords internal
NULL


# cogmod_lognormal ------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_lognormal <- dcogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_lognormal <- log_lik_cogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_lognormal <- posterior_epred_cogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_lognormal <- posterior_predict_cogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
rrt_lognormal <- rcogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
rt_lognormal <- cogmod_lognormal

#' @rdname cogmod-deprecated
#' @export
rt_lognormal_lpdf_expose <- cogmod_lognormal_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_lognormal_stanvars <- cogmod_lognormal_stanvars

# cogmod_loggamma -------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_loggamma <- dcogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_loggamma <- log_lik_cogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_loggamma <- posterior_epred_cogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_loggamma <- posterior_predict_cogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
rrt_loggamma <- rcogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
rt_loggamma <- cogmod_loggamma

#' @rdname cogmod-deprecated
#' @export
rt_loggamma_lpdf_expose <- cogmod_loggamma_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_loggamma_stanvars <- cogmod_loggamma_stanvars

# cogmod_invgaussian ----------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_invgaussian <- dcogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_invgaussian <- log_lik_cogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_invgaussian <- posterior_epred_cogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_invgaussian <- posterior_predict_cogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
prt_invgaussian <- pcogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
rrt_invgaussian <- rcogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
rt_invgaussian <- cogmod_invgaussian

#' @rdname cogmod-deprecated
#' @export
rt_invgaussian_lpdf_expose <- cogmod_invgaussian_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_invgaussian_stanvars <- cogmod_invgaussian_stanvars

# cogmod_gamma ----------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_gamma <- dcogmod_gamma

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_gamma <- log_lik_cogmod_gamma

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_gamma <- posterior_epred_cogmod_gamma

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_gamma <- posterior_predict_cogmod_gamma

#' @rdname cogmod-deprecated
#' @export
rrt_gamma <- rcogmod_gamma

#' @rdname cogmod-deprecated
#' @export
rt_gamma <- cogmod_gamma

#' @rdname cogmod-deprecated
#' @export
rt_gamma_lpdf_expose <- cogmod_gamma_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_gamma_stanvars <- cogmod_gamma_stanvars

# cogmod_invgamma -------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_invgamma <- dcogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_invgamma <- log_lik_cogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_invgamma <- posterior_epred_cogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_invgamma <- posterior_predict_cogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
rrt_invgamma <- rcogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
rt_invgamma <- cogmod_invgamma

#' @rdname cogmod-deprecated
#' @export
rt_invgamma_lpdf_expose <- cogmod_invgamma_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_invgamma_stanvars <- cogmod_invgamma_stanvars

# cogmod_weibull --------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_weibull <- dcogmod_weibull

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_weibull <- log_lik_cogmod_weibull

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_weibull <- posterior_epred_cogmod_weibull

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_weibull <- posterior_predict_cogmod_weibull

#' @rdname cogmod-deprecated
#' @export
rrt_weibull <- rcogmod_weibull

#' @rdname cogmod-deprecated
#' @export
rt_weibull <- cogmod_weibull

#' @rdname cogmod-deprecated
#' @export
rt_weibull_lpdf_expose <- cogmod_weibull_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_weibull_stanvars <- cogmod_weibull_stanvars

# cogmod_invweibull -----------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_invweibull <- dcogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_invweibull <- log_lik_cogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_invweibull <- posterior_epred_cogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_invweibull <- posterior_predict_cogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
rrt_invweibull <- rcogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
rt_invweibull <- cogmod_invweibull

#' @rdname cogmod-deprecated
#' @export
rt_invweibull_lpdf_expose <- cogmod_invweibull_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_invweibull_stanvars <- cogmod_invweibull_stanvars

# cogmod_logweibull -----------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_logweibull <- dcogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_logweibull <- log_lik_cogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_logweibull <- posterior_epred_cogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_logweibull <- posterior_predict_cogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
rrt_logweibull <- rcogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
rt_logweibull <- cogmod_logweibull

#' @rdname cogmod-deprecated
#' @export
rt_logweibull_lpdf_expose <- cogmod_logweibull_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_logweibull_stanvars <- cogmod_logweibull_stanvars

# cogmod_exgaussian -----------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_exgaussian <- dcogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_exgaussian <- log_lik_cogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_exgaussian <- posterior_epred_cogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_exgaussian <- posterior_predict_cogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
rrt_exgaussian <- rcogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
rt_exgaussian <- cogmod_exgaussian

#' @rdname cogmod-deprecated
#' @export
rt_exgaussian_lpdf_expose <- cogmod_exgaussian_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_exgaussian_stanvars <- cogmod_exgaussian_stanvars

# cogmod_lba1 -----------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drt_lba <- dcogmod_lba1

#' @rdname cogmod-deprecated
#' @export
log_lik_rt_lba <- log_lik_cogmod_lba1

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rt_lba <- posterior_epred_cogmod_lba1

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rt_lba <- posterior_predict_cogmod_lba1

#' @rdname cogmod-deprecated
#' @export
rrt_lba <- rcogmod_lba1

#' @rdname cogmod-deprecated
#' @export
rt_lba <- cogmod_lba1

#' @rdname cogmod-deprecated
#' @export
rt_lba_lpdf_expose <- cogmod_lba1_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rt_lba_stanvars <- cogmod_lba1_stanvars

# cogmod_lba2 -----------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
dlba <- dcogmod_lba2

#' @rdname cogmod-deprecated
#' @export
lba <- cogmod_lba2

#' @rdname cogmod-deprecated
#' @export
lba_lpdf_expose <- cogmod_lba2_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
lba_stanvars <- cogmod_lba2_stanvars

#' @rdname cogmod-deprecated
#' @export
log_lik_lba <- log_lik_cogmod_lba2

#' @rdname cogmod-deprecated
#' @export
posterior_epred_lba <- posterior_epred_cogmod_lba2

#' @rdname cogmod-deprecated
#' @export
posterior_predict_lba <- posterior_predict_cogmod_lba2

#' @rdname cogmod-deprecated
#' @export
rlba <- rcogmod_lba2

# cogmod_lnr ------------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
dlnr <- dcogmod_lnr

#' @rdname cogmod-deprecated
#' @export
lnr <- cogmod_lnr

#' @rdname cogmod-deprecated
#' @export
lnr_lpdf_expose <- cogmod_lnr_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
lnr_stanvars <- cogmod_lnr_stanvars

#' @rdname cogmod-deprecated
#' @export
log_lik_lnr <- log_lik_cogmod_lnr

#' @rdname cogmod-deprecated
#' @export
posterior_epred_lnr <- posterior_epred_cogmod_lnr

#' @rdname cogmod-deprecated
#' @export
posterior_predict_lnr <- posterior_predict_cogmod_lnr

#' @rdname cogmod-deprecated
#' @export
rlnr <- rcogmod_lnr

# cogmod_rdm ------------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
drdm <- dcogmod_rdm

#' @rdname cogmod-deprecated
#' @export
log_lik_rdm <- log_lik_cogmod_rdm

#' @rdname cogmod-deprecated
#' @export
posterior_epred_rdm <- posterior_epred_cogmod_rdm

#' @rdname cogmod-deprecated
#' @export
posterior_predict_rdm <- posterior_predict_cogmod_rdm

#' @rdname cogmod-deprecated
#' @export
prdm <- pcogmod_rdm

#' @rdname cogmod-deprecated
#' @export
rdm <- cogmod_rdm

#' @rdname cogmod-deprecated
#' @export
rdm_lpdf_expose <- cogmod_rdm_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
rdm_stanvars <- cogmod_rdm_stanvars

#' @rdname cogmod-deprecated
#' @export
rrdm <- rcogmod_rdm

# cogmod_ddm ------------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
dddm <- dcogmod_ddm

#' @rdname cogmod-deprecated
#' @export
ddm <- cogmod_ddm

#' @rdname cogmod-deprecated
#' @export
ddm_lpdf_expose <- cogmod_ddm_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
ddm_stanvars <- cogmod_ddm_stanvars

#' @rdname cogmod-deprecated
#' @export
log_lik_ddm <- log_lik_cogmod_ddm

#' @rdname cogmod-deprecated
#' @export
posterior_epred_ddm <- posterior_epred_cogmod_ddm

#' @rdname cogmod-deprecated
#' @export
posterior_predict_ddm <- posterior_predict_cogmod_ddm

#' @rdname cogmod-deprecated
#' @export
rddm <- rcogmod_ddm

# cogmod_choco ----------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
choco <- cogmod_choco

#' @rdname cogmod-deprecated
#' @export
choco_lpdf_expose <- cogmod_choco_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
choco_stanvars <- cogmod_choco_stanvars

#' @rdname cogmod-deprecated
#' @export
dchoco <- dcogmod_choco

#' @rdname cogmod-deprecated
#' @export
log_lik_choco <- log_lik_cogmod_choco

#' @rdname cogmod-deprecated
#' @export
posterior_epred_choco <- posterior_epred_cogmod_choco

#' @rdname cogmod-deprecated
#' @export
posterior_predict_choco <- posterior_predict_cogmod_choco

#' @rdname cogmod-deprecated
#' @export
rchoco <- rcogmod_choco

# cogmod_betagate -------------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
betagate <- cogmod_betagate

#' @rdname cogmod-deprecated
#' @export
betagate_lpdf_expose <- cogmod_betagate_lpdf_expose

#' @rdname cogmod-deprecated
#' @export
betagate_stanvars <- cogmod_betagate_stanvars

#' @rdname cogmod-deprecated
#' @export
dbetagate <- dcogmod_betagate

#' @rdname cogmod-deprecated
#' @export
log_lik_betagate <- log_lik_cogmod_betagate

#' @rdname cogmod-deprecated
#' @export
posterior_epred_betagate <- posterior_epred_cogmod_betagate

#' @rdname cogmod-deprecated
#' @export
posterior_predict_betagate <- posterior_predict_cogmod_betagate

#' @rdname cogmod-deprecated
#' @export
rbetagate <- rcogmod_betagate

# cogmod_betadiscrete ---------------------------------------------------

#' @rdname cogmod-deprecated
#' @export
betadiscrete <- cogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
betadiscrete_lpmf_expose <- cogmod_betadiscrete_lpmf_expose

#' @rdname cogmod-deprecated
#' @export
betadiscrete_stanvars <- cogmod_betadiscrete_stanvars

#' @rdname cogmod-deprecated
#' @export
dbetadiscrete <- dcogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
log_lik_betadiscrete <- log_lik_cogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
pbetadiscrete <- pcogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
posterior_epred_betadiscrete <- posterior_epred_cogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
posterior_predict_betadiscrete <- posterior_predict_cogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
qbetadiscrete <- qcogmod_betadiscrete

#' @rdname cogmod-deprecated
#' @export
rbetadiscrete <- rcogmod_betadiscrete
