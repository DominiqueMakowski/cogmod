% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/rt_gamma_brms.R
\name{rt_gamma}
\alias{rt_gamma}
\alias{rt_gamma_lpdf_expose}
\alias{rt_gamma_stanvars}
\alias{log_lik_rt_gamma}
\alias{posterior_predict_rt_gamma}
\alias{posterior_epred_rt_gamma}
\title{Shifted Gamma Model for brms}
\usage{
rt_gamma(
  link_mu = "softplus",
  link_sigma = "softplus",
  link_tau = "logit",
  link_minrt = "identity"
)

rt_gamma_lpdf_expose()

rt_gamma_stanvars()

log_lik_rt_gamma(i, prep)

posterior_predict_rt_gamma(i, prep, ...)

posterior_epred_rt_gamma(prep)
}
\arguments{
\item{link_mu, link_sigma, link_tau, link_minrt}{Link functions for the parameters.}

\item{i, prep}{For brms' functions to run: index of the observation and a \code{brms} preparation object.}

\item{...}{Additional arguments.}
}
\value{
A \code{brms::custom_family} object.
}
\description{
Provides the necessary functions to use a shifted Gamma distribution
as a custom family in \code{brms}. This version is parameterized using
\code{tau} (proportion of non-decision time relative to minimum RT) and \code{minrt}
(minimum possible RT), where the non-decision time \code{ndt = tau * minrt}.
The distribution is parameterized by the shape (\code{mu} in brms, alpha in standard Gamma)
and scale (\code{sigma} in brms, beta in standard Gamma),
as well as \code{tau}, and \code{minrt} (which is typically fixed).
}
