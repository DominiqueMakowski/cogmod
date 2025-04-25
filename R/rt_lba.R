#' Single-Accumulator LBA model
#'
#' @description
#' This function simulates reaction times using a single-accumulator
#' version of the Linear Ballistic Accumulator (LBA) model.
#'
#' @param n A single positive integer giving the number of trials.
#' @param drift Mean drift rate.
#' @param sigma Standard deviation of the drift rate.
#' @param sigmabias The starting-point range (A); must be positive.
#' @param bs The additional offset such that b = sigmabias + bs; must be positive.
#' @param ndt Non-decision time.
#' 
#' @examples
#' # Simulate 1000 trials with specified parameters
#' rts <- rrt_lba(n = 1000, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3)
#' hist(rts, breaks = 100, xlab = "RT (s)", col = "lightblue")
#' 
#' @export
rrt_lba <- function(n, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3) {

  # Input checks
  if (length(n)!=1 || n<=0 || n!=floor(n)) {
    stop("`n` must be a single positive integer.")
  }

  # Validate and recycle parameters
  p <- .prepare_lba(n = n, drift = drift, sigma = sigma, sigmabias = sigmabias, bs = bs, ndt = ndt)

  # Derived threshold
  b <- p$sigmabias + p$bs

  # Sample, for each draw, a drift rate from truncated Normal and a uniform start point.
  v <- mapply(function(mu, s) {
    .rnorm_truncated(1, mean = mu, sd = s, lower = 0, upper = Inf)
  }, mu = p$drift, s = p$sigma)
  
  sp <- mapply(function(a) {
    stats::runif(1, min = 0, max = a)
  }, a = p$sigmabias)
  
  # Compute decision times (RT component) for each draw:
  dt <- (b - sp) / v
  
  # Final simulated RT is non-decision time plus the decision time.
  rt <- p$ndt + dt
  rt

}



#' @rdname rrt_lba
#' @inheritParams lnr
#' @export
drt_lba <- function(x, drift = 3, sigma = 1, sigmabias = 0.5, bs = 0.5, ndt = 0.3, log = FALSE) {
  
  # Validate and recycle parameters
  p <- .prepare_lba(x = x, drift = drift, sigma = sigma, sigmabias = sigmabias, bs = bs, ndt = ndt)
  
  # Derived quantities 
  A <- p$sigmabias
  b <- p$sigmabias + p$bs
  T <- p$x - p$ndt   # decision‐time portion

  # Initialize output vector
  out <- numeric(p$ndraws)
  
  # For indices where x is not above ndt, set density (or log-density) accordingly.
  not_ok <- T <= 0
  if (any(not_ok)) {
    if (log) {
      out[not_ok] <- -Inf
    } else {
      out[not_ok] <- 0
    }
  }
  
  # Process valid cases
  ok <- T > 0L
  if (any(ok)) {
    Tp <- T[ok]
    dr <- p$drift[ok]
    sg <- p$sigma[ok]
    Ao <- A[ok]
    bo <- b[ok]
    
    st <- sg * Tp
    # avoid exactly zero
    st <- pmax(st, 1e-10)
    
    z1 <- (bo - Ao - dr * Tp) / st
    z2 <- (bo - dr * Tp) / st
    
    # defective density: drift * (phi(z2) - phi(z1))  +  sigma * (φ(z1) - φ(z2))
    n2 <- dr * (stats::pnorm(z2) - stats::pnorm(z1)) + sg * (stats::dnorm(z1) - stats::dnorm(z2))
    n2 <- pmax(n2, 1e-10)
    
    # normalization: Pr(drift > 0) = 1 - phi(-drift/sigma)
    # Compute in the log domain.
    log_ppos <- stats::pnorm(-dr/sg, log.p = TRUE)
    # Compute log(1 - phi(-drift/sigma)) stably:
    log_norm <- log1p(-exp(log_ppos))
    
    if (log) {
      out[ok] <- log(n2) - log(Ao) - log_norm
    } else {
      out[ok] <- (1 / Ao) * n2 * exp(-log_norm)
    }
  }
  
  out
}



# Internal helpers --------------------------------------------------

# Helper function to prepare (and recycle) LBA parameters.
# Either x (observed RTs) or n (number of draws) must be supplied.
.prepare_lba <- function(x = NULL, n = NULL, drift, sigma, sigmabias, bs, ndt) {

  # Validate parameters once
  if (any(sigma <= 0))     stop("`sigma` must be positive.")
  if (any(sigmabias <= 0)) stop("`sigmabias` must be positive.")
  if (any(bs <= 0))        stop("`bs` must be positive.")
  if (any(ndt < 0))        stop("`ndt` must be non-negative.")


  if (!is.null(n)) {
    n_out <- as.integer(n)
    out <- list(ndraws = n_out)
  } else if (!is.null(x)) {
    n_out <- n_out <- max(length(x), length(drift))
    out <- list(x = rep(x, length.out = n_out), ndraws = n_out)
  }

  out$drift <- rep(drift, length.out = n_out)
  out$sigma <- rep(sigma, length.out = n_out)
  out$sigmabias <- rep(sigmabias, length.out = n_out)
  out$bs <- rep(bs, length.out = n_out)
  out$ndt <- rep(ndt, length.out = n_out)

  out
}
