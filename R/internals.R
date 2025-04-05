# Stable computation of log(1 - exp(x))
# A stable version of log1m_exp (assumes x is scalar or vector)
.log1m_exp <- function(x) {
  # Computes log(1 - exp(x)) in a numerically stable way.
  # x should be negative; for vectorized input, use ifelse.
  ifelse(x < log(0.5), log1p(-exp(x)), log(-expm1(x)))
}
