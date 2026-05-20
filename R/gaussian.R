#' Random R(t) series from an exponentiated-quadratic Gaussian process
#'
#' @param n      Number of time points in the returned series.
#' @param mu     Median of R(t).
#' @param spread Standard deviation of `log(R(t)/mu)`. Approximates the
#'    coefficient of variation of R(t) when small.
#' @param scale  Inverse length scale: smaller `scale` gives smoother R(t),
#'    larger `scale` gives a rougher trajectory.
#'
#' @return Numeric vector of length `n` containing R(t).
#' @export
random_Rt_expq = function(n, mu, spread, scale)
{
    X = seq(0, by = scale, length.out = n)
    S = exp(-0.5 * outer(X, X, "-")^2)
    L = chol(S + diag(1e-10, n))    # jitter for positive definiteness
    y = crossprod(L, stats::rnorm(n))
    return (mu * exp(spread * as.vector(y)))
}
