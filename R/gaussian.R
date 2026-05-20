#' Get a random Rt series
#' @export
random_Rt_expq = function(n, mu, spread, scale)
{
    X = seq(0, by = scale, length.out = n)
    S = exp(-0.5 * outer(X, X, "-")^2)
    L = chol(S + diag(1e-10, n))    # jitter for numerical PD
    y = crossprod(L, stats::rnorm(n))
    return (mu * exp(spread * as.vector(y)))
}
