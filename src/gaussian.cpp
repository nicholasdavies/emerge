#include <cpp11.hpp>
#include "random.h"

// Random R(t) series from a Matérn-3/2 Gaussian process on a unit-spaced
// grid, mapped to log-normal: R(t) = mu * exp(sigma*y(t) - sigma^2/2).
//    n       number of data points
//    mu      expectation of R(t) over all realisations
//    sigma   standard deviation of R(t) on log scale
//    ell     length constant for Matérn process
//    seed    random seed
[[cpp11::register]]
cpp11::doubles random_Rt_m32(int n, double mu, double sigma, double ell, unsigned int seed)
{
    Engine e(seed);
    std::normal_distribution<double> randn(0.0, 1.0);

    double lam  = std::sqrt(3.0) / ell;
    double ex   = std::exp(-lam);
    double ex2  = ex * ex;
    double lam2 = lam * lam;
    double lam3 = lam2 * lam;

    // Transition matrix A
    double a11 = ex * (1.0 + lam);
    double a12 = ex;
    double a21 = -ex * lam2;
    double a22 = ex * (1.0 - lam);

    // Process-noise covariance Q = P_inf - A P_inf A^T
    double q11 = 1.0  - ex2 * (1.0 + 2.0 * lam + 2.0 * lam2);
    double q12 = 2.0 * ex2 * lam3;
    double q22 = lam2 * (1.0 - ex2 * (2.0 * lam2 - 2.0 * lam + 1.0));

    // Cholesky of Q (lower triangular)
    double L11 = std::sqrt(q11);
    double L21 = q12 / L11;
    double L22 = std::sqrt(q22 - L21 * L21);

    // Bias correction so that E[R(t)] = mu
    double bias = 0.5 * sigma * sigma;

    auto R = cpp11::writable::doubles(n);

    // Initial state from stationary distribution P_inf = diag(1, lam^2)
    double y  = randn(e);
    double yp = lam * randn(e);
    R[0] = mu * std::exp(sigma * y - bias);

    // Generate R[t]
    for (int t = 1; t < n; ++t)
    {
        double n1 = randn(e);
        double n2 = randn(e);
        double w1 = L11 * n1;
        double w2 = L21 * n1 + L22 * n2;

        double y_new  = a11 * y  + a12 * yp + w1;
        double yp_new = a21 * y  + a22 * yp + w2;
        y  = y_new;
        yp = yp_new;

        R[t] = mu * std::exp(sigma * y - bias);
    }

    return R;
}
