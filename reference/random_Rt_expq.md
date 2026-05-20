# Random R(t) series from an exponentiated-quadratic Gaussian process

Random R(t) series from an exponentiated-quadratic Gaussian process

## Usage

``` r
random_Rt_expq(n, mu, spread, scale)
```

## Arguments

- n:

  Number of time points in the returned series.

- mu:

  Median of R(t).

- spread:

  Standard deviation of `log(R(t)/mu)`. Approximates the coefficient of
  variation of R(t) when small.

- scale:

  Inverse length scale: smaller `scale` gives smoother R(t), larger
  `scale` gives a rougher trajectory.

## Value

Numeric vector of length `n` containing R(t).
