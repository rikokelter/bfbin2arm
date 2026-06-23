# Bayesian and frequentist operating characteristics for a fixed-sample single-arm BF design

Computes operating characteristics for a genuine fixed-sample single-arm
binomial design with final efficacy decision based on `BF01 <= k`.

## Usage

``` r
powerbinbf01_fixed(
  n,
  k,
  p0,
  a0 = 1,
  b0 = 1,
  a1 = 1,
  b1 = 1,
  da0 = 1,
  db0 = 1,
  da1 = 1,
  db1 = 1,
  dp = NA_real_,
  type = c("point", "direction"),
  k_ce = NULL,
  grid_size = 801L
)
```

## Arguments

- n:

  Integer scalar. Total sample size.

- k:

  Numeric scalar. Efficacy threshold on the `BF01` scale.

- p0:

  Numeric scalar in `(0, 1)`. Null response probability.

- a0, b0:

  Numeric scalars. Beta analysis-prior parameters under H0.

- a1, b1:

  Numeric scalars. Beta analysis-prior parameters under H1.

- da0, db0:

  Numeric scalars. Beta design-prior parameters under H0.

- da1, db1:

  Numeric scalars. Beta design-prior parameters under H1.

- dp:

  Optional numeric scalar in `(0,1)`. If supplied, frequentist power
  under `H1` is computed at `p = dp`.

- type:

  Character string. One of `\"point\"` or `\"direction\"`.

- k_ce:

  Optional numeric scalar greater than 1. Threshold for compelling
  evidence in favour of `H0` on the `BF01` scale.

- grid_size:

  Integer number of grid points used for numerical averaging.

## Value

A list with Bayesian and frequentist operating characteristics for the
fixed-sample design.

## Details

Bayesian operating characteristics are computed under separate design
priors:

- for `type = "direction"`, Bayesian power averages over `p > p0` under
  the H1 design prior truncated to `(p0, 1]`, Bayesian type-I error
  averages over `p <= p0` under the H0 design prior truncated to
  `[0, p0]`, and CE(H0) is averaged over the same truncated H0 design
  prior;

- for `type = "point"`, Bayesian power averages under the H1 design
  prior on `(0, 1)`, Bayesian type-I error is evaluated at the point
  null `p = p0`, and CE(H0) is also evaluated at `p = p0`.
