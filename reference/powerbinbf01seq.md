# Bayesian and frequentist operating characteristics for a single-arm two-stage BF design

Computes naive fixed-sample and corrected two-stage operating
characteristics for a single-arm binomial design with one interim
analysis for futility. The Bayes factor is oriented as `BF01`, so
efficacy corresponds to small values (`BF01 <= k`) and futility
corresponds to large values (`BF01 >= kf`).

## Usage

``` r
powerbinbf01seq(
  n1,
  n2,
  k,
  kf,
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

- n1:

  Integer scalar. Interim sample size.

- n2:

  Integer scalar. Final sample size, with `n1 < n2`.

- k:

  Numeric scalar. Efficacy threshold on the `BF01` scale.

- kf:

  Numeric scalar. Futility threshold on the `BF01` scale.

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

A list with Bayesian and frequentist operating characteristics.

## Details

Bayesian operating characteristics are computed under separate design
priors:

- for `type = "direction"`, Bayesian power averages over `p > p0` under
  the H1 design prior truncated to `(p0, 1]`, and Bayesian type-I error
  averages over `p <= p0` under the H0 design prior truncated to
  `[0, p0]`;

- for `type = "point"`, Bayesian power averages under the H1 design
  prior on `(0, 1)`, and Bayesian type-I error is evaluated at the point
  null `p = p0`.

If `dp` is supplied, additional frequentist power under `H1` is computed
at the fixed point alternative `p = dp`. Frequentist type-I error is
computed at `p = p0`.
