# Optimal two-stage single-arm Bayes factor design

Searches over admissible two-stage single-arm designs with a binary
endpoint and returns the feasible design with smallest expected sample
size under `H0`.

## Usage

``` r
optimal_twostage_singlearm_bf(
  n1_min,
  n2_max,
  k,
  k_f,
  p0,
  a0 = 1,
  b0 = 1,
  a1 = 1,
  b1 = 1,
  dp = NA_real_,
  da0 = 1,
  db0 = 1,
  da1 = 1,
  db1 = 1,
  type = c("point", "direction"),
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  target_power = 0.8,
  target_type1 = 0.05,
  target_ce_h0 = 0,
  target_freq_power = 0.8,
  target_freq_type1 = 0.05,
  power_cushion = 0
)
```

## Arguments

- n1_min:

  Minimum admissible interim sample size.

- n2_max:

  Maximum admissible final sample size.

- k:

  Efficacy threshold on the `BF01` scale.

- k_f:

  Futility threshold on the `BF01` scale.

- p0:

  Null response probability.

- a0, b0:

  Beta analysis-prior parameters under `H0`.

- a1, b1:

  Beta analysis-prior parameters under `H1`.

- dp:

  Optional fixed point alternative used for frequentist power.

- da0, db0:

  Beta design-prior parameters under `H0`.

- da1, db1:

  Beta design-prior parameters under `H1`.

- type:

  Character string; one of `"point"` or `"direction"`.

- calibration:

  Character string; one of `"Bayesian"`, `"frequentist"`, `"hybrid"`, or
  `"full"`.

- target_power:

  Target corrected Bayesian power.

- target_type1:

  Target corrected Bayesian type-I error.

- target_ce_h0:

  Optional lower bound on corrected Bayesian compelling evidence in
  favour of `H0`.

- target_freq_power:

  Target corrected frequentist power at `dp`.

- target_freq_type1:

  Target corrected frequentist type-I error at `p0`.

- power_cushion:

  Optional additive cushion for the fixed-sample power target in the
  first step of the search.

## Value

A list describing the optimal design and search results.

## Details

Analysis priors are specified separately under `H0` and `H1` via
`a0, b0, a1, b1`. Design priors are specified separately under `H0` and
`H1` via `da0, db0, da1, db1`.
