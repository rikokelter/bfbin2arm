# Internal calibration routine for one-stage single-arm BF designs

Internal calibration routine for one-stage single-arm BF designs

## Usage

``` r
optimal_onestage_singlearm_bf(
  n_min,
  n_max,
  k,
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
  power_cushion = 0,
  k_ce = NULL,
  sustain_n = 10L
)
```

## Arguments

- n_min:

  Integer. Minimum admissible sample size in the search grid.

- n_max:

  Integer. Maximum admissible sample size in the search grid.

- k:

  Numeric scalar greater than 0. Evidence threshold on the \\BF\_{01}\\
  scale used for efficacy.

- p0:

  Numeric scalar in \\(0,1)\\. Null response probability.

- a0, b0:

  Positive numeric scalars. Beta analysis-prior parameters under
  \\H_0\\.

- a1, b1:

  Positive numeric scalars. Beta analysis-prior parameters under
  \\H_1\\.

- dp:

  Optional numeric scalar in \\(0,1)\\. Fixed point alternative used for
  frequentist power calculations under \\H_1\\.

- da0, db0:

  Positive numeric scalars. Beta design-prior parameters under \\H_0\\.

- da1, db1:

  Positive numeric scalars. Beta design-prior parameters under \\H_1\\.

- type:

  Character string specifying the Bayes-factor test. One of `"point"` or
  `"direction"`.

- calibration:

  Character string specifying the calibration mode. One of `"Bayesian"`,
  `"frequentist"`, `"hybrid"`, or `"full"`.

- target_power:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian power.

- target_type1:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian type-I error.

- target_ce_h0:

  Numeric scalar in \\\[0,1)\\. Optional lower bound on the corrected
  Bayesian probability of compelling evidence in favour of \\H_0\\.

- target_freq_power:

  Numeric scalar in \\(0,1)\\. Target corrected frequentist power at
  `dp`.

- target_freq_type1:

  Numeric scalar in \\(0,1)\\. Target corrected frequentist type-I error
  at \\p = p_0\\.

- power_cushion:

  Non-negative numeric scalar. Optional additive cushion applied to the
  power targets during calibration.

- k_ce:

  Optional numeric scalar greater than 1. Threshold on the \\BF\_{01}\\
  scale used for CE(H0) / PCE(H0) calculations.

- sustain_n:

  Non-negative integer. A candidate sample size is declared feasible
  only if the relevant constraints are satisfied at that sample size and
  for the next `sustain_n` larger sample sizes, subject to the search
  range.

## Value

A list with feasibility, selected design, operating characteristics, and
full search results.
