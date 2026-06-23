# Design or evaluate a one-stage two-arm Bayes factor trial

Calibrates or evaluates a one-stage two-arm Bayes factor design for a
binary endpoint with fixed randomisation between the two arms.

## Usage

``` r
design_twoarm_onestage_bf(
  n_min,
  n_max,
  k = 1/3,
  k_f = 3,
  test = c("BF01", "BF+0", "BF-0", "BF+-"),
  a_0_d = 1,
  b_0_d = 1,
  a_0_a = 1,
  b_0_a = 1,
  a_1_d = 1,
  b_1_d = 1,
  a_2_d = 1,
  b_2_d = 1,
  a_1_a = 1,
  b_1_a = 1,
  a_2_a = 1,
  b_2_a = 1,
  a_1_d_Hminus = 1,
  b_1_d_Hminus = 1,
  a_2_d_Hminus = 1,
  b_2_d_Hminus = 1,
  a_1_a_Hminus = 1,
  b_1_a_Hminus = 1,
  a_2_a_Hminus = 1,
  b_2_a_Hminus = 1,
  alloc1 = 0.5,
  alloc2 = 0.5,
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  target_power = 0.8,
  target_type1 = 0.05,
  target_ce_h0 = 0,
  target_freq_power = 0.8,
  target_freq_type1 = 0.05,
  p1_grid = seq(0.01, 0.99, 0.02),
  p2_grid = seq(0.01, 0.99, 0.02),
  p1_power = NULL,
  p2_power = NULL,
  power_cushion = 0,
  sustain_n = 10L,
  report_freq_type1 = FALSE,
  algorithm = c("optimal", "manual"),
  n_total = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- n_min:

  Integer. Minimum admissible total sample size.

- n_max:

  Integer. Maximum admissible total sample size.

- k:

  Numeric scalar greater than 0. Evidence threshold used for power and
  type-I error.

- k_f:

  Numeric scalar greater than 1. Threshold used for CE(H0) / PCE(H0).

- test:

  Character string, one of `"BF01"`, `"BF+0"`, `"BF-0"`, or `"BF+-"`.

- a_0_d, b_0_d, a_0_a, b_0_a:

  Shape parameters for design and analysis priors under \\H_0\\.

- a_1_d, b_1_d, a_2_d, b_2_d:

  Shape parameters for design priors under \\H_1\\ or \\H\_+\\.

- a_1_a, b_1_a, a_2_a, b_2_a:

  Shape parameters for analysis priors under \\H_1\\ or \\H\_+\\.

- a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus:

  Optional design priors under \\H\_-\\ for directional tests.

- a_1_a_Hminus, b_1_a_Hminus, a_2_a_Hminus, b_2_a_Hminus:

  Optional analysis priors under \\H\_-\\ for directional tests.

- alloc1, alloc2:

  Fixed randomisation probabilities for arm 1 and arm 2. Must be
  positive and sum to 1.

- calibration:

  Character string specifying the calibration mode. One of `"Bayesian"`,
  `"frequentist"`, `"hybrid"`, or `"full"`.

- target_power:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian power.

- target_type1:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian type-I error.

- target_ce_h0:

  Numeric scalar in \\\[0,1)\\. Optional lower bound on the corrected
  Bayesian probability of compelling evidence in favour of \\H_0\\ (or
  \\H\_-\\ for `test = "BF+-"`).

- target_freq_power:

  Numeric scalar in \\(0,1)\\. Target frequentist power under
  `p1_power, p2_power`.

- target_freq_type1:

  Numeric scalar in \\(0,1)\\. Target frequentist type-I error.

- p1_grid, p2_grid:

  Grids of true proportions used to compute supremum frequentist type-I
  error.

- p1_power, p2_power:

  Optional true proportions used for frequentist power.

- power_cushion:

  Non-negative numeric scalar. Optional additive cushion applied to the
  power targets during calibration.

- sustain_n:

  Non-negative integer. A candidate total sample size is considered
  feasible only if the relevant target constraints hold at that total
  sample size and for the next `sustain_n` larger total sample sizes in
  the search range.

- report_freq_type1:

  Logical. If `TRUE`, compute and report the frequentist type-I error
  for the final selected design even when the chosen calibration mode
  does not use frequentist criteria. This additional computation has no
  effect on the calibration itself. Defaults to `FALSE`.

- algorithm:

  Character string specifying whether the design should be optimized or
  only evaluated.

- n_total:

  Optional integer total sample size used when `algorithm = "manual"`.

- progress:

  Logical; if `TRUE`, print simple progress information during
  optimization.

- ...:

  Reserved for future extensions.

## Value

An object of class `"twoarm_onestage_bf_design"`.

## Details

The design uses one of the Bayes factor tests implemented in
[`powertwoarmbinbf01()`](https://rikokelter.github.io/bfbin2arm/reference/powertwoarmbinbf01.md).
Small values of the relevant inverted Bayes factor indicate evidence
against the null, so efficacy is concluded when the Bayes factor is
below `k`. Large values indicate evidence in favour of the null (or
\\H\_-\\ for `test = "BF+-"`), and the optional CE(H0) / PCE(H0)
constraint is evaluated using `k_f`.
