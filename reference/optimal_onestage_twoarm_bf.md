# Internal calibration routine for one-stage two-arm BF designs

Internal calibration routine for one-stage two-arm BF designs

## Usage

``` r
optimal_onestage_twoarm_bf(
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
  progress = FALSE
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

  Numeric scalar greater than 1. Threshold used for CE(H0).

- test:

  Character string, one of `"BF01"`, `"BF+0"`, `"BF-0"`, or `"BF+-"`.

- a_0_d, b_0_d, a_0_a, b_0_a:

  Shape parameters for design and analysis priors under \\H_0\\.

- a_1_d, b_1_d, a_2_d, b_2_d:

  Shape parameters for design priors under \\H_1\\ or \\H\_+\\.

- a_1_a, b_1_a, a_2_a, b_2_a:

  Shape parameters for analysis priors under \\H_1\\ or \\H\_+\\.

- a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus:

  Optional design priors under \\H\_-\\.

- a_1_a_Hminus, b_1_a_Hminus, a_2_a_Hminus, b_2_a_Hminus:

  Optional analysis priors under \\H\_-\\.

- alloc1, alloc2:

  Fixed randomisation probabilities for arm 1 and arm 2.

- calibration:

  Character string specifying the calibration mode.

- target_power, target_type1, target_ce_h0, target_freq_power,
  target_freq_type1:

  Target operating characteristics.

- p1_grid, p2_grid:

  Grids for supremum frequentist type-I error.

- p1_power, p2_power:

  Optional true proportions for frequentist power.

- power_cushion:

  Non-negative numeric scalar applied to power targets.

- sustain_n:

  Non-negative integer. Rolling feasibility window size.

- progress:

  Logical; if `TRUE`, emit progress information.

## Value

A list with feasibility, selected design, operating characteristics, and
full search results.
