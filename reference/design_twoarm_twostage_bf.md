# Design an optimal two-stage two-arm Bayes factor trial

Calibrates a two-stage two-arm Bayes factor design for a binary endpoint
by calling
[`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md)
and packaging the result in a user-friendly object of class
`"twoarm_twostage_bf_design"`.

## Usage

``` r
design_twoarm_twostage_bf(
  n1_min,
  n2_max,
  alloc1 = 0.5,
  alloc2 = 0.5,
  power_cushion = 0,
  interim_fraction = c(0.25, 0.75),
  grid_step = 1L,
  coarse_step = 4L,
  max_iter = 40L,
  ncores = getOption("bfbin2arm.ncores", 1L),
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
  calibration = c("Bayesian", "frequentist", "hybrid"),
  calibration_en = c("Bayesian", "frequentist"),
  target_power = 0.8,
  target_type1 = 0.05,
  target_ce_h0 = 0,
  target_freq_power = 0.8,
  target_freq_type1 = 0.05,
  p1_power = NULL,
  p2_power = NULL,
  p1_en_h0 = NULL,
  p2_en_h0 = NULL,
  p_null_grid = NULL,
  progress = FALSE,
  ...
)
```

## Arguments

- n1_min:

  Numeric vector of length 2, minimum interim sample sizes for arms 1
  and 2.

- n2_max:

  Numeric vector of length 2, maximum final sample sizes for arms 1 and
  2.

- alloc1, alloc2:

  Positive numbers, allocation probabilities to arms 1 and 2.

- power_cushion:

  Numeric scalar, optional extra power cushion used in the fixed-sample
  search of step 1.

- interim_fraction:

  Numeric vector of length 2 giving lower and upper bounds for the
  interim sample size in each arm as a fraction of the fixed sample
  size.

- grid_step:

  Positive integer giving the spacing of the interim design grid.

- coarse_step:

  Positive integer giving the spacing of the coarse fixed-sample search
  grid in step 1.

- max_iter:

  Integer, maximum number of total fixed-sample sizes searched in step
  1.

- ncores:

  Integer; number of parallel worker processes to use in the
  calibration. Defaults to `getOption("bfbin2arm.ncores", 1L)`. In
  vignettes and examples, a conservative value (e.g. 1 or 2) is
  recommended for CRAN checks, whereas users can increase this to
  exploit all available cores on their own machines.

- k:

  Numeric scalar, efficacy threshold; evidence against the null
  hypothesis is declared when the corresponding Bayes factor is smaller
  than `k`.

- k_f:

  Numeric scalar, futility threshold; compelling evidence for the null
  hypothesis is declared when the corresponding Bayes factor is at least
  `k_f`.

- test:

  Character string, one of `"BF01"`, `"BF+0"`, `"BF-0"`, `"BF+-"`.

- a_0_d, b_0_d, a_0_a, b_0_a:

  Shape parameters for design and analysis priors under \\H_0\\.

- a_1_d, b_1_d, a_2_d, b_2_d:

  Shape parameters for design priors under \\H_1\\ or \\H\_+\\.

- a_1_a, b_1_a, a_2_a, b_2_a:

  Shape parameters for analysis priors under \\H_1\\ or \\H\_+\\.

- a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus:

  Optional design priors under \\H\_-\\ for directional tests.

- a_1_a_Hminus, b_1_a_Hminus:

  Shape parameters of the analysis prior under the directional null
  hypothesis H0- for arm 1.

- a_2_a_Hminus, b_2_a_Hminus:

  Shape parameters of the analysis prior under the directional null
  hypothesis H0- for arm 2.

- calibration:

  Character string specifying the calibration mode at the wrapper level.
  One of `"Bayesian"`, `"frequentist"`, or `"hybrid"`. This is passed to
  [`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md)
  as `calibration_mode`.

- calibration_en:

  Character string or `NULL` specifying whether the design is ranked by
  Bayesian or frequentist expected sample size under the null
  hypothesis. This is passed to
  [`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md)
  as `calibration_EN`.

- target_power, target_type1, target_ce_h0, target_freq_power,
  target_freq_type1:

  Numeric targets for Bayesian and frequentist operating
  characteristics. These are translated to the `alpha`, `beta`,
  `alpha_freq`, and `beta_freq` arguments of
  [`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md).

- p1_power, p2_power:

  Optional true response probabilities used for frequentist power.
  Passed through to
  [`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md).

- p1_en_h0, p2_en_h0:

  Optional null response probabilities used when
  `calibration_en = "frequentist"` to compute expected sample size under
  the null.

- p_null_grid:

  Optional grid of null response probabilities used for frequentist
  type-I-error maximisation.

- progress:

  Logical; if `TRUE`, print simple progress information during the
  calibration.

- ...:

  Reserved for future extensions; currently ignored.

## Value

An object of class `"twoarm_twostage_bf_design"`.

## Details

The design uses one of the Bayes factor tests implemented in
[`powertwoarmbinbf01()`](https://rikokelter.github.io/bfbin2arm/reference/powertwoarmbinbf01.md).
Small values of the relevant inverted Bayes factor indicate evidence
against the null, so efficacy is concluded when the Bayes factor is
below `k`. Large values indicate evidence in favour of the null (or
\\H\_-\\ for `test = "BF+-"`), and the optional CE(H0) / PCE(H0)
constraint is evaluated using `k_f`.
