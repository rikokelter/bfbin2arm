# Optimal two-stage two-arm Bayes-factor design for binary endpoints

Computes an optimal two-stage two-arm Bayes-factor design for binary
endpoints, minimizing the expected sample size under the null hypothesis
while correcting the operating characteristics for the possibility of
early stopping for futility.

## Usage

``` r
optimal_twostage_2arm_bf(
  alpha = 0.05,
  beta = 0.2,
  k = 1/3,
  k_f = 3,
  n1_min = c(5, 5),
  n2_max = c(50, 50),
  alloc1 = 0.5,
  alloc2 = 0.5,
  power_cushion = 0,
  pceH0 = NULL,
  interim_fraction = c(0, 1),
  grid_step = 1L,
  coarse_step = 10L,
  progress = TRUE,
  max_iter = 10000L,
  ncores = getOption("bfbin2arm.ncores", 1L),
  compute_freq_oc = NULL,
  calibration_mode = c("Bayesian", "frequentist", "hybrid"),
  calibration_EN = NULL,
  p1_EN_H0 = NULL,
  p2_EN_H0 = NULL,
  alpha_freq = alpha,
  beta_freq = beta,
  p1_power = NULL,
  p2_power = NULL,
  p_null_grid = NULL,
  test = "BF01",
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
  b_2_a_Hminus = 1
)
```

## Arguments

- alpha:

  Numeric scalar, Bayesian type-I-error target.

- beta:

  Numeric scalar, 1 minus the minimal Bayesian power target.

- k:

  Numeric scalar, efficacy threshold; evidence against the null
  hypothesis is declared when the corresponding Bayes factor is smaller
  than `k`.

- k_f:

  Numeric scalar, futility threshold; compelling evidence for the null
  hypothesis is declared when the corresponding Bayes factor is at least
  `k_f`.

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

- pceH0:

  Optional numeric scalar in `[0,1]`. If specified, candidate two-stage
  designs must satisfy corrected `CE_H0 >= pceH0`.

- interim_fraction:

  Numeric vector of length 2 giving lower and upper bounds for the
  interim sample size in each arm as a fraction of the fixed sample
  size.

- grid_step:

  Positive integer giving the spacing of the interim design grid.

- coarse_step:

  Positive integer giving the spacing of the coarse fixed-sample search
  grid in step 1.

- progress:

  Logical; if `TRUE`, prints progress information.

- max_iter:

  Integer, maximum number of total fixed-sample sizes searched in step
  1.

- ncores:

  Integer; number of parallel worker processes to use in the
  calibration. Defaults to `getOption("bfbin2arm.ncores", 1L)`. In
  vignettes and examples, a conservative value (e.g. 1 or 2) is
  recommended for CRAN checks, whereas users can increase this to
  exploit all available cores on their own machines.

- compute_freq_oc:

  Logical or `NULL`. Controls whether frequentist operating
  characteristics are computed for candidate two-stage designs during
  the search.

- calibration_mode:

  Character string specifying the calibration mode. Must be one of
  `"Bayesian"`, `"frequentist"`, or `"hybrid"`.

- calibration_EN:

  Character string or `NULL` specifying whether the design is ranked by
  Bayesian or frequentist expected sample size under the null
  hypothesis.

- p1_EN_H0, p2_EN_H0:

  Numeric scalars specifying the null response probabilities in control
  and treatment arm used when `calibration_EN = "frequentist"`.

- alpha_freq:

  Numeric scalar, frequentist type-I error target.

- beta_freq:

  Numeric scalar, 1 minus the frequentist power target.

- p1_power, p2_power:

  Numeric scalars specifying the success probabilities in control and
  treatment arm used for the frequentist power calculation.

- p_null_grid:

  Optional numeric vector giving the grid of null response probabilities
  used for frequentist type-I-error maximization. If `NULL`, a default
  grid is used.

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

## Value

A list with the following components:

- design:

  Four-element integer vector containing the selected two-stage design:
  interim sample sizes in arms 1 and 2 followed by final sample sizes in
  arms 1 and 2.

- naive_oc:

  Named list of uncorrected fixed-sample operating characteristics and
  fixed-sample sizes found in step 1.

- occ:

  Named numeric vector of corrected Bayesian operating characteristics
  for the selected two-stage design.

- priors:

  List storing design hyperparameters and search settings.

- freq_occ:

  Named numeric vector with fixed-sample and two-stage frequentist
  operating characteristics for the final design when frequentist
  calibration or reporting is active; otherwise `NULL`.

- conv:

  Character string describing the search outcome. Typical values include
  `"converged"`, `"no_feasible_fixed"`, `"no_interim_grid"`, and
  `"no_feasible_design"`. In frequentist or hybrid calibration modes,
  additional informative status values may be returned when the best
  available design is returned although all requested constraints were
  not fully satisfied.

## Examples

``` r
## Fast Bayesian example with small search space
res <- optimal_twostage_2arm_bf(
  alpha = 0.10,
  beta = 0.20,
  k = 1 / 3,
  k_f = 3,
  n1_min = c(3, 3),
  n2_max = c(12, 12),
  alloc1 = 0.5,
  alloc2 = 0.5,
  power_cushion = 0,
  pceH0 = NULL,
  interim_fraction = c(0.25, 0.75),
  grid_step = 2L,
  coarse_step = 4L,
  progress = FALSE,
  max_iter = 24L,
  calibration_mode = "Bayesian",
  test = "BF01",
  a_0_d = 1, b_0_d = 1,
  a_0_a = 1, b_0_a = 1,
  a_1_d = 1, b_1_d = 1,
  a_2_d = 1, b_2_d = 1,
  a_1_a = 1, b_1_a = 1,
  a_2_a = 1, b_2_a = 1
)
res$design
#> [1] NA NA NA NA
res$occ
#> NULL

# \donttest{
res2 <- optimal_twostage_2arm_bf(
  alpha = 0.05,
  beta = 0.20,
  k = 1 / 3,
  k_f = 3,
  n1_min = c(5, 5),
  n2_max = c(20, 20),
  alloc1 = 0.5,
  alloc2 = 0.5,
  power_cushion = 0.02,
  pceH0 = 0.50,
  interim_fraction = c(0.25, 0.75),
  grid_step = 1L,
  coarse_step = 4L,
  progress = FALSE,
  max_iter = 40L,
  calibration_mode = "Bayesian",
  test = "BF+0",
  a_0_d = 1, b_0_d = 1,
  a_0_a = 1, b_0_a = 1,
  a_1_d = 1, b_1_d = 2,
  a_2_d = 2, b_2_d = 1,
  a_1_a = 1, b_1_a = 1,
  a_2_a = 1, b_2_a = 1,
  a_1_d_Hminus = 1, b_1_d_Hminus = 1,
  a_2_d_Hminus = 1, b_2_d_Hminus = 1,
  a_1_a_Hminus = 1, b_1_a_Hminus = 1,
  a_2_a_Hminus = 1, b_2_a_Hminus = 1
)
res2$design
#> [1] NA NA NA NA
res2$occ
#> NULL
# }
```
