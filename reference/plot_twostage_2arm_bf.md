# Plot an optimal two-stage two-arm Bayes factor design

Given the result from
[`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md),
this function produces a six-panel base R plot showing the design
schematic, operating characteristics, and the design and analysis priors
under \\H_0\\ and \\H_1\\.

## Usage

``` r
plot_twostage_2arm_bf(
  res,
  main = "Optimal two-stage two-arm Bayes factor design"
)
```

## Arguments

- res:

  A list returned by
  [`optimal_twostage_2arm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_2arm_bf.md),
  containing components `$design`, `$naive_oc`, `$occ` and `$priors`.

- main:

  Character string with the main title of the plot.

## Value

Invisibly returns `NULL`; called for its side effect of producing a
plot.

## Examples

``` r
res <- optimal_twostage_2arm_bf(
  alpha = 0.10, beta = 0.20, k = 1/3, k_f = 3,
  n1_min = c(3, 3), n2_max = c(8, 8),
  alloc1 = 0.5, alloc2 = 0.5,
  power_cushion = 0,
  interim_fraction = c(0.5, 0.5),
  grid_step = 1,
  progress = FALSE,
  max_iter = 16,
  test = "BF01",
  a_0_d = 1, b_0_d = 1,
  a_0_a = 1, b_0_a = 1,
  a_1_d = 1, b_1_d = 1,
  a_2_d = 1, b_2_d = 1,
  a_1_a = 1, b_1_a = 1,
  a_2_a = 1, b_2_a = 1
)
if (is.numeric(res$design) && length(res$design) == 4 && !anyNA(res$design)) {
  plot_twostage_2arm_bf(res)
}
```
