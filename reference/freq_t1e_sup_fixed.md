# Fixed-sample frequentist type-I error supremum over a null grid

Internal helper. For a given fixed-sample design (n1, n2), computes the
supremum of the frequentist type-I error over a grid of null parameter
configurations, using `freq_oc_twoarm_fixed()`.

## Usage

``` r
freq_t1e_sup_fixed(
  n1,
  n2,
  k,
  k_f,
  test,
  p_null_grid = NULL,
  a_0_a,
  b_0_a,
  a_1_a,
  b_1_a,
  a_2_a,
  b_2_a,
  a_1_a_Hminus,
  b_1_a_Hminus,
  a_2_a_Hminus,
  b_2_a_Hminus,
  alpha_target = NULL,
  tol_excess = 1e-04
)
```
