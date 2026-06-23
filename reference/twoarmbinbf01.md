# Two-arm binomial Bayes factor BF01

Computes the Bayes factor BF\\\_{01}\\ comparing the point-null \\H_0:
p_1 = p_2\\ to the alternative \\H_1: p_1 \neq p_2\\ in a two-arm
binomial setting with Beta priors.

## Usage

``` r
twoarmbinbf01(
  y1,
  y2,
  n1,
  n2,
  a_0_a = 1,
  b_0_a = 1,
  a_1_a = 1,
  b_1_a = 1,
  a_2_a = 1,
  b_2_a = 1
)
```

## Arguments

- y1:

  Number of successes in arm 1 (control).

- y2:

  Number of successes in arm 2 (treatment).

- n1:

  Sample size in arm 1.

- n2:

  Sample size in arm 2.

- a_0_a, b_0_a:

  Shape parameters of the Beta prior for the common-\\p\\ under the null
  model (analysis prior).

- a_1_a, b_1_a:

  Shape parameters of the Beta prior for \\p_1\\ under the alternative
  (analysis prior).

- a_2_a, b_2_a:

  Shape parameters of the Beta prior for \\p_2\\ under the alternative
  (analysis prior).

## Value

Numeric scalar, the Bayes factor BF\\\_{01}\\.

## Examples

``` r
twoarmbinbf01(10, 20, 30, 30)
#> [1] 0.1202496
```
