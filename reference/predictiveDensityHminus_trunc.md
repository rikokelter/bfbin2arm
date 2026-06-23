# Predictive density under H-: p2 \<= p1 (truncated prior)

Predictive density under H-: p2 \<= p1 (truncated prior)

## Usage

``` r
predictiveDensityHminus_trunc(
  y1,
  y2,
  n1,
  n2,
  a_1_d = 1,
  b_1_d = 1,
  a_2_d = 1,
  b_2_d = 1
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

- a_1_d, b_1_d:

  Design-prior parameters for p1.

- a_2_d, b_2_d:

  Design-prior parameters for p2.

## Value

Numeric scalar, predictive density under H-.
