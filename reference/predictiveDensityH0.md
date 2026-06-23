# Predictive density under H0: p1 = p2 = p

Beta-binomial predictive density for data (y1,y2) under H0.

## Usage

``` r
predictiveDensityH0(y1, y2, n1, n2, a_0_d = 1, b_0_d = 1)
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

- a_0_d, b_0_d:

  Design-prior parameters for common p under H0.

## Value

Numeric scalar, predictive density.
