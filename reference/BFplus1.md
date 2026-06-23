# Bayes factor BF+1: H+ vs H1

Bayes factor BF+1: H+ vs H1

## Usage

``` r
BFplus1(y1, y2, n1, n2, a_1_a = 1, b_1_a = 1, a_2_a = 1, b_2_a = 1)
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

- a_1_a, b_1_a:

  Shape parameters of the Beta prior for \\p_1\\ under the alternative
  (analysis prior).

- a_2_a, b_2_a:

  Shape parameters of the Beta prior for \\p_2\\ under the alternative
  (analysis prior).

## Value

Numeric scalar, BF\\\_{+1}\\ = posterior odds / prior odds for H+ vs H1.
