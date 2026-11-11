# Exact probability of a single-arm Bayes factor threshold event

Exact probability of a single-arm Bayes factor threshold event

## Usage

``` r
pbinbf01_singlearm_exact(
  k,
  n,
  p0,
  a0 = 1,
  b0 = 1,
  a1 = 1,
  b1 = 1,
  dp = NA_real_,
  da0 = 1,
  db0 = 1,
  da1 = 1,
  db1 = 1,
  type = c("point", "direction"),
  lower.tail = TRUE
)
```

## Details

With lower.tail = TRUE, calculates P(BF01 \< k). With lower.tail =
FALSE, calculates P(BF01 \>= k).

The argument k is a threshold on BF01. To require BF10 \> K, supply k =
1 / K.
