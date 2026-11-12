# Exact threshold probabilities for a single-arm Bayes factor

Calculates the probability that the Bayes factor BF01 falls strictly
below a specified threshold, or the complementary probability that it is
greater than or equal to that threshold. Probabilities are evaluated by
summing over all possible binomial response counts.

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

## Arguments

- k:

  Positive finite threshold on BF01, the Bayes factor comparing H0 with
  H1. To calculate the probability of BF10 \> K, supply `k = 1 / K` and
  `lower.tail = TRUE`.

- n:

  Number of patients or binomial trials; a nonnegative integer.

- p0:

  Null response probability for `type = "point"`, or the boundary
  between the null and alternative parameter regions for
  `type = "direction"`. For directional hypotheses, must lie strictly
  between zero and one.

- a0, b0:

  Positive beta shape parameters for the analysis prior under H0 when
  `type = "direction"`. The prior is truncated and normalized on
  \\\[0,p_0\]\\. These parameters are not used for the point-null
  hypothesis.

- a1, b1:

  Positive beta shape parameters for the analysis prior under H1. The
  prior is supported on \\\[0,1\]\\ for `type = "point"` and truncated
  and normalized on \\(p_0,1\]\\ for `type = "direction"`.

- dp:

  Fixed response probability at which to evaluate the threshold
  probability. If `NA_real_`, probabilities are averaged over the H1
  design prior specified by `da1` and `db1`. Otherwise, binomial
  probabilities at `dp` are used.

- da0, db0:

  Retained for interface compatibility. These arguments are not used by
  this function, which evaluates either a fixed-probability distribution
  or the H1 design-prior predictive distribution.

- da1, db1:

  Positive beta shape parameters for the H1 design prior used when
  `dp = NA_real_`. This prior is supported on \\\[0,1\]\\ for
  `type = "point"` and truncated and normalized on \\(p_0,1\]\\ for
  `type = "direction"`. Ignored when `dp` is supplied.

- type:

  Character string specifying the hypotheses. `"point"` compares the
  point null \\p=p_0\\ with an alternative having a beta prior on
  \\\[0,1\]\\. `"direction"` compares \\p \le p_0\\ with \\p \> p_0\\.

- lower.tail:

  Logical. If `TRUE`, return \\P(BF\_{01}\<k)\\. If `FALSE`, return the
  complementary probability \\P(BF\_{01}\ge k)\\.

## Value

A numeric scalar giving the requested threshold probability.

## Details

Analysis priors determine the Bayes factor at each response count. The
fixed response probability `dp`, or the H1 design prior when `dp` is
missing, determines the probabilities assigned to those counts.

The lower-tail inequality is strict, unlike the usual inclusive
cumulative distribution convention. Equality is included in the
complementary upper tail. Thus, `lower.tail = FALSE` does not calculate
\\P(BF\_{01}\>k)\\.

Comparisons are applied directly to the computed floating-point Bayes
factors without an explicit numerical tolerance.
