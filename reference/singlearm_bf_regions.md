# Construct strict Bayes factor evidence regions

Classifies Bayes factors into an efficacy region and, optionally, a
region of compelling evidence for H0 and an inconclusive region.

## Usage

``` r
singlearm_bf_regions(bfvals, k, k_ce = NULL)
```

## Arguments

- bfvals:

  Numeric vector of nonnegative Bayes factors BF01 comparing H0 with H1,
  without missing values. Zero and positive infinity are permitted.

- k:

  Positive finite efficacy threshold on BF01. Efficacy is defined by the
  strict inequality \\BF\_{01}\<k\\.

- k_ce:

  Optional positive finite threshold greater than `k`. Compelling
  evidence for H0 is defined by the strict inequality
  \\BF\_{01}\>k\_{\mathrm{ce}}\\. If `NULL`, only efficacy and its
  complement are returned.

## Value

A named list of logical vectors, each with the same length as `bfvals`.
If `k_ce = NULL`, the components are:

- efficacy:

  Whether \\BF\_{01}\<k\\.

- not_efficacy:

  Whether \\BF\_{01}\ge k\\.

Otherwise, the components are:

- efficacy:

  Whether \\BF\_{01}\<k\\.

- compelling_h0:

  Whether \\BF\_{01}\>k\_{\mathrm{ce}}\\.

- inconclusive:

  Whether \\k \le BF\_{01} \le k\_{\mathrm{ce}}\\.

## Details

When `k_ce` is supplied, the regions are mutually exclusive and
exhaustive. Values equal to either threshold are classified as
inconclusive: \\k \le BF\_{01} \le k\_{\mathrm{ce}}\\.

When `k_ce = NULL`, `not_efficacy` denotes the complement of efficacy,
\\BF\_{01}\ge k\\. It should not be interpreted as compelling evidence
for H0.

The usual evidential interpretation uses \\0\<k\<1\<k\_{\mathrm{ce}}\\,
although the function only requires `k > 0` and, when supplied,
`k_ce > k`.

Comparisons are applied directly to the supplied floating-point values
without an explicit numerical tolerance.
