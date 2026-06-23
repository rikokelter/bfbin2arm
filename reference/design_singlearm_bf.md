# Design or evaluate a single-arm two-stage Bayes factor trial

Calibrates or evaluates a single-arm two-stage Bayes factor design for a
binary endpoint with one interim analysis for futility.

## Usage

``` r
design_singlearm_bf(
  n1_min,
  n2_max,
  k,
  k_f,
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
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  target_power = 0.8,
  target_type1 = 0.05,
  target_ce_h0 = 0,
  target_freq_power = 0.8,
  target_freq_type1 = 0.05,
  algorithm = c("optimal", "manual"),
  interim = NULL,
  final = NULL,
  power_cushion = 0,
  ...
)
```

## Arguments

- n1_min:

  Integer. Minimum admissible interim sample size.

- n2_max:

  Integer. Maximum admissible final sample size.

- k:

  Numeric scalar greater than 0. Efficacy threshold on the \\BF\_{01}\\
  scale.

- k_f:

  Numeric scalar greater than 1. Futility threshold on the \\BF\_{01}\\
  scale.

- p0:

  Numeric scalar in \\(0,1)\\. Null response probability.

- a0, b0:

  Positive numeric scalars. Beta analysis-prior parameters under
  \\H_0\\.

- a1, b1:

  Positive numeric scalars. Beta analysis-prior parameters under
  \\H_1\\.

- dp:

  Optional numeric scalar in \\(0,1)\\. Fixed point alternative used for
  frequentist power calculations under \\H_1\\.

- da0, db0:

  Positive numeric scalars. Beta design-prior parameters under \\H_0\\.

- da1, db1:

  Positive numeric scalars. Beta design-prior parameters under \\H_1\\.

- type:

  Character string specifying the Bayes-factor test. One of `"point"` or
  `"direction"`.

- calibration:

  Character string specifying the calibration mode. One of `"Bayesian"`,
  `"frequentist"`, `"hybrid"`, or `"full"`.

- target_power:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian power.

- target_type1:

  Numeric scalar in \\(0,1)\\. Target corrected Bayesian type-I error.

- target_ce_h0:

  Numeric scalar in \\\[0,1)\\. Optional lower bound on the corrected
  Bayesian probability of compelling evidence in favour of \\H_0\\.

- target_freq_power:

  Numeric scalar in \\(0,1)\\. Target corrected frequentist power at
  `dp`.

- target_freq_type1:

  Numeric scalar in \\(0,1)\\. Target corrected frequentist type-I error
  at \\p = p_0\\.

- algorithm:

  Character string specifying whether the design should be optimized or
  only evaluated.

- interim:

  Optional integer interim sample size used when `algorithm = "manual"`.

- final:

  Optional integer final sample size used when `algorithm = "manual"`.

- power_cushion:

  Optional additive cushion applied only in the fixed-sample anchor
  search of the first optimization step. This can be useful because
  introducing an interim futility analysis typically reduces corrected
  power relative to the fixed-sample anchor.

- ...:

  Reserved for future extensions.

## Value

An object of class `"singlearm_bf_design"`.

## Details

The design uses the Bayes factor \\BF\_{01}\\. Small values of
\\BF\_{01}\\ indicate evidence against \\H_0\\, so final efficacy is
concluded when \\BF\_{01} \\le k\\. Large values indicate evidence in
favour of \\H_0\\, so interim futility is concluded when \\BF\_{01} \\ge
k_f\\.

Analysis priors are specified separately under \\H_0\\ and \\H_1\\ via
`a0, b0, a1, b1`. Design priors are specified separately under \\H_0\\
and \\H_1\\ via `da0, db0, da1, db1`.
