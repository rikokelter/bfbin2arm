# Design or evaluate a one-stage single-arm Bayes factor trial

Calibrates or evaluates a one-stage single-arm Bayes factor design for a
binary endpoint.

## Usage

``` r
design_singlearm_onestage_bf(
  n_min,
  n_max,
  k,
  k_ce = NULL,
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
  n = NULL,
  power_cushion = 0,
  sustain_n = 10L,
  ...
)
```

## Arguments

- n_min:

  Integer. Minimum admissible sample size.

- n_max:

  Integer. Maximum admissible sample size.

- k:

  Numeric scalar greater than 0. Evidence threshold on the \\BF\_{01}\\
  scale for efficacy, used for power and type-I error.

- k_ce:

  Optional numeric scalar greater than 1. Threshold on the \\BF\_{01}\\
  scale used for CE(H0) / PCE(H0). Must be supplied when
  `target_ce_h0 > 0`.

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

- n:

  Optional integer sample size used when `algorithm = "manual"`.

- power_cushion:

  Optional additive cushion applied to the power targets in the
  optimizer.

- sustain_n:

  Non-negative integer. A candidate design is considered feasible only
  if the relevant operating characteristics satisfy their target
  constraints at the candidate sample size and for the next `sustain_n`
  larger sample sizes, subject to the search range. This also applies to
  the CE(H0) constraint when `target_ce_h0 > 0`.

- ...:

  Reserved for future extensions.

## Value

An object of class `"singlearm_onestage_bf_design"`.

## Details

The design uses the Bayes factor \\BF\_{01}\\. Small values of
\\BF\_{01}\\ indicate evidence against \\H_0\\, so efficacy is concluded
when \\BF\_{01} \le k\\. Large values indicate evidence in favour of
\\H_0\\, and the optional CE(H0) / PCE(H0) constraint is evaluated using
the separate threshold `k_ce`.

Analysis priors are specified separately under \\H_0\\ and \\H_1\\ via
`a0, b0, a1, b1`. Design priors are specified separately under \\H_0\\
and \\H_1\\ via `da0, db0, da1, db1`.
