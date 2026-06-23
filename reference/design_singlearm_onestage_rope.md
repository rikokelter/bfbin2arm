# Calibrate a one-stage single-arm ROPE design for a binomial endpoint

Calibrate a one-stage single-arm ROPE design for a binomial endpoint

## Usage

``` r
design_singlearm_onestage_rope(
  n_min,
  n_max,
  p0,
  delta,
  gamma_eq,
  gamma_diff = gamma_eq,
  direction = c("equivalence", "noninferiority", "superiority"),
  a = 1,
  b = 1,
  da0,
  db0,
  da1,
  db1,
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  dp = NULL,
  target_power = NULL,
  target_type1 = NULL,
  target_pce_h0 = NULL,
  target_freq_power = NULL,
  target_freq_type1 = NULL,
  sustain_n = 1,
  return_grid = TRUE
)
```

## Arguments

- n_min:

  Minimum sample size.

- n_max:

  Maximum sample size.

- p0:

  Benchmark response probability.

- delta:

  ROPE half-width (equivalence), NI margin, or superiority margin.

- gamma_eq:

  Posterior probability threshold for accepting H1.

- gamma_diff:

  Posterior probability threshold for compelling evidence for H0.
  Defaults to gamma_eq.

- direction:

  Decision type: "equivalence", "noninferiority", or "superiority".

- a, b:

  Analysis prior parameters for Beta(a,b).

- da0, db0:

  Design prior parameters under H0.

- da1, db1:

  Design prior parameters under H1.

- calibration:

  Calibration mode: "Bayesian", "frequentist", "hybrid", or "full".

- dp:

  Point alternative in the favorable H1 region at which frequentist
  power is computed.

- target_power:

  Target Bayesian predictive power under H1.

- target_type1:

  Target Bayesian predictive type-I error under H0.

- target_pce_h0:

  Optional target for predictive compelling evidence for H0 under H0.

- target_freq_power:

  Target frequentist power at dp.

- target_freq_type1:

  Target worst-case frequentist type-I error at the null boundary.

- sustain_n:

  Number of consecutive feasible sample sizes required.

- return_grid:

  Return the full evaluation grid.

## Value

An object of class `bfbin2arm_rope_design`.
