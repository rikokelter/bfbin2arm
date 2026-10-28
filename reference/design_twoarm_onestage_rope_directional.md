# Calibrate a one-stage two-arm directional ROPE design for binary endpoints

Calculates equal-allocation sample sizes for a parallel-group one-stage
study with binary endpoints, where the treatment effect is the risk
difference `pT - pC`. Depending on `direction`, the design targets
non-inferiority, superiority, or a nested procedure in which superiority
is assessed only after non-inferiority has been established. Superiority
is declared when the posterior probability of `pT - pC > delta_sup` is
at least `gamma_sup`; non-inferiority is declared when the posterior
probability of `pT - pC > -delta_ni` is at least `gamma_ni`; harmful
inferiority is declared when the posterior probability of
`pT - pC <= -delta_ni` is at least `gamma_inf`; otherwise the decision
is inconclusive.

## Usage

``` r
design_twoarm_onestage_rope_directional(
  nmin,
  nmax,
  nstep = 1L,
  direction = c("noninferiority", "superiority", "nested"),
  delta_ni,
  delta_sup = 0,
  gamma_ni = 0.975,
  gamma_sup = 0.975,
  gamma_inf = 0.975,
  analysis_prior_C = c(1, 1),
  analysis_prior_T = c(1, 1),
  design_prior_inf = NULL,
  design_prior_ni = NULL,
  design_prior_sup = NULL,
  design_prior_not_sup = NULL,
  target_ni_power = NULL,
  target_sup_power = NULL,
  target_false_ni = NULL,
  target_false_sup = NULL,
  sustainn = 1L,
  rel.tol = 1e-10,
  integration = c("adaptive", "quantile"),
  quad_nodes = 128L,
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  freq_ni_pC = NULL,
  freq_ni_delta = NULL,
  freq_sup_pC = NULL,
  freq_sup_delta = NULL,
  target_freq_ni_power = NULL,
  target_freq_sup_power = NULL,
  target_freq_false_ni = NULL,
  target_freq_false_sup = NULL,
  freq_grid_n = 101L,
  freq_pC_range = c(0, 1),
  parallel = FALSE,
  ncores = 1L,
  parallel_backend = c("auto", "psock", "none"),
  returngrid = TRUE,
  returnmatrices = FALSE,
  progress = interactive(),
  progress_every = 1L
)
```

## Arguments

- nmin, nmax:

  Minimum and maximum per-arm sample size.

- nstep:

  Positive per-arm sample-size increment.

- direction:

  One of `"noninferiority"`, `"superiority"`, or `"nested"`.

- delta_ni:

  Non-inferiority margin on the risk-difference scale.

- delta_sup:

  Superiority margin on the risk-difference scale, `>= 0` by convention,
  with `delta_sup > -delta_ni`.

- gamma_ni:

  Posterior threshold for non-inferiority.

- gamma_sup:

  Posterior threshold for superiority; must be at least `gamma_ni`.

- gamma_inf:

  Posterior threshold for harmful inferiority.

- analysis_prior_C, analysis_prior_T:

  Length-two Beta analysis-prior shapes.

- design_prior_inf, design_prior_ni, design_prior_sup,
  design_prior_not_sup:

  Baseline-risk/risk-difference design priors under the
  harmful-inferiority, non-inferiority, superiority, and not-superiority
  design distributions, respectively. Required depending on `direction`
  and the requested targets.

- target_ni_power:

  Target predictive probability of a non-inferiority or superiority
  conclusion, evaluated under `design_prior_ni`.

- target_sup_power:

  Target predictive probability of a superiority conclusion, evaluated
  under `design_prior_sup`.

- target_false_ni:

  Maximum predictive probability of a non-inferiority or superiority
  conclusion, evaluated under `design_prior_inf`.

- target_false_sup:

  Maximum predictive probability of a superiority conclusion, evaluated
  under `design_prior_not_sup`.

- sustainn:

  Number of consecutive feasible evaluated sample sizes.

- rel.tol:

  Adaptive-integration relative tolerance.

- integration:

  Either `"adaptive"` or `"quantile"`.

- quad_nodes:

  Number of Gauss-Legendre nodes for quantile quadrature.

- calibration:

  One of `"Bayesian"`, `"frequentist"`, `"hybrid"`, or `"full"`.

- freq_ni_pC, freq_ni_delta:

  Control response and true difference for pointwise frequentist
  non-inferiority power.

- freq_sup_pC, freq_sup_delta:

  Control response and true difference for pointwise frequentist
  superiority power.

- target_freq_ni_power, target_freq_sup_power:

  Target pointwise frequentist powers.

- target_freq_false_ni, target_freq_false_sup:

  Maximum boundary-grid frequentist error rates.

- freq_grid_n:

  Boundary-grid size.

- freq_pC_range:

  Control-risk interval for boundary-grid evaluation.

- parallel:

  Logical; if `TRUE`, evaluate candidate per-arm sample sizes in
  parallel.

- ncores:

  Positive integer giving the maximum number of worker processes used
  when `parallel = TRUE`.

- parallel_backend:

  Parallel backend. `"auto"` and `"psock"` use a PSOCK cluster. `"none"`
  evaluates sequentially.

- returngrid:

  Return the full operating-characteristic grid.

- returnmatrices:

  Retain matrices for the selected design.

- progress:

  Logical; if `TRUE`, show a progress bar during sequential evaluation.

- progress_every:

  Positive integer specifying the interval, in evaluated candidate
  sample sizes, at which sequential progress is updated.

## Value

An object of class `"twoarm_rope_directional_onestage_design"`.
