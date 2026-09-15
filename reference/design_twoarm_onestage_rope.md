# Calibrate a one-stage two-arm ROPE equivalence design for binary endpoints

Calculates equal-allocation sample sizes for a parallel-group one-stage
equivalence study with binary endpoints. The treatment effect is the
risk difference `pT - pC`. Equivalence is declared when the posterior
probability of `abs(pT - pC) < delta` is at least `gammaeq`; meaningful
non-equivalence is declared when the posterior probability outside the
ROPE is at least `gammadiff`; otherwise the decision is inconclusive.

## Usage

``` r
design_twoarm_onestage_rope(
  nmin,
  nmax,
  nstep = 1L,
  delta,
  gammaeq = 0.9,
  gammadiff = gammaeq,
  analysis_prior_C = c(1, 1),
  analysis_prior_T = c(1, 1),
  design_prior_eq_C = NULL,
  design_prior_eq_T = NULL,
  design_prior_eq = NULL,
  design_prior_ne_C,
  design_prior_ne_T,
  targetpower = NULL,
  targettype1 = NULL,
  targetpce_h0 = NULL,
  sustainn = 1L,
  rel.tol = 1e-10,
  integration = c("adaptive", "quantile"),
  quad_nodes = 128L,
  calibration = c("Bayesian", "frequentist", "hybrid", "full"),
  freq_pC = NULL,
  freq_delta = NULL,
  targetfreqpower = NULL,
  targetfreqtype1 = NULL,
  freq_grid_n = 101L,
  freq_pC_range = c(0, 1),
  parallel = FALSE,
  ncores = 1L,
  parallel_backend = c("auto", "psock", "multicore", "none"),
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

- delta:

  ROPE half-width on the risk-difference scale.

- gammaeq:

  Posterior threshold for equivalence.

- gammadiff:

  Posterior threshold for meaningful non-equivalence.

- analysis_prior_C, analysis_prior_T:

  Length-two Beta analysis-prior shapes.

- design_prior_eq_C, design_prior_eq_T:

  Legacy equivalence-prior Beta kernels.

- design_prior_eq:

  Optional joint baseline-risk/risk-difference design prior.

- design_prior_ne_C, design_prior_ne_T:

  Non-equivalence-prior Beta kernels.

- targetpower:

  Target Bayesian predictive probability of declaring practical
  equivalence under the equivalence design distribution `G1`.

- targettype1:

  Maximum Bayesian predictive probability of falsely declaring practical
  equivalence under the non-equivalence design distribution `G0`.

- targetpce_h0:

  Optional minimum Bayesian predictive probability of compelling
  evidence for meaningful non-equivalence under `G0`. More precisely,
  this is the predictive probability of declaring meaningful difference,
  i.e. the posterior probability of `abs(pT - pC) >= delta` is at least
  `gammadiff`. If supplied, this criterion is enforced for `"Bayesian"`,
  `"hybrid"`, and `"full"` calibration.

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

- freq_pC:

  Control response probability for frequentist pointwise power.

- freq_delta:

  Treatment-control difference for frequentist pointwise power.

- targetfreqpower:

  Target pointwise frequentist equivalence power.

- targetfreqtype1:

  Maximum boundary-grid frequentist type-I error.

- freq_grid_n:

  Boundary-grid size.

- freq_pC_range:

  Control-risk interval for boundary-grid evaluation.

- parallel:

  Logical; if `TRUE`, evaluate candidate per-arm sample sizes in
  parallel.

- ncores:

  Number of workers if `parallel = TRUE`.

- parallel_backend:

  Parallel backend. `"auto"` and `"psock"` use a PSOCK cluster, which is
  supported on Windows, macOS, and Linux. `"none"` evaluates candidates
  sequentially. `"multicore"` is retained for compatibility but
  currently resolves to `"psock"` because forked processes can be unsafe
  with threaded numerical libraries.

- returngrid:

  Return the full operating-characteristic grid.

- returnmatrices:

  Retain matrices for the selected design.

- progress:

  Logical; if `TRUE`, show a progress bar during sequential evaluation.
  Progress display is disabled for parallel evaluation.

- progress_every:

  Positive integer specifying the interval, in evaluated candidate
  sample sizes, at which sequential progress is updated.

## Value

An object of class `"bfbin2armrope2armdesign"`.
