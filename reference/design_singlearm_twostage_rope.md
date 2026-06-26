# Calibrate an optimal single-arm two-stage ROPE design

Finds a single-arm two-stage Bayesian design based on the region of
practical equivalence (ROPE) for a binary endpoint, with a single
interim analysis allowing early stopping for futility. The design covers
three decision types via the `direction` argument:

- `"equivalence"`:

  Posterior mass inside the two-sided ROPE \\\[p_0 - \delta,\\ p_0 +
  \delta\]\\ must exceed \\\gamma\_{\mathrm{eq}}\\.

- `"noninferiority"`:

  Posterior probability \\\Pr(p \ge p_0 - \delta \mid Y)\\ must exceed
  \\\gamma\_{\mathrm{eq}}\\.

- `"superiority"`:

  Posterior probability \\\Pr(p \> p_0 + \delta \mid Y)\\ must exceed
  \\\gamma\_{\mathrm{eq}}\\.

## Usage

``` r
design_singlearm_twostage_rope(
  p0,
  delta,
  analysis_prior = c(1, 1),
  design_prior_h0,
  design_prior_h1,
  gamma_1 = 0.5,
  gamma_eq = 0.9,
  gamma_diff = gamma_eq,
  alpha = 0.1,
  power = 0.8,
  nmax = 300L,
  direction = c("equivalence", "noninferiority", "superiority"),
  minimax = FALSE,
  progress = TRUE
)
```

## Arguments

- p0:

  Benchmark response probability.

- delta:

  ROPE half-width (`"equivalence"`), non-inferiority margin
  (`"noninferiority"`), or superiority margin (`"superiority"`). Must be
  a single positive number.

- analysis_prior:

  Numeric vector `c(a, b)` for the \\\mathrm{Beta}(a, b)\\ analysis
  prior on \\p\\. Defaults to `c(1, 1)` (uniform).

- design_prior_h0:

  Numeric vector `c(a, b)` for the null design prior.

- design_prior_h1:

  Numeric vector `c(a, b)` for the alternative design prior.

- gamma_1:

  Interim futility threshold in \\(0, 1)\\: continuation requires the
  interim posterior ROPE probability to exceed `gamma_1`.

- gamma_eq:

  Final evidence threshold in \\(0.5, 1)\\: the appropriate posterior
  ROPE probability must exceed `gamma_eq` to declare equivalence,
  non-inferiority, or superiority.

- gamma_diff:

  Threshold for compelling evidence for \\H_0\\: the complementary
  posterior ROPE probability must exceed `gamma_diff`. Defaults to
  `gamma_eq`.

- alpha:

  Target predictive type-I error level (upper bound).

- power:

  Target predictive power (lower bound).

- nmax:

  Upper bound on the fixed-sample size \\n^\*\\ searched in step 1. An
  informative error is raised if no feasible size is found.

- direction:

  Character string specifying the decision type. One of `"equivalence"`
  (default), `"noninferiority"`, or `"superiority"`.

- minimax:

  Logical. If `TRUE`, minimise \\n\\ (minimax criterion); if `FALSE`
  (default), minimise \\\mathrm{EN}\_0\\ (optimal criterion).

- progress:

  Logical. If `TRUE` (default), print progress messages.

## Value

An object of class `"singlearm_rope_twostage_design"` with components:

- `call`:

  The matched call.

- `p0`, `delta`, `direction`:

  Design parameters.

- `analysis_prior`, `design_prior_h0`, `design_prior_h1`:

  Prior specifications.

- `alpha`, `target_power`:

  Constraint levels.

- `gamma_1`, `gamma_eq`, `gamma_diff`:

  Evidence thresholds.

- `optimality`:

  Either `"optimal"` or `"minimax"`.

- `design`:

  A one-row data frame with the optimal design and its operating
  characteristics.

- `continuation_region`:

  Integer vector of stage-1 response counts in \\\mathcal{C}\_1\\.

- `candidates`:

  Data frame of all feasible two-stage designs, sorted by the primary
  optimality criterion.

## Details

The search proceeds in two steps: (1) find the minimum fixed-sample size
\\n^\*\\ at which the one-stage constraints are satisfied; (2) enumerate
all two-stage splits \\n_1 + n_2 = n^\*\\ and retain those satisfying
the two-stage constraints. The optimal design minimises
\\\mathrm{EN}\_0\\ (or \\n^\*\\ under the minimax criterion) among all
feasible splits.

## See also

[`print.singlearm_rope_twostage_design`](https://rikokelter.github.io/bfbin2arm/reference/print.singlearm_rope_twostage_design.md),
[`summary.singlearm_rope_twostage_design`](https://rikokelter.github.io/bfbin2arm/reference/summary.singlearm_rope_twostage_design.md),
[`plot.singlearm_rope_twostage_design`](https://rikokelter.github.io/bfbin2arm/reference/plot.singlearm_rope_twostage_design.md)
