# Check asymptotic feasibility of a ROPE-based design

Computes the asymptotic (n -\> Inf) ceiling on Bayesian predictive power
under H1 and the asymptotic floor on Bayesian predictive type-I error
under H0, implied by the specified beta design priors and decision
boundary. These quantities provide a necessary (but not sufficient)
condition for the feasibility of a calibrated design, and can be used to
diagnose "No feasible design found" results from
[`design_singlearm_onestage_rope`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_onestage_rope.md)
before running the full numerical root-finding search. See Lemma 1
(non-inferiority / superiority) and Lemma 2 (equivalence) in the
accompanying manuscript.

## Usage

``` r
check_feasibility_rope(
  p0,
  delta,
  da0,
  db0,
  da1,
  db1,
  direction = c("equivalence", "noninferiority", "superiority"),
  target_power = NULL,
  target_type1 = NULL,
  verbose = TRUE
)
```

## Arguments

- p0:

  Benchmark response probability.

- delta:

  ROPE half-width (equivalence), NI margin, or superiority margin.

- da0, db0:

  Design prior parameters under H0, Beta(da0, db0).

- da1, db1:

  Design prior parameters under H1, Beta(da1, db1).

- direction:

  Decision type: "equivalence", "noninferiority", or "superiority".

- target_power:

  Target Bayesian predictive power under H1 (optional).

- target_type1:

  Target Bayesian predictive type-I error under H0 (optional).

- verbose:

  Logical; if TRUE (default), prints a human-readable feasibility
  report.

## Value

A list with components `power_ceiling`, `type1_floor`, `power_feasible`
(logical or NA if no target given), and `type1_feasible` (logical or NA
if no target given).

## Examples

``` r
check_feasibility_rope(
  p0 = 0.30, delta = 0.10,
  da0 = 45, db0 = 105, da1 = 44, db1 = 36,
  direction = "superiority",
  target_power = 0.80, target_type1 = 0.10
)
#> ROPE design feasibility check (asymptotic n -> Inf)
#> Direction: superiority 
#> Power ceiling: 0.9966 (target 0.80) -> OK
#> Type-I floor: 0.0052 (target 0.10) -> OK
```
