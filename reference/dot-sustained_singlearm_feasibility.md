# Check sustained feasibility over future n

Given vectors of operating characteristics over n, check whether power,
type-I-error, and CE(H0) satisfy their thresholds at n and for at least
sustain_n subsequent sample sizes.

## Usage

``` r
.sustained_singlearm_feasibility(
  n_vec,
  power_vec,
  type1_vec,
  ce_vec,
  target_power,
  target_type1,
  target_ce,
  sustain_n
)
```

## Arguments

- n_vec:

  Integer vector of sample sizes.

- power_vec:

  Numeric vector of power values (same length as n_vec).

- type1_vec:

  Numeric vector of type-I-error values.

- ce_vec:

  Numeric vector of CE(H0) values (may contain NA).

- target_power:

  Numeric target power.

- target_type1:

  Numeric target type-I-error.

- target_ce:

  Numeric target CE(H0); if \<= 0, CE constraint is ignored.

- sustain_n:

  Integer, number of subsequent sample sizes that must also satisfy the
  constraints.

## Value

Logical vector of length length(n_vec): TRUE if n_i and the next
sustain_n sample sizes satisfy all active constraints.
