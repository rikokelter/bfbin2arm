# Create a prior for the standardized effect

Create a prior for the standardized effect

## Usage

``` r
make_prior(
  family = c("normal", "student", "cauchy", "jzs", "custom"),
  location = 0,
  scale = 1,
  df = 1,
  support = c(-Inf, Inf),
  log_density = NULL,
  density = NULL,
  name = NULL
)
```

## Arguments

- family:

  Prior family.

- location:

  Location parameter.

- scale:

  Scale parameter (\> 0).

- df:

  Degrees of freedom for Student prior.

- support:

  Numeric vector `c(lower, upper)` defining support.

- log_density:

  Optional custom log-density.

- density:

  Optional custom density.

- name:

  Optional display name.

## Value

Object of class `bf_ttest_prior`.
