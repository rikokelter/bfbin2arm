# Plot method for calibrated Bayesian t-test designs

Plot method for calibrated Bayesian t-test designs

## Usage

``` r
# S3 method for class 'bf_ttest_design'
plot(
  x,
  which = c("decision", "analysis_prior", "design_prior", "n_curve"),
  annotate = TRUE,
  prior_points = 801L,
  ...
)
```

## Arguments

- x:

  Object of class `bfttestdesign`.

- which:

  Character vector choosing panels among `"decision"`,
  `"analysisprior"`, `"designprior"`, `"effect"`, and `"ncurve"`.

- annotate:

  Logical; add global title/subtitle.

- prior_points:

  Number of grid points for prior density curves.

- ...:

  Ignored.

## Value

A `ggplot` object or a patchwork object.
