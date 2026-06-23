# Plot a single-arm Bayes factor design

Produces a diagnostic plot for a fitted single-arm two-stage Bayes
factor design. Depending on the available information in the object, the
plot shows the interim-search results, selected operating
characteristics, and the design and analysis priors under \\H_0\\ and
\\H_1\\.

## Usage

``` r
# S3 method for class 'singlearm_bf_design'
plot(x, ...)
```

## Arguments

- x:

  An object of class `"singlearm_bf_design"`.

- ...:

  Currently unused.

## Value

The input object `x`, invisibly.

## See also

[`summary.singlearm_bf_design()`](https://rikokelter.github.io/bfbin2arm/reference/summary.singlearm_bf_design.md),
[`design_singlearm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_bf.md),
[`optimal_twostage_singlearm_bf()`](https://rikokelter.github.io/bfbin2arm/reference/optimal_twostage_singlearm_bf.md)
