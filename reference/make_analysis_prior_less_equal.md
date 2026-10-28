# Create the default non-positive analysis prior

Constructs a zero-centred normal prior for the standardized effect,
truncated to non-positive effects. This is the default H0 prior for
`hypothesis = "greater.composite"`.

## Usage

``` r
make_analysis_prior_less_equal(sd = 0.5)
```

## Arguments

- sd:

  Positive prior standard deviation.

## Value

An object of class `"bf_ttest_prior"` with support `c(-Inf, 0)`.

## Examples

``` r
make_analysis_prior_less_equal(sd = 0.5)
#> Prior: TN(0, 0.5^2; (-Inf, 0]) 
#>   family:   normal
#>   location: 0
#>   scale:    0.5
#>   support:  [-Inf, 0]
```
