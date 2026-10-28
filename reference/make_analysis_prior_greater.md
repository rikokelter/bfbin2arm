# Create the default positive-direction analysis prior

Constructs a zero-centred normal prior for the standardized effect,
truncated to positive effects. This is the default alternative prior for
`hypothesis = "greater.pointnull"` and the H1 prior for
`hypothesis = "greater.composite"`.

## Usage

``` r
make_analysis_prior_greater(sd = 0.5)
```

## Arguments

- sd:

  Positive prior standard deviation.

## Value

An object of class `"bf_ttest_prior"` with support `c(0, Inf)`.

## Examples

``` r
make_analysis_prior_greater(sd = 0.5)
#> Prior: TN(0, 0.5^2; [0, Inf)) 
#>   family:   normal
#>   location: 0
#>   scale:    0.5
#>   support:  [0, Inf]
```
