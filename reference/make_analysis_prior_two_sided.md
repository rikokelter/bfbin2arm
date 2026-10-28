# Create the default two-sided analysis prior

Constructs a zero-centred normal prior for the standardized effect with
support on the entire real line.

## Usage

``` r
make_analysis_prior_two_sided(sd = 0.5)
```

## Arguments

- sd:

  Positive prior standard deviation.

## Value

An object of class `"bf_ttest_prior"`.

## Examples

``` r
make_analysis_prior_two_sided()
#> Prior: N(0, 0.5^2) 
#>   family:   normal
#>   location: 0
#>   scale:    0.5
#>   support:  [-Inf, Inf]
```
