# Create a zero-centred normal design-prior object

Convenience wrapper around
[`make_prior()`](https://rikokelter.github.io/bfbin2arm/reference/make_prior.md)
for a normal prior on the standardized effect, with support on the
entire real line.

## Usage

``` r
make_design_prior_normal(sd = 0.5)
```

## Arguments

- sd:

  Positive prior standard deviation.

## Value

An object of class `"bf_ttest_prior"`.

## Examples

``` r
make_design_prior_normal(sd = 0.5)
#> Prior: N(0, 0.5^2) 
#>   family:   normal
#>   location: 0
#>   scale:    0.5
#>   support:  [-Inf, Inf]
```
