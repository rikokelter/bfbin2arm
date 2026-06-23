# Plot a one-stage single-arm BF design

Plot a one-stage single-arm BF design

## Usage

``` r
# S3 method for class 'singlearm_onestage_bf_design'
plot(x, what = c("all", "oc"), legend_pos = "right", legend_inset = 0, ...)
```

## Arguments

- x:

  An object of class `"singlearm_onestage_bf_design"`.

- what:

  Character string; currently one of `"all"` or `"oc"`.

- legend_pos:

  Position passed to [`legend`](https://rdrr.io/r/graphics/legend.html).
  Either a keyword such as `"topright"` or a numeric vector `c(x, y)`.

- legend_inset:

  Numeric inset passed to
  [`legend()`](https://rdrr.io/r/graphics/legend.html) when `legend_pos`
  is a keyword.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.
