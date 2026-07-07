# Plot a one-stage single-arm BF design

Produces a 2x2 figure: top-left: operating characteristic curves;
top-right: table-like summary; bottom-left: design priors under H0 and
H1; bottom-right: analysis priors under H0 and H1.

## Usage

``` r
# S3 method for class 'singlearm_onestage_bf_design'
plot(
  x,
  what = c("all", "oc"),
  legend_pos = "right",
  legend_inset = 0,
  col_h0 = "#0072B2",
  col_h1 = "#D55E00",
  prior_lwd = 2,
  ...
)
```

## Arguments

- x:

  An object of class `"singlearm_onestage_bf_design"`.

- what:

  Character string; one of `"all"` or `"oc"`.

- legend_pos:

  Position passed to
  [`legend()`](https://rdrr.io/r/graphics/legend.html).

- legend_inset:

  Numeric inset for
  [`legend()`](https://rdrr.io/r/graphics/legend.html).

- col_h0:

  Colour for H0 priors in bottom panels.

- col_h1:

  Colour for H1 priors in bottom panels.

- prior_lwd:

  Line width for prior density curves.

- ...:

  Currently unused.

## Value

Invisibly returns `x`.
