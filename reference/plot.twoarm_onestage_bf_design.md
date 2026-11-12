# Plot a one-stage two-arm Bayes factor design

Displays operating characteristics, feasibility, or an overview of
design and analysis priors and operating characteristics.

## Usage

``` r
# S3 method for class 'twoarm_onestage_bf_design'
plot(
  x,
  type = c("old", "oc", "feasibility"),
  h0_design = NULL,
  h0_analysis = NULL,
  ...
)
```

## Arguments

- x:

  An object of class `twoarm_onestage_bf_design`.

- type:

  Plot type: `"old"` for the combined layout, `"oc"` for operating
  characteristics, or `"feasibility"` for feasibility.

- h0_design, h0_analysis:

  Optional named numeric vectors `c(a = ..., b = ...)` specifying beta
  priors for the common response probability under the equality null.
  Defaults are obtained from `x$inputs$a_0_d`, `x$inputs$b_0_d`,
  `x$inputs$a_0_a`, and `x$inputs$b_0_a`. Explicit values affect the
  displayed priors only, not the calibrated design. Not used for
  `test = "BF+-"`.

- ...:

  Additional arguments; currently unused.

## Value

A `ggplot` object for `"oc"` or `"feasibility"`, or a `patchwork` object
for `"old"`.

## Details

The combined layout displays a shared mathematical title above both
prior rows. For `"BF+-"`, the title labels the comparison as H+: p2 \>
p1 versus H-: p1 \>= p2. This label does not change the underlying
calculations.

A single prior legend appears below the second prior row.
Alternative-prior panels display beta basis densities rather than
conditional marginal densities induced by directional restrictions on
the joint parameter space.

Operating-characteristic curves are solid, horizontal target lines are
dashed, and a dotted vertical line marks the first jointly feasible
design in the stored search results. Individual target-attainment sample
sizes appear at the right edge and need not equal the jointly calibrated
sample size.
