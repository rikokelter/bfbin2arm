# Plot a single-arm two-stage ROPE design

Produces a six-panel figure in a 2-row by 3-column layout:

- Top-left:

  Predictive type-I error and power as functions of the stage-1 sample
  size \\n_1\\ for the fixed maximum \\n^\*\\, with the optimal \\n_1\\
  marked.

- Top-centre:

  Predictive PCE under \\H_0\\ as a function of \\n_1\\ for the fixed
  maximum \\n^\*\\, with the optimal \\n_1\\ marked.

- Top-right:

  Textual summary of the optimal design and its operating
  characteristics, enclosed in a full border.

- Bottom-left:

  Null design prior density with the direction-appropriate ROPE region
  shaded.

- Bottom-centre:

  Alternative design prior density with the direction-appropriate ROPE
  region shaded.

- Bottom-right:

  Analysis prior density with the direction-appropriate ROPE region
  shaded.

## Usage

``` r
# S3 method for class 'singlearm_rope_twostage_design'
plot(x, type = c("default", "interim", "final"), ...)
```

## Arguments

- x:

  An object of class `"singlearm_rope_twostage_design"`.

- type:

  Character string selecting the plot type. One of `"default"`
  (default), `"interim"`, or `"final"`.

- ...:

  Further graphical parameters passed to
  [`par()`](https://rdrr.io/r/graphics/par.html) for the `"interim"` and
  `"final"` types.

## Value

Invisibly returns `x`.

## Details

Additional plot types are available via the `type` argument:

- `"default"`:

  The 2x3 summary layout described above (default).

- `"interim"`:

  Interim posterior probability vs. stage-1 responses, showing the
  continuation region \\C_1\\.

- `"final"`:

  Final posterior probability vs. total responses, showing the
  acceptance region.
