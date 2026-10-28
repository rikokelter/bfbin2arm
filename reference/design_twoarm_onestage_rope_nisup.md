# Calibrate a one-stage two-arm nested non-inferiority/superiority ROPE design

Convenience wrapper for
[`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md)
with `direction = "nested"`. Superiority is only assessed within the
non-inferiority acceptance region, so that a superiority conclusion
always also satisfies the non-inferiority evidence requirement. This
reproduces the hierarchical testing logic used in trials such as SINGLE,
where a design may be powered chiefly for non-inferiority while
superiority is evaluated as a prespecified secondary claim.

## Usage

``` r
design_twoarm_onestage_rope_nisup(...)
```

## Arguments

- ...:

  Arguments passed to
  [`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md).
  The `direction` argument is set internally to `"nested"` and must not
  be supplied.

## Value

An object of class `"twoarm_rope_directional_onestage_design"`.

## See also

[`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md)
for all supported arguments and details of the design calculation.
