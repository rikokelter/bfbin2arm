# Calibrate a one-stage two-arm ROPE non-inferiority design

Convenience wrapper for
[`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md)
with `direction = "noninferiority"`.

## Usage

``` r
design_twoarm_onestage_rope_ni(...)
```

## Arguments

- ...:

  Arguments passed to
  [`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md).
  The `direction` argument is set internally to `"noninferiority"` and
  must not be supplied.

## Value

An object of class `"twoarm_rope_directional_onestage_design"`.

## See also

[`design_twoarm_onestage_rope_directional()`](https://rikokelter.github.io/bfbin2arm/reference/design_twoarm_onestage_rope_directional.md)
for all supported arguments and details of the design calculation.
