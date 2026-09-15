# Calibrate a two-stage ROPE superiority design

Convenience wrapper for
[`design_singlearm_twostage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_twostage_rope.md)
with `direction = "superiority"`.

## Usage

``` r
design_singlearm_twostage_rope_sup(...)
```

## Arguments

- ...:

  Arguments passed to
  [`design_singlearm_twostage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_twostage_rope.md).
  The `direction` argument is set internally to `"superiority"` and must
  not be supplied.

## Value

An object of class `"singlearm_rope_twostage_design"`.

## See also

[`design_singlearm_twostage_rope()`](https://rikokelter.github.io/bfbin2arm/reference/design_singlearm_twostage_rope.md)
for all supported arguments and details of the design calculation.
