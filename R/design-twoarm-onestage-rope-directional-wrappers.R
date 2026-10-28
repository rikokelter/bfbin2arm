#' Calibrate a one-stage two-arm ROPE non-inferiority design
#'
#' Convenience wrapper for [design_twoarm_onestage_rope_directional()] with
#' `direction = "noninferiority"`.
#'
#' @param ... Arguments passed to [design_twoarm_onestage_rope_directional()].
#'   The `direction` argument is set internally to `"noninferiority"` and must
#'   not be supplied.
#' @return An object of class `"twoarm_rope_directional_onestage_design"`.
#' @seealso [design_twoarm_onestage_rope_directional()] for all supported
#'   arguments and details of the design calculation.
#' @export
design_twoarm_onestage_rope_ni <- function(...) {
  design_twoarm_onestage_rope_directional(
    ...,
    direction = "noninferiority"
  )
}

#' Calibrate a one-stage two-arm ROPE superiority design
#'
#' Convenience wrapper for [design_twoarm_onestage_rope_directional()] with
#' `direction = "superiority"`.
#'
#' @param ... Arguments passed to [design_twoarm_onestage_rope_directional()].
#'   The `direction` argument is set internally to `"superiority"` and must
#'   not be supplied.
#' @return An object of class `"twoarm_rope_directional_onestage_design"`.
#' @seealso [design_twoarm_onestage_rope_directional()] for all supported
#'   arguments and details of the design calculation.
#' @export
design_twoarm_onestage_rope_sup <- function(...) {
  design_twoarm_onestage_rope_directional(
    ...,
    direction = "superiority"
  )
}

#' Calibrate a one-stage two-arm nested non-inferiority/superiority ROPE design
#'
#' Convenience wrapper for [design_twoarm_onestage_rope_directional()] with
#' `direction = "nested"`. Superiority is only assessed within the
#' non-inferiority acceptance region, so that a superiority conclusion always
#' also satisfies the non-inferiority evidence requirement. This reproduces
#' the hierarchical testing logic used in trials such as SINGLE, where a
#' design may be powered chiefly for non-inferiority while superiority is
#' evaluated as a prespecified secondary claim.
#'
#' @param ... Arguments passed to [design_twoarm_onestage_rope_directional()].
#'   The `direction` argument is set internally to `"nested"` and must not be
#'   supplied.
#' @return An object of class `"twoarm_rope_directional_onestage_design"`.
#' @seealso [design_twoarm_onestage_rope_directional()] for all supported
#'   arguments and details of the design calculation.
#' @export
design_twoarm_onestage_rope_nisup <- function(...) {
  design_twoarm_onestage_rope_directional(
    ...,
    direction = "nested"
  )
}