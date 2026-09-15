#' Calibrate a two-stage ROPE equivalence design
#'
#' Convenience wrapper for [design_singlearm_twostage_rope()] with
#' `direction = "equivalence"`.
#'
#' @param ... Arguments passed to [design_singlearm_twostage_rope()]. The
#'   `direction` argument is set internally to `"equivalence"` and must not be
#'   supplied.
#' @return An object of class `"singlearm_rope_twostage_design"`.
#' @seealso [design_singlearm_twostage_rope()] for all supported arguments and
#'   details of the design calculation.
#' @export
design_singlearm_twostage_rope_equiv <- function(...) {
  design_singlearm_twostage_rope(
    ...,
    direction = "equivalence"
  )
}

#' Calibrate a two-stage ROPE non-inferiority design
#'
#' Convenience wrapper for [design_singlearm_twostage_rope()] with
#' `direction = "noninferiority"`.
#'
#' @param ... Arguments passed to [design_singlearm_twostage_rope()]. The
#'   `direction` argument is set internally to `"noninferiority"` and must not
#'   be supplied.
#' @return An object of class `"singlearm_rope_twostage_design"`.
#' @seealso [design_singlearm_twostage_rope()] for all supported arguments and
#'   details of the design calculation.
#' @export
design_singlearm_twostage_rope_ni <- function(...) {
  design_singlearm_twostage_rope(
    ...,
    direction = "noninferiority"
  )
}

#' Calibrate a two-stage ROPE superiority design
#'
#' Convenience wrapper for [design_singlearm_twostage_rope()] with
#' `direction = "superiority"`.
#'
#' @param ... Arguments passed to [design_singlearm_twostage_rope()]. The
#'   `direction` argument is set internally to `"superiority"` and must not be
#'   supplied.
#' @return An object of class `"singlearm_rope_twostage_design"`.
#' @seealso [design_singlearm_twostage_rope()] for all supported arguments and
#'   details of the design calculation.
#' @export
design_singlearm_twostage_rope_sup <- function(...) {
  design_singlearm_twostage_rope(
    ...,
    direction = "superiority"
  )
}