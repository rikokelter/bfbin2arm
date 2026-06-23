#' Predictive density under H0: p1 = p2 = p
#'
#' Beta-binomial predictive density for data (y1,y2) under H0.
#'
#' @inheritParams twoarmbinbf01
#' @param a_0_d,b_0_d Design-prior parameters for common p under H0.
#'
#' @return Numeric scalar, predictive density.
#' @export
predictiveDensityH0 <- function(y1, y2, n1, n2, a_0_d = 1, b_0_d = 1) {
  exp(lchoose(n1, y1) + lchoose(n2, y2) +
        lbeta(a_0_d + y1 + y2,
                     b_0_d + n1 + n2 - y1 - y2) -
        lbeta(a_0_d, b_0_d))
}

#' Predictive density under H1: p1 != p2
#'
#' Product of two independent Beta-binomial predictive densities.
#'
#' @inheritParams twoarmbinbf01
#' @param a_1_d,b_1_d Design-prior parameters for p1.
#' @param a_2_d,b_2_d Design-prior parameters for p2.
#'
#' @return Numeric scalar, predictive density.
#' @export
predictiveDensityH1 <- function(y1, y2, n1, n2,
                                a_1_d = 1, b_1_d = 1,
                                a_2_d = 1, b_2_d = 1) {
  exp(VGAM::dbetabinom.ab(y1, n1, a_1_d, b_1_d, log = TRUE) +
        VGAM::dbetabinom.ab(y2, n2, a_2_d, b_2_d, log = TRUE))
}

# internal normalizing constant for truncated H+ under the design prior
C_trunc_plus <- function(a_1_d, b_1_d, a_2_d, b_2_d) {
  stats::integrate(function(p2)
    stats::dbeta(p2, a_2_d, b_2_d) *
      stats::pbeta(p2, a_1_d, b_1_d),
    0, 1, rel.tol = 1e-4)$value
}

#' Predictive density under H+: p2 > p1 (truncated prior)
#'
#' @inheritParams predictiveDensityH1
#'
#' @return Numeric scalar, predictive density under H+.
#' @export
predictiveDensityHplus_trunc <- function(y1, y2, n1, n2,
                                         a_1_d = 1, b_1_d = 1,
                                         a_2_d = 1, b_2_d = 1) {
  raw_int <- stats::integrate(
    function(p1)
      stats::dbeta(p1, y1 + a_1_d, n1 - y1 + b_1_d) *
      (1 - stats::pbeta(p1, y2 + a_2_d, n2 - y2 + b_2_d)),
    0, 1, rel.tol = 1e-4
  )$value
  
  pred_untr <- exp(
    VGAM::dbetabinom.ab(y1, n1, a_1_d, b_1_d, log = TRUE) +
      VGAM::dbetabinom.ab(y2, n2, a_2_d, b_2_d, log = TRUE)
  ) * raw_int
  
  pred_untr / C_trunc_plus(a_1_d, b_1_d, a_2_d, b_2_d)
}

# internal normalizing constant for truncated H- under the design prior
C_trunc_minus <- function(a_1_d, b_1_d, a_2_d, b_2_d) {
  stats::integrate(function(p1)
    stats::dbeta(p1, a_1_d, b_1_d) *
      stats::pbeta(p1, a_2_d, b_2_d),
    0, 1, rel.tol = 1e-4)$value
}

#' Predictive density under H-: p2 <= p1 (truncated prior)
#'
#' @inheritParams predictiveDensityH1
#'
#' @return Numeric scalar, predictive density under H-.
#' @export
predictiveDensityHminus_trunc <- function(y1, y2, n1, n2,
                                          a_1_d = 1, b_1_d = 1,
                                          a_2_d = 1, b_2_d = 1) {
  raw_int <- stats::integrate(
    function(p1)
      stats::dbeta(p1, y1 + a_1_d, n1 - y1 + b_1_d) *
      stats::pbeta(p1, y2 + a_2_d, n2 - y2 + b_2_d),
    0, 1, rel.tol = 1e-4
  )$value
  
  pred_untr <- exp(
    VGAM::dbetabinom.ab(y1, n1, a_1_d, b_1_d, log = TRUE) +
      VGAM::dbetabinom.ab(y2, n2, a_2_d, b_2_d, log = TRUE)
  ) * raw_int
  
  pred_untr / C_trunc_minus(a_1_d, b_1_d, a_2_d, b_2_d)
}
