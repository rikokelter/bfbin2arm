#' Bayesian power, type-I error, and PCE(H0) for two-arm binomial Bayes factors
#'
#' Computes Bayesian power, Bayesian type-I error, and the probability of
#' compelling evidence under H_0 (or H_- for BF+-), for a given
#' sample size and Bayes factor test. Optionally, frequentist type-I error
#' and frequentist power are computed by summing over the rejection region.
#'
#' @param n1,n2 Sample sizes in arms 1 and 2.
#' @param k Evidence threshold for rejecting the null (inverted BF).
#' @param k_f Evidence threshold for "compelling evidence" in favour of the null.
#' @param test Character string, one of \code{"BF01"}, \code{"BF+0"},
#'   \code{"BF-0"}, \code{"BF+-"}.
#' @param a_0_d,b_0_d,a_0_a,b_0_a Shape parameters for design and analysis priors
#'   under \eqn{H_0}.
#' @param a_1_d,b_1_d,a_2_d,b_2_d Shape parameters for design priors under
#'   \eqn{H_1} or \eqn{H_+}.
#' @param a_1_a,b_1_a,a_2_a,b_2_a Shape parameters for analysis priors under
#'   \eqn{H_1} or \eqn{H_+}.
#' @param a_1_d_Hminus,b_1_d_Hminus,a_2_d_Hminus,b_2_d_Hminus Optional design
#'   priors under \eqn{H_-} for directional tests.
#' @param output One of \code{"numeric"}, \code{"predDensmatrix"},
#'   \code{"t1ematrix"}, \code{"ceH0matrix"}, \code{"frequentist_t1e"}.
#' @param compute_freq_t1e Logical; if \code{TRUE}, compute frequentist
#'   type-I error over a grid.
#' @param p1_grid,p2_grid Grids of true proportions for frequentist T1E.
#' @param p1_power,p2_power Optional true proportions for frequentist power.
#' @param a_1_a_Hminus,b_1_a_Hminus,a_2_a_Hminus,b_2_a_Hminus Shape parameters
#'   for analysis priors under \eqn{H_-} (directional tests).
#'
#' @return Depending on \code{output}, either a named numeric vector with
#'   components \code{Power}, \code{Type1_Error}, \code{CE_H0} (and optionally
#'   frequentist metrics) or matrices of predictive densities.
#' @examples
#' # Basic Bayesian power for BF01 test
#' powertwoarmbinbf01(n1 = 30, n2 = 30, k = 1/3, test = "BF01")
#'
#' # Directional test BF+0 with frequentist type-I error
#' powertwoarmbinbf01(n1 = 40, n2 = 40, k = 1/3, k_f = 3,
#'                    test = "BF+0", compute_freq_t1e = TRUE)
#'
#' # Predictive density matrices (advanced)
#' powertwoarmbinbf01(n1 = 25, n2 = 25, output = "predDensmatrix")
#' @export
powertwoarmbinbf01 <- function(
    n1, n2,
    k = 1/3, k_f = 1/3,
    test = c("BF01", "BF+0", "BF-0", "BF+-"),
    a_0_d = 1, b_0_d = 1,
    a_0_a = 1, b_0_a = 1,
    a_1_d = 1, b_1_d = 1,
    a_2_d = 1, b_2_d = 1,
    a_1_a = 1, b_1_a = 1,
    a_2_a = 1, b_2_a = 1,
    output = c("numeric", "predDensmatrix", "t1ematrix", "ceH0matrix", "frequentist_t1e"),
    a_1_d_Hminus = 1, b_1_d_Hminus = 1,
    a_2_d_Hminus = 1, b_2_d_Hminus = 1,
    compute_freq_t1e = FALSE,
    p1_grid = seq(0.01, 0.99, 0.02),
    p2_grid = seq(0.01, 0.99, 0.02),
    p1_power = NULL, p2_power = NULL,
    a_1_a_Hminus = 1, b_1_a_Hminus = 1,
    a_2_a_Hminus = 1, b_2_a_Hminus = 1
) {
  test   <- match.arg(test)
  output <- match.arg(output)
  
  # ============================================================
  #  Hypothesis description text (full mathematical form)
  # ============================================================
  test_description <- switch(
    test,
    "BF01" = "H[1]:~p[1] != p[2] ~~ vs ~~ H[0]:~p[1] == p[2]",
    "BF+0" = "H[+]:~p[2] > p[1] ~~ vs ~~ H[0]:~p[1] == p[2]",
    "BF-0" = "H[-]:~p[2] <= p[1] ~~ vs ~~ H[0]:~p[1] == p[2]",
    "BF+-" = "H[+]:~p[2] > p[1] ~~ vs ~~ H[-]:~p[2] <= p[1]"
  )
  
  # -------------------------------------------------------------------------
  # internal Bayes factor helper functions (analysis priors only)
  # -------------------------------------------------------------------------
  twoarmbinbf01_internal <- function(y1, y2, n1, n2,
                                     a_0_a, b_0_a,
                                     a_1_a, b_1_a,
                                     a_2_a, b_2_a) {
    numerator <- beta(a_0_a + y1 + y2,
                      b_0_a + n1 + n2 - y1 - y2) / beta(a_0_a, b_0_a)
    denominator <- (beta(a_1_a + y1, b_1_a + n1 - y1) *
                      beta(a_2_a + y2, b_2_a + n2 - y2)) /
      (beta(a_1_a, b_1_a) * beta(a_2_a, b_2_a))
    numerator / denominator
  }
  
  priorProbHplus_internal <- function(a_1_a, b_1_a, a_2_a, b_2_a)
    stats::integrate(
      function(p2) stats::dbeta(p2, a_2_a, b_2_a) *
        stats::pbeta(p2, a_1_a, b_1_a),
      lower = 0, upper = 1, rel.tol = 1e-4
    )$value
  
  priorProbHminus_internal <- function(a_1_a, b_1_a, a_2_a, b_2_a)
    1 - priorProbHplus_internal(a_1_a, b_1_a, a_2_a, b_2_a)
  
  postProbHplus_internal <- function(y1, y2, n1, n2,
                                     a_1_a, b_1_a, a_2_a, b_2_a)
    stats::integrate(
      function(p2)
        stats::dbeta(p2, a_2_a + y2, b_2_a + n2 - y2) *
        stats::pbeta(p2, a_1_a + y1, b_1_a + n1 - y1),
      lower = 0, upper = 1, rel.tol = 1e-4
    )$value
  
  postProbHminus_internal <- function(y1, y2, n1, n2,
                                      a_1_a, b_1_a, a_2_a, b_2_a)
    1 - postProbHplus_internal(y1, y2, n1, n2, a_1_a, b_1_a, a_2_a, b_2_a)
  
  BFplus1_internal <- function(y1, y2, n1, n2,
                               a_1_a, b_1_a, a_2_a, b_2_a)
    postProbHplus_internal(y1, y2, n1, n2, a_1_a, b_1_a, a_2_a, b_2_a) /
    priorProbHplus_internal(a_1_a, b_1_a, a_2_a, b_2_a)
  
  BFminus1_internal <- function(y1, y2, n1, n2,
                                a_1_a, b_1_a, a_2_a, b_2_a)
    postProbHminus_internal(y1, y2, n1, n2, a_1_a, b_1_a, a_2_a, b_2_a) /
    priorProbHminus_internal(a_1_a, b_1_a, a_2_a, b_2_a)
  
  # -------------------------------------------------------------------------
  # predictive densities (design priors)
  # -------------------------------------------------------------------------
  predictiveDensityH0 <- function(y1, y2, n1, n2, a_0_d, b_0_d) {
    exp(lchoose(n1, y1) + lchoose(n2, y2) +
          lbeta(a_0_d + y1 + y2, b_0_d + n1 + n2 - y1 - y2) -
          lbeta(a_0_d, b_0_d))
  }
  
  predictiveDensityH1 <- function(y1, y2, n1, n2, a_1_d, b_1_d, a_2_d, b_2_d) {
    exp(VGAM::dbetabinom.ab(y1, n1, a_1_d, b_1_d, log = TRUE) +
          VGAM::dbetabinom.ab(y2, n2, a_2_d, b_2_d, log = TRUE))
  }
  
  C_trunc <- function(a_1_d, b_1_d, a_2_d, b_2_d)
    stats::integrate(
      function(p2) stats::dbeta(p2, a_2_d, b_2_d) *
        stats::pbeta(p2, a_1_d, b_1_d),
      0, 1, rel.tol = 1e-4
    )$value
  
  predictiveDensityHplus_trunc <- function(y1, y2, n1, n2,
                                           a_1_d, b_1_d, a_2_d, b_2_d) {
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
    pred_untr / C_trunc(a_1_d, b_1_d, a_2_d, b_2_d)
  }
  
  C_trunc_Hminus <- function(a_1_d, b_1_d, a_2_d, b_2_d)
    stats::integrate(
      function(p1) stats::dbeta(p1, a_1_d, b_1_d) *
        stats::pbeta(p1, a_2_d, b_2_d),
      0, 1, rel.tol = 1e-4
    )$value
  
  predictiveDensityHminus_trunc <- function(y1, y2, n1, n2,
                                            a_1_d, b_1_d, a_2_d, b_2_d) {
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
    pred_untr / C_trunc_Hminus(a_1_d, b_1_d, a_2_d, b_2_d)
  }
  
  # -------------------------------------------------------------------------
  # Frequentist Type-I error and power
  # -------------------------------------------------------------------------
  sup_freq_t1e <- NA_real_
  freq_power   <- NA_real_
  
  reject_region <- matrix(FALSE, n1 + 1, n2 + 1)
  
  # First pass: identify rejection region
  for (i in 1:(n1 + 1)) {
    for (j in 1:(n2 + 1)) {
      y1 <- i - 1; y2 <- j - 1
      
      if (test == "BF01") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        reject_region[i, j] <- (BF01 < k)
        
      } else if (test == "BF+0") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFp1 <- BFplus1_internal(
          y1, y2, n1, n2,
          a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFp0 <- BFp1 / BF01
        BF0p <- 1 / BFp0
        reject_region[i, j] <- (BF0p < k)
        
      } else if (test == "BF-0") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFm1 <- BFminus1_internal(
          y1, y2, n1, n2,
          a_1_a_Hminus, b_1_a_Hminus,
          a_2_a_Hminus, b_2_a_Hminus
        )
        BFm0 <- BFm1 / BF01
        BF0m <- 1 / BFm0
        reject_region[i, j] <- (BF0m < k)
        
      } else if (test == "BF+-") {
        BFp1 <- BFplus1_internal(
          y1, y2, n1, n2,
          a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFm1 <- BFminus1_internal(
          y1, y2, n1, n2,
          a_1_a_Hminus, b_1_a_Hminus,
          a_2_a_Hminus, b_2_a_Hminus
        )
        BFpm <- BFp1 / BFm1
        BFmp <- 1 / BFpm
        reject_region[i, j] <- (BFmp < k)  # reject H- in favour of H+
      }
    }
  }
  
  # Second pass: supremum Type-I error over null region
  if (compute_freq_t1e) {
    freq_t1e_grid <- matrix(0, length(p1_grid), length(p2_grid))
    for (pi in seq_along(p1_grid)) {
      for (pj in seq_along(p2_grid)) {
        p1_true <- p1_grid[pi]
        p2_true <- p2_grid[pj]
        
        if (test %in% c("BF01", "BF+0", "BF-0")) {
          if (abs(p1_true - p2_true) < 1e-6) {
            freq_t1e_grid[pi, pj] <- sum(
              reject_region *
                stats::dbinom(0:n1, n1, p1_true) %o%
                stats::dbinom(0:n2, n2, p2_true)
            )
          }
        } else if (test == "BF+-") {
          if (p2_true <= p1_true) {
            freq_t1e_grid[pi, pj] <- sum(
              reject_region *
                stats::dbinom(0:n1, n1, p1_true) %o%
                stats::dbinom(0:n2, n2, p2_true)
            )
          }
        }
      }
    }
    sup_freq_t1e <- max(freq_t1e_grid, na.rm = TRUE)
  }
  
  # Frequentist power for specified p1_power, p2_power
  if (!is.null(p1_power) && !is.null(p2_power)) {
    freq_power <- sum(
      reject_region *
        stats::dbinom(0:n1, n1, p1_power) %o%
        stats::dbinom(0:n2, n2, p2_power)
    )
  }
  
  # -------------------------------------------------------------------------
  # matrices for Bayesian calculations
  # -------------------------------------------------------------------------
  BFmat       <- matrix(0, n1 + 1, n2 + 1)
  BFmat_t1e   <- matrix(0, n1 + 1, n2 + 1)
  BFmat_ceH0  <- matrix(0, n1 + 1, n2 + 1)
  
  # main grid loop (Bayesian power)
  for (i in 1:(n1 + 1)) {
    for (j in 1:(n2 + 1)) {
      y1 <- i - 1; y2 <- j - 1
      
      if (test == "BF01") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        if (BF01 < k) {
          BFmat[i, j]      <- predictiveDensityH1(y1, y2, n1, n2, a_1_d, b_1_d, a_2_d, b_2_d)
          BFmat_t1e[i, j]  <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
        if (BF01 > k_f) {
          BFmat_ceH0[i, j] <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
      }
      
      if (test == "BF+0") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFp1 <- BFplus1_internal(
          y1, y2, n1, n2,
          a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFp0 <- BFp1 / BF01
        BF0p <- 1 / BFp0
        if (BF0p < k) {
          BFmat[i, j]     <- predictiveDensityHplus_trunc(y1, y2, n1, n2, a_1_d, b_1_d, a_2_d, b_2_d)
          BFmat_t1e[i, j] <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
        if (BF0p > k_f) {
          BFmat_ceH0[i, j] <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
      }
      
      if (test == "BF-0") {
        BF01 <- twoarmbinbf01_internal(
          y1, y2, n1, n2,
          a_0_a, b_0_a, a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFm1 <- BFminus1_internal(
          y1, y2, n1, n2,
          a_1_a_Hminus, b_1_a_Hminus,
          a_2_a_Hminus, b_2_a_Hminus
        )
        BFm0 <- BFm1 / BF01
        BF0m <- 1 / BFm0
        if (BF0m < k) {
          BFmat[i, j]     <- predictiveDensityHminus_trunc(
            y1, y2, n1, n2,
            a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus
          )
          BFmat_t1e[i, j] <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
        if (BF0m > k_f) {
          BFmat_ceH0[i, j] <- predictiveDensityH0(y1, y2, n1, n2, a_0_d, b_0_d)
        }
      }
      
      if (test == "BF+-") {
        BFp1 <- BFplus1_internal(
          y1, y2, n1, n2,
          a_1_a, b_1_a, a_2_a, b_2_a
        )
        BFm1 <- BFminus1_internal(
          y1, y2, n1, n2,
          a_1_a_Hminus, b_1_a_Hminus,
          a_2_a_Hminus, b_2_a_Hminus
        )
        BFpm <- BFp1 / BFm1
        BFmp <- 1 / BFpm
        if (BFmp < k) {
          BFmat[i, j]     <- predictiveDensityHplus_trunc(
            y1, y2, n1, n2, a_1_d, b_1_d, a_2_d, b_2_d
          )
          BFmat_t1e[i, j] <- predictiveDensityHminus_trunc(
            y1, y2, n1, n2,
            a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus
          )
        }
        if (BFmp > k_f) {
          BFmat_ceH0[i, j] <- predictiveDensityHminus_trunc(
            y1, y2, n1, n2,
            a_1_d_Hminus, b_1_d_Hminus, a_2_d_Hminus, b_2_d_Hminus
          )
        }
      }
    }
  }
  
  vals <- c(sum(BFmat), sum(BFmat_t1e), sum(BFmat_ceH0))
  names(vals) <- c("Power", "Type1_Error", "CE_H0")
  
  if (!is.na(sup_freq_t1e)) {
    vals <- c(vals, Frequentist_Type1_Error = sup_freq_t1e)
  }
  if (!is.null(p1_power) && !is.null(p2_power) && !is.na(freq_power)) {
    vals <- c(vals, Frequentist_Power = freq_power)
  }
  
  if (output == "numeric") {
    attr(vals, "hypothesis")        <- test_description
    attr(vals, "compute_freq_t1e")  <- compute_freq_t1e
    attr(vals, "p1_power")          <- p1_power
    attr(vals, "p2_power")          <- p2_power
    return(vals)
  }
  
  if (output == "predDensmatrix") return(round(BFmat,      4))
  if (output == "t1ematrix")      return(round(BFmat_t1e,  4))
  if (output == "ceH0matrix")     return(round(BFmat_ceH0, 4))
  if (output == "frequentist_t1e" && compute_freq_t1e) {
    return(list(
      sup_freq_t1e    = sup_freq_t1e,
      p1_grid         = p1_grid,
      p2_grid         = p2_grid,
      test_description = test_description
    ))
  }
}


#' Sample size calibration for two-arm binomial Bayes factor designs
#'
#' Searches over a grid of total sample sizes n to find the smallest n such that
#' Bayesian power, Bayesian type-I error, and probability of compelling evidence
#' under H0 meet specified design criteria. Optionally, frequentist type-I error
#' and power constraints are also evaluated. Unequal fixed randomisation between
#' the two arms is allowed via alloc1 and alloc2.
#'
#' @inheritParams powertwoarmbinbf01
#' @param power Target Bayesian power.
#' @param alpha Target Bayesian type-I error.
#' @param pce_H0 Target probability of compelling evidence under \eqn{H_0}.
#' @param nrange Integer vector of length 2 giving the search range for total n.
#' @param n_step Step size for n.
#' @param progress Logical; if \code{TRUE}, print progress to the console.
#' @param output \code{"plot"} or \code{"numeric"}.
#' @param a_1_a_Hminus,b_1_a_Hminus,a_2_a_Hminus,b_2_a_Hminus Shape parameters for analysis priors under H-.
#' @param alloc1,alloc2 Fixed randomisation probabilities for arm 1 and arm 2;
#'   must be positive and sum to 1. Defaults are 0.5 and 0.5.
#'
#' @return If \code{output = "plot"}, returns invisibly a list with recommended
#'   sample sizes and a ggplot object printed to the device. If
#'   \code{output = "numeric"}, returns a list with recommended n and summary.
#' @examples
#' # Standard calibration with equal allocation: power 80%, type-I 5%, CE(H0) 80%
#' \donttest{
#' ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.8, output = "numeric")
#' }
#'
#' # 1:2 allocation (control:treatment) via alloc1 = 1/3, alloc2 = 2/3
#' \donttest{
#' ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.8,
#'                alloc1 = 1/3, alloc2 = 2/3, output = "numeric")
#' }
#'
#' # BF+0 directional test with plot
#' \donttest{
#' ntwoarmbinbf01(power = 0.8, alpha = 0.05, pce_H0 = 0.9,
#'                test = "BF+0", output = "plot")
#' }
#' @export
ntwoarmbinbf01 <- function(
    k = 1/3, k_f = 1/3,
    power = 0.8, alpha = 0.05, pce_H0 = 0.9,
    test = c("BF01", "BF+0", "BF-0", "BF+-"),
    nrange = c(10, 150), n_step = 1,
    progress = TRUE, compute_freq_t1e = FALSE,
    p1_grid = seq(0.01, 0.99, 0.02),
    p2_grid = seq(0.01, 0.99, 0.02),
    p1_power = NULL, p2_power = NULL,
    a_0_d = 1, b_0_d = 1, a_0_a = 1, b_0_a = 1,
    a_1_d = 1, b_1_d = 1, a_2_d = 1, b_2_d = 1,
    a_1_a = 1, b_1_a = 1, a_2_a = 1, b_2_a = 1,
    output = c("plot", "numeric"),
    a_1_d_Hminus = 1, b_1_d_Hminus = 1,
    a_2_d_Hminus = 1, b_2_d_Hminus = 1,
    a_1_a_Hminus = 1, b_1_a_Hminus = 1,
    a_2_a_Hminus = 1, b_2_a_Hminus = 1,
    alloc1 = 0.5, alloc2 = 0.5
) {
  test   <- match.arg(test)
  output <- match.arg(output)
  
  ## --- allocation checks ---
  if (!is.numeric(alloc1) || !is.numeric(alloc2) ||
      length(alloc1) != 1L || length(alloc2) != 1L ||
      alloc1 <= 0 || alloc2 <= 0) {
    stop("alloc1 and alloc2 must be positive scalars.")
  }
  alloc_sum <- alloc1 + alloc2
  if (!isTRUE(all.equal(alloc_sum, 1))) {
    warning(sprintf(
      "alloc1 + alloc2 != 1 (%.4f); normalizing to sum to 1.", alloc_sum
    ))
    alloc1 <- alloc1 / alloc_sum
    alloc2 <- alloc2 / alloc_sum
  }
  
  if (!is.null(p1_power) && !is.null(p2_power)) {
    message(sprintf("Frequentist power computation: p1=%.2f, p2=%.2f\n",
                    p1_power, p2_power))
  }
  
  hyp_desc <- switch(
    test,
    "BF01" = "BF01 test: H0: p1 = p2 vs H1: p1 != p2",
    "BF+0" = "BF+0 test: H0: p1 = p2 vs H+: p2 > p1",
    "BF-0" = "BF-0 test: H0: p1 = p2 vs H-: p1 > p2",
    "BF+-" = "BF+- test: H+: p2 > p1 vs H-: p1 > p2"
  )
  
  n_seq  <- seq(nrange[1], nrange[2], by = n_step)
  n_vals <- length(n_seq)
  message(sprintf(
    "Computing for total n = %d to %d (step = %d, %d values)\n",
    nrange[1], nrange[2], n_step, n_vals
  ))
  message(sprintf("Allocation: alloc1 = %.3f, alloc2 = %.3f\n\n", alloc1, alloc2))
  
  if (compute_freq_t1e) message("Frequentist Type-I error computation: ENABLED\n")
  if (!is.null(p1_power) && !is.null(p2_power))
    message("Frequentist power computation: ENABLED\n")
  
  results <- NULL
  freq_t1e_results <- NULL
  freq_power_results <- NULL
  
  # Simulation loop over total n
  for (n_tot in n_seq) {
    n1 <- round(n_tot * alloc1)
    n2 <- n_tot - n1  # ensure n1 + n2 = n_tot
    
    res <- powertwoarmbinbf01(
      n1 = n1, n2 = n2, k = k, k_f = k_f, test = test,
      a_0_d = a_0_d, b_0_d = b_0_d, a_0_a = a_0_a, b_0_a = b_0_a,
      a_1_d = a_1_d, b_1_d = b_1_d, a_2_d = a_2_d, b_2_d = b_2_d,
      a_1_a = a_1_a, b_1_a = b_1_a, a_2_a = a_2_a, b_2_a = b_2_a,
      output = "numeric",
      a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
      a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
      a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
      a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus,
      compute_freq_t1e = compute_freq_t1e,
      p1_grid = p1_grid, p2_grid = p2_grid,
      p1_power = p1_power, p2_power = p2_power
    )
    
    results <- rbind(results, c(n_tot, res[1:3]))
    if (compute_freq_t1e && length(res) > 3 &&
        "Frequentist_Type1_Error" %in% names(res)) {
      freq_t1e_results <- c(freq_t1e_results, res["Frequentist_Type1_Error"])
    }
    if (!is.null(p1_power) && !is.null(p2_power) &&
        "Frequentist_Power" %in% names(res)) {
      freq_power_results <- c(freq_power_results, res["Frequentist_Power"])
    }
    
    if (progress && interactive()) {
      cat(sprintf(
        "\rSimulating: n_total=%-4d (n1=%-3d, n2=%-3d) | Power=%.3f | Type-I=%.3f | P(CE|H0)=%.3f",
        n_tot, n1, n2, res["Power"], res["Type1_Error"], res["CE_H0"]
      ))
      if (compute_freq_t1e && length(freq_t1e_results) > 0)
        cat(sprintf(" | FreqT1e=%.3f", utils::tail(freq_t1e_results, 1)))
      if (!is.null(p1_power) && !is.null(p2_power) &&
          length(freq_power_results) > 0)
        cat(sprintf(" | FreqPow=%.3f", utils::tail(freq_power_results, 1)))
      utils::flush.console()
    }
  }
  message("\n\nSimulation complete.\n\n")
  
  powervec_k    <- results[,2]
  t1evec_k      <- results[,3]
  pceH0_vec_k_f <- results[,4]
  ns            <- results[,1]
  
  power_met <- any(powervec_k >= power)
  t1e_met   <- any(t1evec_k <= alpha)
  pce_met   <- any(pceH0_vec_k_f >= pce_H0)
  
  power_idx <- which(powervec_k >= power)
  n_power   <- ifelse(power_met, ns[min(power_idx)], NA)
  t1e_idx   <- which(t1evec_k <= alpha)
  n_t1e     <- ifelse(t1e_met, ns[min(t1e_idx)], NA)
  pce_idx   <- which(pceH0_vec_k_f >= pce_H0)
  n_pceH0   <- ifelse(pce_met, ns[min(pce_idx)], NA)
  
  # Frequentist power
  n_freq_power <- NA
  freq_power_met <- FALSE
  if (!is.null(p1_power) && !is.null(p2_power) &&
      length(freq_power_results) > 0) {
    freq_power_idx <- which(freq_power_results >= power)
    freq_power_met <- any(freq_power_results >= power)
    if (freq_power_met) n_freq_power <- ns[min(freq_power_idx)]
  }
  
  # ========================
  # PRINT FINAL SUMMARY
  # ========================
  message(sprintf("SUMMARY for %s:\n", test))
  message(sprintf("  Hypotheses: %s\n", hyp_desc))
  message(sprintf("  k = %.3f, k_f = %.3f\n", k, k_f))
  message(sprintf("  Allocation: alloc1 = %.3f, alloc2 = %.3f\n",
                  alloc1, alloc2))
  message(sprintf("  Target power = %.3f, alpha = %.3f, P(CE|H0) = %.3f\n\n",
                  power, alpha, pce_H0))
  
  if (!power_met)
    message(sprintf("    POWER not reached: max=%.3f at n_total=%d\n",
                    max(powervec_k), max(ns)))
  else
    message(sprintf("    Power >= %.3f achieved at n_total=%d\n",
                    power, n_power))
  
  if (!t1e_met)
    message(sprintf("    Bayesian Type-I error too high: min=%.3f at n_total=%d\n",
                    min(t1evec_k), max(ns)))
  else
    message(sprintf("    Bayesian Type-I error <= %.3f achieved at n_total=%d\n",
                    alpha, n_t1e))
  
  if (!pce_met)
    message(sprintf("    P(CE|H0) not reached: max=%.3f at n_total=%d\n",
                    max(pceH0_vec_k_f), max(ns)))
  else
    message(sprintf("    P(CE|H0) >= %.3f achieved at n_total=%d\n",
                    pce_H0, n_pceH0))
  
  if (compute_freq_t1e && length(freq_t1e_results) > 0) {
    max_freq_t1e <- max(freq_t1e_results)
    freq_t1e_met <- max_freq_t1e <= alpha
    freq_t1e_idx <- which(freq_t1e_results <= alpha)
    n_freq_t1e <- ifelse(freq_t1e_met, ns[min(freq_t1e_idx)], NA)
    if (freq_t1e_met)
      message(sprintf(
        "    Frequentist Type-I error <= %.3f achieved (max(sup)=%.3f)\n",
        alpha, max_freq_t1e
      ))
    else
      message(sprintf(
        "    FREQUENTIST Type-I error TOO HIGH: max(sup)=%.3f > %.3f\n",
        max_freq_t1e, alpha
      ))
  }
  
  if (!is.null(p1_power) && !is.null(p2_power) &&
      length(freq_power_results) > 0) {
    if (freq_power_met)
      message(sprintf(
        "    Frequentist power >= %.3f achieved at n_total=%d (p1=%.3f, p2=%.3f)\n",
        power, n_freq_power, p1_power, p2_power
      ))
    else
      message(sprintf(
        "    Frequentist power not reached: max=%.3f at n_total=%d (p1=%.3f, p2=%.3f)\n",
        max(freq_power_results), max(ns), p1_power, p2_power
      ))
  }
  
  # plotting output prep
  df <- data.frame(
    n = ns,
    power = powervec_k,
    t1e = t1evec_k,
    pceH0 = pceH0_vec_k_f
  )
  if (!is.null(p1_power) && !is.null(p2_power) &&
      length(freq_power_results) > 0) {
    df$freq_power <- freq_power_results
  }
  
  # Annotations
  x_annot <- min(ns) + (max(ns) - min(ns)) * 0.02
  y_base <- 0.9
  y_step <- -0.08
  n_power_text <- ifelse(power_met, sprintf("n=%d", n_power), "not reached")
  n_t1e_text <- ifelse(t1e_met, sprintf("n=%d", n_t1e), "not reached")
  n_freq_power_text <- ifelse(freq_power_met,
                              sprintf("n=%d", n_freq_power), "not reached")
  n_pce_text <- ifelse(pce_met, sprintf("n=%d", n_pceH0), "not reached")
  
  xpos_right <- max(ns) * 0.98
  
  # ========================
  # MIDDLE PLOT
  # ========================
  p1_plot <- ggplot2::ggplot(df, ggplot2::aes(x = n)) +
    ggplot2::geom_line(
      ggplot2::aes(y = power, color = "Bayes Power"),
      linewidth = 0.9
    ) +
    ggplot2::geom_line(
      ggplot2::aes(y = t1e, color = "Type I Error"),
      linewidth = 0.9
    )
  
  if ("freq_power" %in% colnames(df))
    p1_plot <- p1_plot +
    ggplot2::geom_line(
      ggplot2::aes(y = freq_power, color = "Frequentist power"),
      linewidth = 0.9
    )
  
  p1_plot <- p1_plot +
    ggplot2::scale_color_manual(
      values = c(
        "Bayes Power" = "black",
        "Type I Error" = "blue",
        "Frequentist power" = "green"
      )
    ) +
    ggplot2::geom_hline(
      yintercept = power,
      linetype = "dashed",
      color = "black",
      alpha = 0.7
    ) +
    ggplot2::geom_hline(
      yintercept = alpha,
      linetype = "dashed",
      color = "blue",
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = ifelse(power_met, n_power, max(ns)),
      color = "black",
      linewidth = 1.1
    ) +
    ggplot2::geom_vline(
      xintercept = ifelse(t1e_met, n_t1e, max(ns)),
      color = "blue",
      linewidth = 1.1
    )
  if (freq_power_met)
    p1_plot <- p1_plot +
    ggplot2::geom_vline(
      xintercept = n_freq_power,
      color = "green",
      linewidth = 1.1
    )
  
  labels_to_draw <- list(
    list(text = n_power_text, color = "black"),
    list(text = n_t1e_text, color = "blue")
  )
  if (freq_power_met && "freq_power" %in% colnames(df))
    labels_to_draw <- append(
      labels_to_draw,
      list(list(text = n_freq_power_text, color = "green"))
    )
  for (i in seq_along(labels_to_draw)) {
    p1_plot <- p1_plot +
      ggplot2::annotate(
        "text",
        x = x_annot,
        y = y_base + (i - 1) * y_step,
        label = labels_to_draw[[i]]$text,
        color = labels_to_draw[[i]]$color,
        size = 4,
        fontface = "bold",
        hjust = 0
      )
  }
  
  p1_plot <- p1_plot +
    ggplot2::labs(
      title = "Power and Type-I Error Rate",
      y = "Probability"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(
        hjust = 0.5, face = "bold", size = 14
      ),
      legend.position = "top",
      legend.title = ggplot2::element_blank()
    )
  
  # ========================
  # BOTTOM PLOT
  # ========================
  bottom_title_expr <- if (test == "BF+-")
    expression("Probability of compelling evidence for " * H["-"])
  else expression("Probability of compelling evidence for H0")
  
  p2_plot <- ggplot2::ggplot(df, ggplot2::aes(x = n, y = pceH0)) +
    ggplot2::coord_cartesian(ylim = c(0, 1)) +
    ggplot2::geom_line(color = "red", linewidth = 0.9) +
    ggplot2::geom_hline(
      yintercept = pce_H0,
      linetype = "dashed",
      color = "red",
      alpha = 0.7
    ) +
    ggplot2::geom_vline(
      xintercept = ifelse(pce_met, n_pceH0, max(ns)),
      color = "red",
      linewidth = 1.1
    ) +
    ggplot2::annotate(
      "text",
      x = x_annot,
      y = 0.9,
      label = n_pce_text,
      color = "red",
      size = 4,
      fontface = "bold",
      hjust = 0
    ) +
    ggplot2::labs(title = bottom_title_expr, y = "Probability") +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14)
    )
  
  # ========================
  # PRIOR PLOT
  # ========================
  xseq <- seq(0, 1, 0.01)
  
  if (test == "BF01") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a_1_d, b_1_d),
        stats::dbeta(xseq, a_1_a, b_1_a),
        stats::dbeta(xseq, a_2_d, b_2_d),
        stats::dbeta(xseq, a_2_a, b_2_a)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"),
                  each = length(xseq)),
      hypothesis = "H[1]"
    )
  } else if (test == "BF+0") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a_1_d, b_1_d),
        stats::dbeta(xseq, a_1_a, b_1_a),
        stats::dbeta(xseq, a_2_d, b_2_d),
        stats::dbeta(xseq, a_2_a, b_2_a)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"),
                  each = length(xseq)),
      hypothesis = "H[+]"
    )
  } else if (test == "BF-0") {
    prior_df <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a_1_d_Hminus, b_1_d_Hminus),
        stats::dbeta(xseq, a_1_a_Hminus, b_1_a_Hminus),
        stats::dbeta(xseq, a_2_d_Hminus, b_2_d_Hminus),
        stats::dbeta(xseq, a_2_a_Hminus, b_2_a_Hminus)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"),
                  each = length(xseq)),
      hypothesis = "H[-]"
    )
  } else {  # BF+-
    prior_df_Hplus <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a_1_d, b_1_d),
        stats::dbeta(xseq, a_1_a, b_1_a),
        stats::dbeta(xseq, a_2_d, b_2_d),
        stats::dbeta(xseq, a_2_a, b_2_a)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"),
                  each = length(xseq)),
      hypothesis = "H[+]"
    )
    prior_df_Hminus <- data.frame(
      x = rep(xseq, 4),
      density = c(
        stats::dbeta(xseq, a_1_d_Hminus, b_1_d_Hminus),
        stats::dbeta(xseq, a_1_a_Hminus, b_1_a_Hminus),
        stats::dbeta(xseq, a_2_d_Hminus, b_2_d_Hminus),
        stats::dbeta(xseq, a_2_a_Hminus, b_2_a_Hminus)
      ),
      prior_type = rep(
        c("Design", "Analysis", "Design", "Analysis"),
        each = length(xseq)
      ),
      param = rep(c("p[1]", "p[1]", "p[2]", "p[2]"),
                  each = length(xseq)),
      hypothesis = "H[-]"
    )
    prior_df <- rbind(prior_df_Hplus, prior_df_Hminus)
  }
  
  prior_title_expr <- switch(
    test,
    "BF01" = expression(
      "Design and analysis priors under " * H[1] * ": " * p[1] != p[2]
    ),
    "BF+0" = expression(
      "Design and analysis priors under " * H["+"] * ": " * p[2] > p[1]
    ),
    "BF-0" = expression(
      "Design and analysis priors under " * H["-"] * ": " * p[1] > p[2]
    ),
    "BF+-" = expression(
      "Design and analysis priors under " * H["+"] * " and " * H["-"]
    )
  )
  
  p_priors <- ggplot2::ggplot(
    prior_df,
    ggplot2::aes(x = x, y = density, linetype = prior_type)
  ) +
    ggplot2::geom_line(linewidth = 1.1) +
    ggplot2::facet_grid(
      hypothesis ~ param,
      scales = "free_y",
      labeller = ggplot2::labeller(param = ggplot2::label_parsed)
    ) +
    ggplot2::scale_linetype_manual(
      breaks = c("Design", "Analysis"),
      values = c("Design" = "solid", "Analysis" = "dotted"),
      labels = c("Design prior (solid)", "Analysis prior (dotted)")
    ) +
    ggplot2::guides(
      linetype = ggplot2::guide_legend(
        override.aes = list(
          linetype = c("solid", "dotted"),
          color = c("black", "black")
        )
      )
    ) +
    ggplot2::labs(
      title = prior_title_expr,
      linetype = "Prior type",
      y = "Density",
      x = "Probability"
    ) +
    ggplot2::theme_minimal(base_size = 11) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(hjust = 0.5, size = 14),
      legend.position = "bottom",
      strip.text = ggplot2::element_text(size = 12, face = "bold"),
      legend.key.width = grid::unit(1.5, "cm")
    )
  
  prior_height <- ifelse(test == "BF+-", 2.0, 1.4)
  
  if (output == "plot" &&
      !requireNamespace("ggplot2", quietly = TRUE)) {
    message("ggplot2 not available for plotting.")
    output <- "numeric"
  }
  
  if (output == "plot") {
    print(
      p_priors / p1_plot / p2_plot +
        patchwork::plot_layout(heights = c(prior_height, 2.5, 1.3))
    )
  }
  
  invisible(list(
    n_power = n_power,
    n_t1e = n_t1e,
    n_pceH0 = n_pceH0,
    n_freq_power = if (freq_power_met) n_freq_power else NA,
    hypotheses = hyp_desc,
    p1_power = p1_power,
    p2_power = p2_power,
    alloc1 = alloc1,
    alloc2 = alloc2,
    results = data.frame(
      n_total = ns,
      power = powervec_k,
      bayes_t1 = t1evec_k,
      pceH0 = pceH0_vec_k_f
    )
  ))
}





#' Sample size calibration for two-arm binomial Bayes factor designs
#'
#' Backward-compatible wrapper around \code{design_twoarm_onestage_bf()}.
#'
#' @inheritParams powertwoarmbinbf01
#' @param power Target Bayesian power.
#' @param alpha Target Bayesian type-I error.
#' @param pce_H0 Target probability of compelling evidence under \eqn{H_0}.
#' @param nrange Integer vector of length 2 giving the search range for total n.
#' @param n_step Step size for n. Currently only \code{n_step = 1} is supported
#'   in the object-based calibration workflow.
#' @param progress Logical; if \code{TRUE}, print progress to the console.
#' @param output \code{"plot"} or \code{"numeric"}.
#' @param alloc1,alloc2 Fixed randomisation probabilities for arm 1 and arm 2;
#'   must be positive and sum to 1.
#' @param sustain_n Non-negative integer. A candidate total sample size is
#'   considered feasible only if the relevant target constraints hold at that
#'   total sample size and for the next \code{sustain_n} larger total sample
#'   sizes in the search range.
#'
#' @return If \code{output = "numeric"}, returns a
#'   \code{"twoarm_onestage_bf_design"} object. If \code{output = "plot"},
#'   the plot is printed and the design object is returned invisibly.
#' @export
ntwoarmbinbf01 <- function(
    k = 1/3, k_f = 3,
    power = 0.8, alpha = 0.05, pce_H0 = 0.9,
    test = c("BF01", "BF+0", "BF-0", "BF+-"),
    nrange = c(10, 150), n_step = 1,
    progress = TRUE, compute_freq_t1e = FALSE,
    p1_grid = seq(0.01, 0.99, 0.02),
    p2_grid = seq(0.01, 0.99, 0.02),
    p1_power = NULL, p2_power = NULL,
    a_0_d = 1, b_0_d = 1, a_0_a = 1, b_0_a = 1,
    a_1_d = 1, b_1_d = 1, a_2_d = 1, b_2_d = 1,
    a_1_a = 1, b_1_a = 1, a_2_a = 1, b_2_a = 1,
    output = c("plot", "numeric"),
    a_1_d_Hminus = 1, b_1_d_Hminus = 1,
    a_2_d_Hminus = 1, b_2_d_Hminus = 1,
    a_1_a_Hminus = 1, b_1_a_Hminus = 1,
    a_2_a_Hminus = 1, b_2_a_Hminus = 1,
    alloc1 = 0.5, alloc2 = 0.5,
    sustain_n = 10L
) {
  test <- match.arg(test)
  output <- match.arg(output)
  
  if (!identical(n_step, 1)) {
    stop("Currently, 'n_step' must be 1 in ntwoarmbinbf01().", call. = FALSE)
  }
  
  calibration <- if (isTRUE(compute_freq_t1e) && !is.null(p1_power) && !is.null(p2_power)) {
    "full"
  } else if (isTRUE(compute_freq_t1e)) {
    "hybrid"
  } else {
    "Bayesian"
  }
  
  obj <- design_twoarm_onestage_bf(
    n_min = nrange[1],
    n_max = nrange[2],
    k = k,
    k_f = k_f,
    test = test,
    a_0_d = a_0_d, b_0_d = b_0_d,
    a_0_a = a_0_a, b_0_a = b_0_a,
    a_1_d = a_1_d, b_1_d = b_1_d,
    a_2_d = a_2_d, b_2_d = b_2_d,
    a_1_a = a_1_a, b_1_a = b_1_a,
    a_2_a = a_2_a, b_2_a = b_2_a,
    a_1_d_Hminus = a_1_d_Hminus, b_1_d_Hminus = b_1_d_Hminus,
    a_2_d_Hminus = a_2_d_Hminus, b_2_d_Hminus = b_2_d_Hminus,
    a_1_a_Hminus = a_1_a_Hminus, b_1_a_Hminus = b_1_a_Hminus,
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus,
    alloc1 = alloc1,
    alloc2 = alloc2,
    calibration = calibration,
    target_power = power,
    target_type1 = alpha,
    target_ce_h0 = pce_H0,
    target_freq_power = power,
    target_freq_type1 = alpha,
    p1_grid = p1_grid,
    p2_grid = p2_grid,
    p1_power = p1_power,
    p2_power = p2_power,
    sustain_n = sustain_n,
    progress = progress
  )
  
  if (output == "plot") {
    print(plot(obj, type = "old"))
    return(invisible(obj))
  }
  
  obj
}
