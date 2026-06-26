# =============================================================================
# Utility helpers for single-arm two-stage ROPE designs
#
# All computations reuse the shared helpers already defined in utils_rope.R:
#   beta_binom_pmf_rope(y, n, a, b)
#   posterior_rope_prob(y, n, p0, delta, analysis_prior)
#   .validate_beta_prior(), .validate_probability(), .validate_count()
# =============================================================================

# -----------------------------------------------------------------------------
#' Continuation region for stage 1
#'
#' Returns y1 values for which the interim posterior ROPE probability exceeds
#' the futility threshold, so the trial proceeds to stage 2.
#'
#' @keywords internal
#' @noRd
.continuation_region_twostage <- function(n1, p0, delta, gamma_1,
                                          analysis_prior) {
  y1_vals <- 0:n1
  post <- vapply(
    y1_vals,
    posterior_rope_prob,
    numeric(1L),
    n              = n1,
    p0             = p0,
    delta          = delta,
    analysis_prior = analysis_prior
  )
  y1_vals[post > gamma_1]
}

# -----------------------------------------------------------------------------
#' Predictive probability of declaring equivalence — two-stage design
#'
#' @keywords internal
#' @noRd
.predictive_equiv_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_eq, analysis_prior, design_prior
) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  if (length(cont) == 0L) return(0)
  n  <- n1 + n2
  aD <- design_prior[1]; bD <- design_prior[2]
  
  out <- 0
  for (y1 in cont) {
    p_y1    <- beta_binom_pmf_rope(y1, n1, aD, bD)
    dp_post <- c(aD + y1, bD + n1 - y1)
    for (y2 in 0:n2) {
      y         <- y1 + y2
      p_y2_y1   <- beta_binom_pmf_rope(y2, n2, dp_post[1], dp_post[2])
      post_rope <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      if (post_rope >= gamma_eq)
        out <- out + p_y1 * p_y2_y1
    }
  }
  out
}

# -----------------------------------------------------------------------------
#' Predictive probability of compelling evidence against equivalence — two-stage
#'
#' @keywords internal
#' @noRd
.predictive_pce_twostage <- function(
    n1, n2, p0, delta, gamma_1, gamma_diff, analysis_prior, design_prior
) {
  cont <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  if (length(cont) == 0L) return(0)
  n  <- n1 + n2
  aD <- design_prior[1]; bD <- design_prior[2]
  
  out <- 0
  for (y1 in cont) {
    p_y1    <- beta_binom_pmf_rope(y1, n1, aD, bD)
    dp_post <- c(aD + y1, bD + n1 - y1)
    for (y2 in 0:n2) {
      y         <- y1 + y2
      p_y2_y1   <- beta_binom_pmf_rope(y2, n2, dp_post[1], dp_post[2])
      post_rope <- posterior_rope_prob(y, n, p0, delta, analysis_prior)
      if ((1 - post_rope) >= gamma_diff)
        out <- out + p_y1 * p_y2_y1
    }
  }
  out
}

# -----------------------------------------------------------------------------
#' Expected sample size under a design prior — two-stage design
#'
#' @keywords internal
#' @noRd
.expected_n_twostage <- function(n1, n2, p0, delta, gamma_1,
                                 analysis_prior, design_prior) {
  cont   <- .continuation_region_twostage(n1, p0, delta, gamma_1, analysis_prior)
  aD     <- design_prior[1]; bD <- design_prior[2]
  p_cont <- if (length(cont) == 0L) 0 else
    sum(vapply(cont, beta_binom_pmf_rope, numeric(1L), n = n1, a = aD, b = bD))
  n1 + n2 * p_cont
}

# -----------------------------------------------------------------------------
#' One-stage predictive probability of declaring equivalence
#'
#' @keywords internal
#' @noRd
.predictive_equiv_onestage <- function(n, p0, delta, gamma_eq,
                                       analysis_prior, design_prior) {
  y_vals    <- 0:n
  post_rope <- vapply(y_vals, posterior_rope_prob, numeric(1L),
                      n = n, p0 = p0, delta = delta,
                      analysis_prior = analysis_prior)
  aD <- design_prior[1]; bD <- design_prior[2]
  sum(vapply(y_vals[post_rope >= gamma_eq], beta_binom_pmf_rope, numeric(1L),
             n = n, a = aD, b = bD))
}

# -----------------------------------------------------------------------------
#' One-stage predictive probability of compelling evidence against equivalence
#'
#' @keywords internal
#' @noRd
.predictive_pce_onestage <- function(n, p0, delta, gamma_diff,
                                     analysis_prior, design_prior) {
  y_vals    <- 0:n
  post_rope <- vapply(y_vals, posterior_rope_prob, numeric(1L),
                      n = n, p0 = p0, delta = delta,
                      analysis_prior = analysis_prior)
  aD <- design_prior[1]; bD <- design_prior[2]
  sum(vapply(y_vals[(1 - post_rope) >= gamma_diff], beta_binom_pmf_rope,
             numeric(1L), n = n, a = aD, b = bD))
}

# -----------------------------------------------------------------------------
#' Evaluate all operating characteristics for one candidate two-stage design
#'
#' @keywords internal
#' @noRd
.evaluate_rope_twostage_oc <- function(n1, n2, p0, delta, gamma_1, gamma_eq,
                                       gamma_diff, analysis_prior,
                                       design_prior_h0, design_prior_h1) {
  n <- n1 + n2
  list(
    n1        = n1,
    n2        = n2,
    n         = n,
    type1_1st = .predictive_equiv_onestage(n, p0, delta, gamma_eq,
                                           analysis_prior, design_prior_h0),
    power_1st = .predictive_equiv_onestage(n, p0, delta, gamma_eq,
                                           analysis_prior, design_prior_h1),
    pce_1st   = .predictive_pce_onestage(n, p0, delta, gamma_diff,
                                         analysis_prior, design_prior_h0),
    type1_2st = .predictive_equiv_twostage(n1, n2, p0, delta, gamma_1, gamma_eq,
                                           analysis_prior, design_prior_h0),
    power_2st = .predictive_equiv_twostage(n1, n2, p0, delta, gamma_1, gamma_eq,
                                           analysis_prior, design_prior_h1),
    pce_2st   = .predictive_pce_twostage(n1, n2, p0, delta, gamma_1, gamma_diff,
                                         analysis_prior, design_prior_h0),
    EN0       = .expected_n_twostage(n1, n2, p0, delta, gamma_1,
                                     analysis_prior, design_prior_h0),
    EN1       = .expected_n_twostage(n1, n2, p0, delta, gamma_1,
                                     analysis_prior, design_prior_h1),
    cont_region = .continuation_region_twostage(n1, p0, delta, gamma_1,
                                                analysis_prior)
  )
}