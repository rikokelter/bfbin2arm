#' Internal calibration routine for one-stage single-arm BF designs
#'
#' @param n_min Integer. Minimum admissible sample size in the search grid.
#' @param n_max Integer. Maximum admissible sample size in the search grid.
#' @param k Numeric scalar greater than 0. Evidence threshold on the
#'   \eqn{BF_{01}} scale used for efficacy.
#' @param p0 Numeric scalar in \eqn{(0,1)}. Null response probability.
#' @param a0,b0 Positive numeric scalars. Beta analysis-prior parameters under
#'   \eqn{H_0}.
#' @param a1,b1 Positive numeric scalars. Beta analysis-prior parameters under
#'   \eqn{H_1}.
#' @param dp Optional numeric scalar in \eqn{(0,1)}. Fixed point alternative
#'   used for frequentist power calculations under \eqn{H_1}.
#' @param da0,db0 Positive numeric scalars. Beta design-prior parameters under
#'   \eqn{H_0}.
#' @param da1,db1 Positive numeric scalars. Beta design-prior parameters under
#'   \eqn{H_1}.
#' @param type Character string specifying the Bayes-factor test. One of
#'   \code{"point"} or \code{"direction"}.
#' @param calibration Character string specifying the calibration mode. One of
#'   \code{"Bayesian"}, \code{"frequentist"}, \code{"hybrid"}, or
#'   \code{"full"}.
#' @param target_power Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#'   power.
#' @param target_type1 Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#'   type-I error.
#' @param target_ce_h0 Numeric scalar in \eqn{[0,1)}. Optional lower bound on
#'   the corrected Bayesian probability of compelling evidence in favour of
#'   \eqn{H_0}.
#' @param target_freq_power Numeric scalar in \eqn{(0,1)}. Target corrected
#'   frequentist power at \code{dp}.
#' @param target_freq_type1 Numeric scalar in \eqn{(0,1)}. Target corrected
#'   frequentist type-I error at \eqn{p = p_0}.
#' @param power_cushion Non-negative numeric scalar. Optional additive cushion
#'   applied to the power targets during calibration.
#' @param k_ce Optional numeric scalar greater than 1. Threshold on the
#'   \eqn{BF_{01}} scale used for CE(H0) / PCE(H0) calculations.
#' @param sustain_n Non-negative integer. A candidate sample size is declared
#'   feasible only if the relevant constraints are satisfied at that sample size
#'   and for the next \code{sustain_n} larger sample sizes, subject to the
#'   search range.
#'
#' @return A list with feasibility, selected design, operating characteristics,
#'   and full search results.
optimal_onestage_singlearm_bf <- function(
    n_min,
    n_max,
    k,
    p0,
    a0 = 1, b0 = 1,
    a1 = 1, b1 = 1,
    dp = NA_real_,
    da0 = 1, db0 = 1,
    da1 = 1, db1 = 1,
    type = c("point", "direction"),
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    target_power = 0.8,
    target_type1 = 0.05,
    target_ce_h0 = 0,
    target_freq_power = 0.8,
    target_freq_type1 = 0.05,
    power_cushion = 0,
    k_ce = NULL,
    sustain_n = 10L
) {
  type <- match.arg(type)
  calibration <- match.arg(calibration)
  sustain_n <- as.integer(sustain_n)
  
  grid_n <- seq.int(from = n_min, to = n_max, by = 1L)
  
  rows <- lapply(grid_n, function(n) {
    oc <- singlearm_fixed_oc(
      n = n,
      k = k,
      p0 = p0,
      a0 = a0, b0 = b0,
      a1 = a1, b1 = b1,
      da0 = da0, db0 = db0,
      da1 = da1, db1 = db1,
      dp = dp,
      type = type,
      k_ce = k_ce
    )
    
    bayes_ok <- !is.na(oc$pfineff) &&
      !is.na(oc$pfineff0) &&
      oc$pfineff >= (target_power + power_cushion) &&
      oc$pfineff0 <= target_type1
    
    freq_ok <- !is.na(oc$pfineff_freq) &&
      !is.na(oc$pfineff_freq0) &&
      oc$pfineff_freq >= (target_freq_power + power_cushion) &&
      oc$pfineff_freq0 <= target_freq_type1
    
    hybrid_ok <- !is.na(oc$pfineff) &&
      !is.na(oc$pfineff_freq0) &&
      oc$pfineff >= (target_power + power_cushion) &&
      oc$pfineff_freq0 <= target_freq_type1
    
    feasible_ce <- if (target_ce_h0 > 0) {
      !is.na(oc$pce0_corr) && oc$pce0_corr >= target_ce_h0
    } else {
      TRUE
    }
    
    feasible_pointwise <- switch(
      calibration,
      Bayesian = bayes_ok && feasible_ce,
      frequentist = freq_ok && feasible_ce,
      hybrid = hybrid_ok && feasible_ce,
      full = bayes_ok && freq_ok && feasible_ce
    )
    
    data.frame(
      n = n,
      power = oc$pfineff,
      type1 = oc$pfineff0,
      ce_h0 = oc$pce0_corr,
      power_naive = oc$pnaive,
      type1_naive = oc$pnaive0,
      erased_power = oc$perased,
      erased_type1 = oc$perased0,
      freq_power = oc$pfineff_freq,
      freq_type1 = oc$pfineff_freq0,
      freq_power_naive = oc$pnaive_freq,
      freq_type1_naive = oc$pnaive_freq0,
      erased_freq_power = oc$perased_freq,
      erased_freq_type1 = oc$perased_freq0,
      bayes_ok = bayes_ok,
      freq_ok = freq_ok,
      hybrid_ok = hybrid_ok,
      feasible_ce = feasible_ce,
      feasible_pointwise = feasible_pointwise,
      feasible = FALSE,
      calibration = calibration,
      stringsAsFactors = FALSE
    )
  })
  
  search_results <- do.call(rbind, rows)
  
  sustained_ok <- logical(nrow(search_results))
  
  for (i in seq_len(nrow(search_results))) {
    idx_end <- min(nrow(search_results), i + sustain_n)
    idx <- i:idx_end
    
    sustained_ok[i] <- switch(
      calibration,
      Bayesian = all(search_results$bayes_ok[idx] & search_results$feasible_ce[idx]),
      frequentist = all(search_results$freq_ok[idx] & search_results$feasible_ce[idx]),
      hybrid = all(search_results$hybrid_ok[idx] & search_results$feasible_ce[idx]),
      full = all(
        search_results$bayes_ok[idx] &
          search_results$freq_ok[idx] &
          search_results$feasible_ce[idx]
      )
    )
  }
  
  search_results$feasible <- sustained_ok
  
  feasible_rows <- search_results[search_results$feasible, , drop = FALSE]
  
  if (nrow(feasible_rows) == 0) {
    return(list(
      feasible = FALSE,
      status = "No feasible one-stage design found.",
      design = c(n = NA_integer_),
      operating_characteristics = NULL,
      search_results = search_results,
      calibration = calibration
    ))
  }
  
  best_row <- feasible_rows[1L, , drop = FALSE]
  
  best_oc <- singlearm_fixed_oc(
    n = best_row$n,
    k = k,
    p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    da0 = da0, db0 = db0,
    da1 = da1, db1 = db1,
    dp = dp,
    type = type,
    k_ce = k_ce
  )
  
  list(
    feasible = TRUE,
    status = "Smallest feasible one-stage design found.",
    design = c(n = best_row$n),
    operating_characteristics = best_oc,
    search_results = search_results,
    calibration = calibration
  )
}