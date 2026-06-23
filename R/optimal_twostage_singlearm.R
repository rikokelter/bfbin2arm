#' Optimal two-stage single-arm Bayes factor design
#'
#' Searches over admissible two-stage single-arm designs with a binary endpoint
#' and returns the feasible design with smallest expected sample size under
#' `H0`.
#'
#' Analysis priors are specified separately under `H0` and `H1` via
#' `a0, b0, a1, b1`. Design priors are specified separately under `H0` and `H1`
#' via `da0, db0, da1, db1`.
#'
#' @param n1_min Minimum admissible interim sample size.
#' @param n2_max Maximum admissible final sample size.
#' @param k Efficacy threshold on the `BF01` scale.
#' @param k_f Futility threshold on the `BF01` scale.
#' @param p0 Null response probability.
#' @param a0,b0 Beta analysis-prior parameters under `H0`.
#' @param a1,b1 Beta analysis-prior parameters under `H1`.
#' @param dp Optional fixed point alternative used for frequentist power.
#' @param da0,db0 Beta design-prior parameters under `H0`.
#' @param da1,db1 Beta design-prior parameters under `H1`.
#' @param type Character string; one of `"point"` or `"direction"`.
#' @param calibration Character string; one of `"Bayesian"`,
#'   `"frequentist"`, `"hybrid"`, or `"full"`.
#' @param target_power Target corrected Bayesian power.
#' @param target_type1 Target corrected Bayesian type-I error.
#' @param target_ce_h0 Optional lower bound on corrected Bayesian compelling
#'   evidence in favour of `H0`.
#' @param target_freq_power Target corrected frequentist power at `dp`.
#' @param target_freq_type1 Target corrected frequentist type-I error at `p0`.
#' @param power_cushion Optional additive cushion for the fixed-sample power
#'   target in the first step of the search.
#'
#' @return A list describing the optimal design and search results.
#' @export
optimal_twostage_singlearm_bf <- function(
    n1_min,
    n2_max,
    k,
    k_f,
    p0,
    a0 = 1,
    b0 = 1,
    a1 = 1,
    b1 = 1,
    dp = NA_real_,
    da0 = 1,
    db0 = 1,
    da1 = 1,
    db1 = 1,
    type = c("point", "direction"),
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    target_power = 0.8,
    target_type1 = 0.05,
    target_ce_h0 = 0,
    target_freq_power = 0.8,
    target_freq_type1 = 0.05,
    power_cushion = 0
) {
  type <- match.arg(type)
  calibration <- match.arg(calibration)
  
  k_ce_use <- if (isTRUE(target_ce_h0 > 0)) k_f else NULL
  
  fixed_ok <- function(n) {
    tmp <- singlearm_fixed_oc(
      n = n,
      k = k,
      p0 = p0,
      a0 = a0, b0 = b0,
      a1 = a1, b1 = b1,
      da0 = da0, db0 = db0,
      da1 = da1, db1 = db1,
      dp = dp,
      type = type,
      k_ce = k_ce_use
    )
    
    bayes_ok <- isTRUE(
      !is.na(tmp$pfineff) &&
        !is.na(tmp$pfineff0) &&
        tmp$pfineff >= (target_power + power_cushion) &&
        tmp$pfineff0 <= target_type1
    )
    
    freq_ok <- isTRUE(
      !is.na(tmp$pfineff_freq) &&
        !is.na(tmp$pfineff_freq0) &&
        tmp$pfineff_freq >= (target_freq_power + power_cushion) &&
        tmp$pfineff_freq0 <= target_freq_type1
    )
    
    hybrid_ok <- isTRUE(
      !is.na(tmp$pfineff) &&
        !is.na(tmp$pfineff_freq0) &&
        tmp$pfineff >= (target_power + power_cushion) &&
        tmp$pfineff_freq0 <= target_freq_type1
    )
    
    feasible_ce <- if (target_ce_h0 > 0) {
      !is.na(tmp$pce0_corr) && tmp$pce0_corr >= target_ce_h0
    } else {
      TRUE
    }
    
    switch(
      calibration,
      Bayesian   = bayes_ok && feasible_ce,
      frequentist = freq_ok,
      hybrid     = hybrid_ok,
      full       = bayes_ok && freq_ok && feasible_ce
    )
  }
  
  n2_anchor <- NA_integer_
  for (n in seq.int(from = n1_min + 1L, to = n2_max, by = 1L)) {
    if (fixed_ok(n)) {
      n2_anchor <- n
      break
    }
  }
  
  if (is.na(n2_anchor)) {
    return(list(
      feasible = FALSE,
      status = "No feasible fixed-sample anchor found.",
      design = c(n1 = NA_integer_, n2 = NA_integer_),
      operating_characteristics = NULL,
      search_results = NULL,
      fixed_anchor = NA_integer_,
      calibration = calibration
    ))
  }
  
  grid_n1 <- seq.int(from = n1_min, to = n2_anchor - 1L, by = 1L)
  
  rows <- lapply(grid_n1, function(n1) {
    oc <- powerbinbf01seq(
      n1 = n1,
      n2 = n2_anchor,
      k = k,
      kf = k_f,
      p0 = p0,
      a0 = a0, b0 = b0,
      a1 = a1, b1 = b1,
      da0 = da0, db0 = db0,
      da1 = da1, db1 = db1,
      dp = dp,
      type = type,
      k_ce = k_ce_use
    )
    
    bayes_ok <- !is.na(oc$pfineff) &&
      !is.na(oc$pfineff0) &&
      oc$pfineff >= target_power &&
      oc$pfineff0 <= target_type1
    
    freq_ok <- !is.na(oc$pfineff_freq) &&
      !is.na(oc$pfineff_freq0) &&
      oc$pfineff_freq >= target_freq_power &&
      oc$pfineff_freq0 <= target_freq_type1
    
    hybrid_ok <- !is.na(oc$pfineff) &&
      !is.na(oc$pfineff_freq0) &&
      oc$pfineff >= target_power &&
      oc$pfineff_freq0 <= target_freq_type1
    
    feasible_ce <- if (target_ce_h0 > 0) {
      !is.na(oc$pce0_corr) && oc$pce0_corr >= target_ce_h0
    } else {
      TRUE
    }
    
    feasible <- switch(
      calibration,
      Bayesian = bayes_ok && feasible_ce,
      frequentist = freq_ok,
      hybrid = hybrid_ok,
      full = bayes_ok && freq_ok && feasible_ce
    )
    
    data.frame(
      n1 = n1,
      n2 = n2_anchor,
      power = oc$pfineff,
      type1 = oc$pfineff0,
      ce_h0 = oc$pce0_corr,
      en_h0 = oc$nexp0,
      en_h1 = oc$nexp,
      power_naive = oc$pnaive,
      type1_naive = oc$pnaive0,
      erased_power = oc$perased,
      erased_type1 = oc$perased0,
      pintfut = oc$pintfut,
      pintfut0 = oc$pintfut0,
      freq_power = oc$pfineff_freq,
      freq_type1 = oc$pfineff_freq0,
      freq_en_h1 = oc$nexp_freq,
      freq_en_h0 = oc$nexp_freq0,
      freq_power_naive = oc$pnaive_freq,
      freq_type1_naive = oc$pnaive_freq0,
      erased_freq_power = oc$perased_freq,
      erased_freq_type1 = oc$perased_freq0,
      pintfut_freq = oc$pintfut_freq,
      pintfut_freq0 = oc$pintfut_freq0,
      freq0_argmax = oc$freq0_argmax,
      bayes_ok = bayes_ok,
      freq_ok = freq_ok,
      hybrid_ok = hybrid_ok,
      feasible_ce = feasible_ce,
      feasible = feasible,
      calibration = calibration,
      stringsAsFactors = FALSE
    )
  })
  
  search_results <- do.call(rbind, rows)
  feasible_rows <- search_results[search_results$feasible, , drop = FALSE]
  
  if (nrow(feasible_rows) == 0) {
    return(list(
      feasible = FALSE,
      status = "No feasible two-stage design found.",
      design = c(n1 = NA_integer_, n2 = n2_anchor),
      operating_characteristics = NULL,
      search_results = search_results,
      fixed_anchor = n2_anchor,
      calibration = calibration
    ))
  }
  
  best_idx <- which.min(feasible_rows$en_h0)
  best_row <- feasible_rows[best_idx, , drop = FALSE]
  
  best_oc <- powerbinbf01seq(
    n1 = best_row$n1,
    n2 = best_row$n2,
    k = k,
    kf = k_f,
    p0 = p0,
    a0 = a0, b0 = b0,
    a1 = a1, b1 = b1,
    da0 = da0, db0 = db0,
    da1 = da1, db1 = db1,
    dp = dp,
    type = type,
    k_ce = k_ce_use
  )
  
  list(
    feasible = TRUE,
    status = "Optimal feasible design found.",
    design = c(n1 = best_row$n1, n2 = best_row$n2),
    operating_characteristics = standardize_singlearm_oc(best_oc),
    search_results = search_results,
    fixed_anchor = n2_anchor,
    calibration = calibration
  )
}