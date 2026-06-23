#' Internal calibration routine for one-stage two-arm BF designs
#'
#' @param n_min Integer. Minimum admissible total sample size.
#' @param n_max Integer. Maximum admissible total sample size.
#' @param k Numeric scalar greater than 0. Evidence threshold used for power and
#'   type-I error.
#' @param k_f Numeric scalar greater than 1. Threshold used for CE(H0).
#' @param test Character string, one of \code{"BF01"}, \code{"BF+0"},
#'   \code{"BF-0"}, or \code{"BF+-"}.
#' @param a_0_d,b_0_d,a_0_a,b_0_a Shape parameters for design and analysis priors
#'   under \eqn{H_0}.
#' @param a_1_d,b_1_d,a_2_d,b_2_d Shape parameters for design priors under
#'   \eqn{H_1} or \eqn{H_+}.
#' @param a_1_a,b_1_a,a_2_a,b_2_a Shape parameters for analysis priors under
#'   \eqn{H_1} or \eqn{H_+}.
#' @param a_1_d_Hminus,b_1_d_Hminus,a_2_d_Hminus,b_2_d_Hminus Optional design
#'   priors under \eqn{H_-}.
#' @param a_1_a_Hminus,b_1_a_Hminus,a_2_a_Hminus,b_2_a_Hminus Optional analysis
#'   priors under \eqn{H_-}.
#' @param alloc1,alloc2 Fixed randomisation probabilities for arm 1 and arm 2.
#' @param calibration Character string specifying the calibration mode.
#' @param target_power,target_type1,target_ce_h0,target_freq_power,target_freq_type1
#'   Target operating characteristics.
#' @param p1_grid,p2_grid Grids for supremum frequentist type-I error.
#' @param p1_power,p2_power Optional true proportions for frequentist power.
#' @param power_cushion Non-negative numeric scalar applied to power targets.
#' @param sustain_n Non-negative integer. Rolling feasibility window size.
#' @param progress Logical; if \code{TRUE}, emit progress information.
#'
#' @return A list with feasibility, selected design, operating characteristics,
#'   and full search results.
optimal_onestage_twoarm_bf <- function(
    n_min,
    n_max,
    k = 1/3,
    k_f = 3,
    test = c("BF01", "BF+0", "BF-0", "BF+-"),
    a_0_d = 1, b_0_d = 1,
    a_0_a = 1, b_0_a = 1,
    a_1_d = 1, b_1_d = 1,
    a_2_d = 1, b_2_d = 1,
    a_1_a = 1, b_1_a = 1,
    a_2_a = 1, b_2_a = 1,
    a_1_d_Hminus = 1, b_1_d_Hminus = 1,
    a_2_d_Hminus = 1, b_2_d_Hminus = 1,
    a_1_a_Hminus = 1, b_1_a_Hminus = 1,
    a_2_a_Hminus = 1, b_2_a_Hminus = 1,
    alloc1 = 0.5,
    alloc2 = 0.5,
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    target_power = 0.8,
    target_type1 = 0.05,
    target_ce_h0 = 0,
    target_freq_power = 0.8,
    target_freq_type1 = 0.05,
    p1_grid = seq(0.01, 0.99, 0.02),
    p2_grid = seq(0.01, 0.99, 0.02),
    p1_power = NULL,
    p2_power = NULL,
    power_cushion = 0,
    sustain_n = 10L,
    progress = FALSE
) {
  test <- match.arg(test)
  calibration <- match.arg(calibration)
  sustain_n <- as.integer(sustain_n)
  
  grid_n <- seq.int(from = n_min, to = n_max, by = 1L)
  
  rows <- lapply(seq_along(grid_n), function(ii) {
    n_total <- grid_n[ii]
    n1 <- round(n_total * alloc1)
    n2 <- n_total - n1
    
    oc <- powertwoarmbinbf01(
      n1 = n1,
      n2 = n2,
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
      compute_freq_t1e = calibration %in% c("frequentist", "hybrid", "full"),
      p1_grid = p1_grid,
      p2_grid = p2_grid,
      p1_power = p1_power,
      p2_power = p2_power,
      output = "numeric"
    )
    
    bayes_power <- unname(oc["Power"])
    bayes_type1 <- unname(oc["Type1_Error"])
    ce_h0 <- unname(oc["CE_H0"])
    
    freq_type1 <- if ("Frequentist_Type1_Error" %in% names(oc)) {
      unname(oc["Frequentist_Type1_Error"])
    } else {
      NA_real_
    }
    
    freq_power <- if ("Frequentist_Power" %in% names(oc)) {
      unname(oc["Frequentist_Power"])
    } else {
      NA_real_
    }
    
    bayes_ok <- !is.na(bayes_power) &&
      !is.na(bayes_type1) &&
      bayes_power >= (target_power + power_cushion) &&
      bayes_type1 <= target_type1
    
    freq_ok <- !is.na(freq_power) &&
      !is.na(freq_type1) &&
      freq_power >= (target_freq_power + power_cushion) &&
      freq_type1 <= target_freq_type1
    
    hybrid_ok <- !is.na(bayes_power) &&
      !is.na(freq_type1) &&
      bayes_power >= (target_power + power_cushion) &&
      freq_type1 <= target_freq_type1
    
    feasible_ce <- if (target_ce_h0 > 0) {
      !is.na(ce_h0) && ce_h0 >= target_ce_h0
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
    
    if (isTRUE(progress)) {
      message(
        sprintf(
          "n_total=%d (n1=%d, n2=%d): power=%.3f, type1=%.3f, ce=%.3f",
          n_total, n1, n2, bayes_power, bayes_type1, ce_h0
        )
      )
    }
    
    data.frame(
      n = n_total,
      n_total = n_total,
      n1 = n1,
      n2 = n2,
      power = bayes_power,
      type1 = bayes_type1,
      t1e = bayes_type1,
      bayes_t1 = bayes_type1,
      ce_h0 = ce_h0,
      pceH0 = ce_h0,
      freq_power = freq_power,
      freq_type1 = freq_type1,
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
  
  if (nrow(feasible_rows) == 0L) {
    return(list(
      feasible = FALSE,
      status = "No feasible one-stage two-arm design found.",
      design = c(n_total = NA_integer_, n1 = NA_integer_, n2 = NA_integer_),
      operating_characteristics = NULL,
      search_results = search_results,
      calibration = calibration
    ))
  }
  
  best_row <- feasible_rows[1L, , drop = FALSE]
  
  best_oc <- powertwoarmbinbf01(
    n1 = best_row$n1,
    n2 = best_row$n2,
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
    compute_freq_t1e = calibration %in% c("frequentist", "hybrid", "full"),
    p1_grid = p1_grid,
    p2_grid = p2_grid,
    p1_power = p1_power,
    p2_power = p2_power,
    output = "numeric"
  )
  
  list(
    feasible = TRUE,
    status = "Smallest feasible one-stage two-arm design found.",
    design = c(n_total = best_row$n_total, n1 = best_row$n1, n2 = best_row$n2),
    operating_characteristics = standardize_twoarm_oc(best_oc),
    search_results = search_results,
    calibration = calibration
  )
}




eval_twoarm_onestage_constraints <- function(
    oc,
    calibration,
    target_power,
    target_type1,
    target_ce_h0,
    target_freq_power,
    target_freq_type1,
    power_cushion = 0
) {
  bayes_power <- unname(oc["Power"])
  bayes_type1 <- unname(oc["Type1_Error"])
  ce_h0 <- unname(oc["CE_H0"])
  
  freq_type1 <- if ("Frequentist_Type1_Error" %in% names(oc)) {
    unname(oc["Frequentist_Type1_Error"])
  } else {
    NA_real_
  }
  
  freq_power <- if ("Frequentist_Power" %in% names(oc)) {
    unname(oc["Frequentist_Power"])
  } else {
    NA_real_
  }
  
  feasible_ce <- if (target_ce_h0 > 0) {
    !is.na(ce_h0) && ce_h0 >= target_ce_h0
  } else {
    TRUE
  }
  
  bayes_ok <- !is.na(bayes_power) &&
    !is.na(bayes_type1) &&
    bayes_power >= (target_power + power_cushion) &&
    bayes_type1 <= target_type1
  
  freq_ok <- !is.na(freq_power) &&
    !is.na(freq_type1) &&
    freq_power >= (target_freq_power + power_cushion) &&
    freq_type1 <= target_freq_type1
  
  hybrid_ok <- !is.na(bayes_power) &&
    !is.na(freq_type1) &&
    bayes_power >= (target_power + power_cushion) &&
    freq_type1 <= target_freq_type1
  
  switch(
    calibration,
    Bayesian = bayes_ok && feasible_ce,
    frequentist = freq_ok && feasible_ce,
    hybrid = hybrid_ok && feasible_ce,
    full = bayes_ok && freq_ok && feasible_ce
  )
}