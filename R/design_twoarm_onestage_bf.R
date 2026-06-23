#' Design or evaluate a one-stage two-arm Bayes factor trial
#'
#' Calibrates or evaluates a one-stage two-arm Bayes factor design for a binary
#' endpoint with fixed randomisation between the two arms.
#'
#' The design uses one of the Bayes factor tests implemented in
#' \code{powertwoarmbinbf01()}. Small values of the relevant inverted Bayes
#' factor indicate evidence against the null, so efficacy is concluded when the
#' Bayes factor is below \code{k}. Large values indicate evidence in favour of
#' the null (or \eqn{H_-} for \code{test = "BF+-"}), and the optional
#' CE(H0) / PCE(H0) constraint is evaluated using \code{k_f}.
#'
#' @param n_min Integer. Minimum admissible total sample size.
#' @param n_max Integer. Maximum admissible total sample size.
#' @param k Numeric scalar greater than 0. Evidence threshold used for power and
#' type-I error.
#' @param k_f Numeric scalar greater than 1. Threshold used for CE(H0) /
#' PCE(H0).
#' @param test Character string, one of \code{"BF01"}, \code{"BF+0"},
#' \code{"BF-0"}, or \code{"BF+-"}.
#' @param a_0_d,b_0_d,a_0_a,b_0_a Shape parameters for design and analysis priors
#' under \eqn{H_0}.
#' @param a_1_d,b_1_d,a_2_d,b_2_d Shape parameters for design priors under
#' \eqn{H_1} or \eqn{H_+}.
#' @param a_1_a,b_1_a,a_2_a,b_2_a Shape parameters for analysis priors under
#' \eqn{H_1} or \eqn{H_+}.
#' @param a_1_d_Hminus,b_1_d_Hminus,a_2_d_Hminus,b_2_d_Hminus Optional design
#' priors under \eqn{H_-} for directional tests.
#' @param a_1_a_Hminus,b_1_a_Hminus,a_2_a_Hminus,b_2_a_Hminus Optional analysis
#' priors under \eqn{H_-} for directional tests.
#' @param alloc1,alloc2 Fixed randomisation probabilities for arm 1 and arm 2.
#' Must be positive and sum to 1.
#' @param calibration Character string specifying the calibration mode. One of
#' \code{"Bayesian"}, \code{"frequentist"}, \code{"hybrid"}, or
#' \code{"full"}.
#' @param target_power Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#' power.
#' @param target_type1 Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#' type-I error.
#' @param target_ce_h0 Numeric scalar in \eqn{[0,1)}. Optional lower bound on
#' the corrected Bayesian probability of compelling evidence in favour of
#' \eqn{H_0} (or \eqn{H_-} for \code{test = "BF+-"}).
#' @param target_freq_power Numeric scalar in \eqn{(0,1)}. Target frequentist
#' power under \code{p1_power, p2_power}.
#' @param target_freq_type1 Numeric scalar in \eqn{(0,1)}. Target frequentist
#' type-I error.
#' @param p1_grid,p2_grid Grids of true proportions used to compute supremum
#' frequentist type-I error.
#' @param p1_power,p2_power Optional true proportions used for frequentist power.
#' @param power_cushion Non-negative numeric scalar. Optional additive cushion
#' applied to the power targets during calibration.
#' @param sustain_n Non-negative integer. A candidate total sample size is
#' considered feasible only if the relevant target constraints hold at that
#' total sample size and for the next \code{sustain_n} larger total sample
#' sizes in the search range.
#' @param report_freq_type1 Logical. If \code{TRUE}, compute and report the
#' frequentist type-I error for the final selected design even when the chosen
#' calibration mode does not use frequentist criteria. This additional
#' computation has no effect on the calibration itself. Defaults to
#' \code{FALSE}.
#' @param algorithm Character string specifying whether the design should be
#' optimized or only evaluated.
#' @param n_total Optional integer total sample size used when
#' \code{algorithm = "manual"}.
#' @param progress Logical; if \code{TRUE}, print simple progress information
#' during optimization.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{"twoarm_onestage_bf_design"}.
#' @export
design_twoarm_onestage_bf <- function(
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
    report_freq_type1 = FALSE,
    algorithm = c("optimal", "manual"),
    n_total = NULL,
    progress = FALSE,
    ...
) {
  test <- match.arg(test)
  calibration <- match.arg(calibration)
  algorithm <- match.arg(algorithm)
  
  validate_twoarm_onestage_inputs(
    n_min = n_min,
    n_max = n_max,
    k = k,
    k_f = k_f,
    alloc1 = alloc1,
    alloc2 = alloc2,
    target_power = target_power,
    target_type1 = target_type1,
    target_ce_h0 = target_ce_h0,
    target_freq_power = target_freq_power,
    target_freq_type1 = target_freq_type1,
    power_cushion = power_cushion,
    sustain_n = sustain_n,
    calibration = calibration,
    p1_power = p1_power,
    p2_power = p2_power,
    report_freq_type1 = report_freq_type1
  )
  
  need_freq_for_calibration <- calibration %in% c("frequentist", "hybrid", "full")
  need_freq_for_reporting <- isTRUE(report_freq_type1)
  
  if (algorithm == "manual") {
    if (is.null(n_total)) {
      stop("For algorithm = 'manual', 'n_total' must be supplied.", call. = FALSE)
    }
    
    if (n_total < n_min || n_total > n_max) {
      stop("Require n_min <= n_total <= n_max.", call. = FALSE)
    }
    
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
      compute_freq_t1e = need_freq_for_calibration || need_freq_for_reporting,
      p1_grid = p1_grid,
      p2_grid = p2_grid,
      p1_power = p1_power,
      p2_power = p2_power,
      output = "numeric"
    )
    
    feasible <- eval_twoarm_onestage_constraints(
      oc = oc,
      calibration = calibration,
      target_power = target_power,
      target_type1 = target_type1,
      target_ce_h0 = target_ce_h0,
      target_freq_power = target_freq_power,
      target_freq_type1 = target_freq_type1,
      power_cushion = power_cushion
    )
    
    out <- list(
      call = match.call(),
      mode = "manual",
      status = "Manual one-stage two-arm design evaluation.",
      feasible = feasible,
      calibration = calibration,
      inputs = list(
        n_min = n_min,
        n_max = n_max,
        k = k,
        k_f = k_f,
        test = test,
        alloc1 = alloc1,
        alloc2 = alloc2,
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
        calibration = calibration,
        target_power = target_power,
        target_type1 = target_type1,
        target_ce_h0 = target_ce_h0,
        target_freq_power = target_freq_power,
        target_freq_type1 = target_freq_type1,
        p1_grid = p1_grid,
        p2_grid = p2_grid,
        p1_power = p1_power,
        p2_power = p2_power,
        power_cushion = power_cushion,
        sustain_n = sustain_n,
        report_freq_type1 = report_freq_type1
      ),
      design = c(n_total = n_total, n1 = n1, n2 = n2),
      operating_characteristics = standardize_twoarm_oc(oc),
      search_results = NULL
    )
    
    class(out) <- "twoarm_onestage_bf_design"
    return(out)
  }
  
  opt_res <- optimal_onestage_twoarm_bf(
    n_min = n_min,
    n_max = n_max,
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
    target_power = target_power,
    target_type1 = target_type1,
    target_ce_h0 = target_ce_h0,
    target_freq_power = target_freq_power,
    target_freq_type1 = target_freq_type1,
    p1_grid = p1_grid,
    p2_grid = p2_grid,
    p1_power = p1_power,
    p2_power = p2_power,
    power_cushion = power_cushion,
    sustain_n = sustain_n,
    progress = progress
  )
  
  final_oc <- opt_res$operating_characteristics
  
  if (isTRUE(opt_res$feasible) &&
      !is.null(opt_res$design) &&
      isTRUE(report_freq_type1) &&
      calibration == "Bayesian") {
    
    n1_final <- unname(opt_res$design["n1"])
    n2_final <- unname(opt_res$design["n2"])
    
    oc_final_report <- powertwoarmbinbf01(
      n1 = n1_final,
      n2 = n2_final,
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
      compute_freq_t1e = TRUE,
      p1_grid = p1_grid,
      p2_grid = p2_grid,
      p1_power = p1_power,
      p2_power = p2_power,
      output = "numeric"
    )
    
    final_oc <- standardize_twoarm_oc(oc_final_report)
  }
  
  out <- list(
    call = match.call(),
    mode = "optimal",
    status = opt_res$status,
    feasible = opt_res$feasible,
    calibration = opt_res$calibration,
    inputs = list(
      n_min = n_min,
      n_max = n_max,
      k = k,
      k_f = k_f,
      test = test,
      alloc1 = alloc1,
      alloc2 = alloc2,
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
      calibration = calibration,
      target_power = target_power,
      target_type1 = target_type1,
      target_ce_h0 = target_ce_h0,
      target_freq_power = target_freq_power,
      target_freq_type1 = target_freq_type1,
      p1_grid = p1_grid,
      p2_grid = p2_grid,
      p1_power = p1_power,
      p2_power = p2_power,
      power_cushion = power_cushion,
      sustain_n = sustain_n,
      report_freq_type1 = report_freq_type1
    ),
    design = opt_res$design,
    operating_characteristics = final_oc,
    search_results = opt_res$search_results
  )
  
  class(out) <- "twoarm_onestage_bf_design"
  out
}

validate_twoarm_onestage_inputs <- function(
    n_min, n_max, k, k_f, alloc1, alloc2,
    target_power, target_type1, target_ce_h0,
    target_freq_power, target_freq_type1,
    power_cushion, sustain_n, calibration,
    p1_power, p2_power, report_freq_type1
) {
  if (!is.numeric(n_min) || length(n_min) != 1L || is.na(n_min) || n_min < 2) {
    stop("'n_min' must be a positive integer >= 2.", call. = FALSE)
  }
  
  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) || n_max < n_min) {
    stop("'n_max' must be an integer greater than or equal to 'n_min'.", call. = FALSE)
  }
  
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0) {
    stop("'k' must be a numeric value > 0.", call. = FALSE)
  }
  
  if (!is.numeric(k_f) || length(k_f) != 1L || is.na(k_f) || k_f <= 1) {
    stop("'k_f' must be a numeric value > 1.", call. = FALSE)
  }
  
  if (!is.numeric(alloc1) || !is.numeric(alloc2) ||
      length(alloc1) != 1L || length(alloc2) != 1L ||
      is.na(alloc1) || is.na(alloc2) ||
      alloc1 <= 0 || alloc2 <= 0) {
    stop("'alloc1' and 'alloc2' must be positive scalars.", call. = FALSE)
  }
  
  if (!isTRUE(all.equal(alloc1 + alloc2, 1, tolerance = 1e-8))) {
    stop("'alloc1' and 'alloc2' must sum to 1.", call. = FALSE)
  }
  
  if (target_power <= 0 || target_power >= 1) {
    stop("'target_power' must be in (0, 1).", call. = FALSE)
  }
  
  if (target_type1 <= 0 || target_type1 >= 1) {
    stop("'target_type1' must be in (0, 1).", call. = FALSE)
  }
  
  if (target_ce_h0 < 0 || target_ce_h0 >= 1) {
    stop("'target_ce_h0' must be in [0, 1).", call. = FALSE)
  }
  
  if (target_freq_power <= 0 || target_freq_power >= 1) {
    stop("'target_freq_power' must be in (0, 1).", call. = FALSE)
  }
  
  if (target_freq_type1 <= 0 || target_freq_type1 >= 1) {
    stop("'target_freq_type1' must be in (0, 1).", call. = FALSE)
  }
  
  if (power_cushion < 0) {
    stop("'power_cushion' must be non-negative.", call. = FALSE)
  }
  
  if (!is.numeric(sustain_n) || length(sustain_n) != 1L || is.na(sustain_n) ||
      sustain_n < 0 || abs(sustain_n - round(sustain_n)) > sqrt(.Machine$double.eps)) {
    stop("'sustain_n' must be a single non-negative integer.", call. = FALSE)
  }
  
  if (!is.character(calibration) || length(calibration) != 1L ||
      !calibration %in% c("Bayesian", "frequentist", "hybrid", "full")) {
    stop(
      "'calibration' must be one of 'Bayesian', 'frequentist', 'hybrid', or 'full'.",
      call. = FALSE
    )
  }
  
  if (!is.logical(report_freq_type1) || length(report_freq_type1) != 1L || is.na(report_freq_type1)) {
    stop("'report_freq_type1' must be TRUE or FALSE.", call. = FALSE)
  }
  
  if (calibration %in% c("frequentist", "full")) {
    if (is.null(p1_power) || is.null(p2_power)) {
      stop(
        "'p1_power' and 'p2_power' must be supplied for calibration = 'frequentist' or 'full'.",
        call. = FALSE
      )
    }
  }
  
  invisible(TRUE)
}

standardize_twoarm_oc <- function(x) {
  if (is.null(x)) return(NULL)
  
  list(
    power = unname(x["Power"]),
    type1 = unname(x["Type1_Error"]),
    ce_h0 = unname(x["CE_H0"]),
    freq_type1 = if ("Frequentist_Type1_Error" %in% names(x)) {
      unname(x["Frequentist_Type1_Error"])
    } else {
      NA_real_
    },
    freq_power = if ("Frequentist_Power" %in% names(x)) {
      unname(x["Frequentist_Power"])
    } else {
      NA_real_
    }
  )
}