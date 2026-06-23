#' Design or evaluate a one-stage single-arm Bayes factor trial
#'
#' Calibrates or evaluates a one-stage single-arm Bayes factor design for a
#' binary endpoint.
#'
#' The design uses the Bayes factor \eqn{BF_{01}}. Small values of
#' \eqn{BF_{01}} indicate evidence against \eqn{H_0}, so efficacy is concluded
#' when \eqn{BF_{01} \le k}. Large values indicate evidence in favour of
#' \eqn{H_0}, and the optional CE(H0) / PCE(H0) constraint is evaluated using
#' the separate threshold \code{k_ce}.
#'
#' Analysis priors are specified separately under \eqn{H_0} and \eqn{H_1} via
#' \code{a0, b0, a1, b1}. Design priors are specified separately under
#' \eqn{H_0} and \eqn{H_1} via \code{da0, db0, da1, db1}.
#'
#' @param n_min Integer. Minimum admissible sample size.
#' @param n_max Integer. Maximum admissible sample size.
#' @param k Numeric scalar greater than 0. Evidence threshold on the
#' \eqn{BF_{01}} scale for efficacy, used for power and type-I error.
#' @param k_ce Optional numeric scalar greater than 1. Threshold on the
#' \eqn{BF_{01}} scale used for CE(H0) / PCE(H0). Must be supplied when
#' \code{target_ce_h0 > 0}.
#' @param p0 Numeric scalar in \eqn{(0,1)}. Null response probability.
#' @param a0,b0 Positive numeric scalars. Beta analysis-prior parameters under
#' \eqn{H_0}.
#' @param a1,b1 Positive numeric scalars. Beta analysis-prior parameters under
#' \eqn{H_1}.
#' @param dp Optional numeric scalar in \eqn{(0,1)}. Fixed point alternative
#' used for frequentist power calculations under \eqn{H_1}.
#' @param da0,db0 Positive numeric scalars. Beta design-prior parameters under
#' \eqn{H_0}.
#' @param da1,db1 Positive numeric scalars. Beta design-prior parameters under
#' \eqn{H_1}.
#' @param type Character string specifying the Bayes-factor test. One of
#' \code{"point"} or \code{"direction"}.
#' @param calibration Character string specifying the calibration mode. One of
#' \code{"Bayesian"}, \code{"frequentist"}, \code{"hybrid"}, or
#' \code{"full"}.
#' @param target_power Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#' power.
#' @param target_type1 Numeric scalar in \eqn{(0,1)}. Target corrected Bayesian
#' type-I error.
#' @param target_ce_h0 Numeric scalar in \eqn{[0,1)}. Optional lower bound on
#' the corrected Bayesian probability of compelling evidence in favour of
#' \eqn{H_0}.
#' @param target_freq_power Numeric scalar in \eqn{(0,1)}. Target corrected
#' frequentist power at \code{dp}.
#' @param target_freq_type1 Numeric scalar in \eqn{(0,1)}. Target corrected
#' frequentist type-I error at \eqn{p = p_0}.
#' @param algorithm Character string specifying whether the design should be
#' optimized or only evaluated.
#' @param n Optional integer sample size used when \code{algorithm = "manual"}.
#' @param power_cushion Optional additive cushion applied to the power targets
#' in the optimizer.
#' @param sustain_n Non-negative integer. A candidate design is considered
#' feasible only if the relevant operating characteristics satisfy their target
#' constraints at the candidate sample size and for the next \code{sustain_n}
#' larger sample sizes, subject to the search range. This also applies to the
#' CE(H0) constraint when \code{target_ce_h0 > 0}.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{"singlearm_onestage_bf_design"}.
#' @export
design_singlearm_onestage_bf <- function(
    n_min,
    n_max,
    k,
    k_ce = NULL,
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
    algorithm = c("optimal", "manual"),
    n = NULL,
    power_cushion = 0,
    sustain_n = 10L,
    ...
) {
  algorithm <- match.arg(algorithm)
  type <- match.arg(type)
  calibration <- match.arg(calibration)
  
  validate_singlearm_onestage_inputs(
    n_min = n_min,
    n_max = n_max,
    k = k,
    k_ce = k_ce,
    p0 = p0,
    a0 = a0,
    b0 = b0,
    a1 = a1,
    b1 = b1,
    dp = dp,
    da0 = da0,
    db0 = db0,
    da1 = da1,
    db1 = db1,
    target_power = target_power,
    target_type1 = target_type1,
    target_ce_h0 = target_ce_h0,
    target_freq_power = target_freq_power,
    target_freq_type1 = target_freq_type1,
    power_cushion = power_cushion,
    sustain_n = sustain_n,
    calibration = calibration
  )
  
  prior_spec <- normalize_singlearm_onestage_design_prior(
    dp = dp, da0 = da0, db0 = db0, da1 = da1, db1 = db1
  )
  
  if (calibration %in% c("frequentist", "full") && is.na(prior_spec$dp)) {
    stop(
      "'dp' must be supplied for calibration = 'frequentist' or calibration = 'full'.",
      call. = FALSE
    )
  }
  
  k_ce_use <- if (isTRUE(target_ce_h0 > 0)) k_ce else NULL
  
  eval_singlearm_onestage_constraints <- function(
    eval_res,
    calibration,
    target_power,
    target_type1,
    target_ce_h0,
    target_freq_power,
    target_freq_type1,
    power_cushion = 0
  ) {
    ce_ok <- if (target_ce_h0 > 0) {
      !is.null(eval_res$pce0_corr) &&
        !is.na(eval_res$pce0_corr) &&
        eval_res$pce0_corr >= target_ce_h0
    } else {
      TRUE
    }
    
    bayes_ok <- isTRUE(
      !is.na(eval_res$pfineff) &&
        !is.na(eval_res$pfineff0) &&
        eval_res$pfineff >= (target_power + power_cushion) &&
        eval_res$pfineff0 <= target_type1 &&
        ce_ok
    )
    
    freq_ok <- isTRUE(
      !is.na(eval_res$pfineff_freq) &&
        !is.na(eval_res$pfineff_freq0) &&
        eval_res$pfineff_freq >= (target_freq_power + power_cushion) &&
        eval_res$pfineff_freq0 <= target_freq_type1 &&
        ce_ok
    )
    
    hybrid_ok <- isTRUE(
      !is.na(eval_res$pfineff) &&
        !is.na(eval_res$pfineff_freq0) &&
        eval_res$pfineff >= (target_power + power_cushion) &&
        eval_res$pfineff_freq0 <= target_freq_type1 &&
        ce_ok
    )
    
    switch(
      calibration,
      Bayesian = bayes_ok,
      frequentist = freq_ok,
      hybrid = hybrid_ok,
      full = bayes_ok && freq_ok
    )
  }
  
  if (algorithm == "manual") {
    if (is.null(n)) {
      stop(
        "For algorithm = 'manual', 'n' must be supplied.",
        call. = FALSE
      )
    }
    
    if (n < n_min || n > n_max) {
      stop("Require n_min <= n <= n_max.", call. = FALSE)
    }
    
    eval_res <- singlearm_fixed_oc(
      n = n,
      k = k,
      p0 = p0,
      a0 = a0, b0 = b0,
      a1 = a1, b1 = b1,
      da0 = prior_spec$da0,
      db0 = prior_spec$db0,
      da1 = prior_spec$da1,
      db1 = prior_spec$db1,
      dp = prior_spec$dp,
      type = type,
      k_ce = k_ce_use
    )
    
    out <- list(
      call = match.call(),
      mode = "manual",
      status = "Manual one-stage design evaluation.",
      calibration = calibration,
      inputs = list(
        n_min = n_min,
        n_max = n_max,
        k = k,
        k_ce = k_ce_use,
        p0 = p0,
        a0 = a0,
        b0 = b0,
        a1 = a1,
        b1 = b1,
        prior = prior_spec,
        dp = dp,
        da0 = da0,
        db0 = db0,
        da1 = da1,
        db1 = db1,
        type = type,
        calibration = calibration,
        target_power = target_power,
        target_type1 = target_type1,
        target_ce_h0 = target_ce_h0,
        target_freq_power = target_freq_power,
        target_freq_type1 = target_freq_type1,
        power_cushion = power_cushion,
        sustain_n = sustain_n
      ),
      design = c(n = n),
      operating_characteristics = eval_res,
      feasible = eval_singlearm_onestage_constraints(
        eval_res = eval_res,
        calibration = calibration,
        target_power = target_power,
        target_type1 = target_type1,
        target_ce_h0 = target_ce_h0,
        target_freq_power = target_freq_power,
        target_freq_type1 = target_freq_type1,
        power_cushion = power_cushion
      )
    )
    
    class(out) <- c("singlearm_onestage_bf_design", "singlearm_bf_design")
    return(out)
  }
  
  opt_res <- optimal_onestage_singlearm_bf(
    n_min = n_min,
    n_max = n_max,
    k = k,
    p0 = p0,
    a0 = a0,
    b0 = b0,
    a1 = a1,
    b1 = b1,
    dp = prior_spec$dp,
    da0 = prior_spec$da0,
    db0 = prior_spec$db0,
    da1 = prior_spec$da1,
    db1 = prior_spec$db1,
    type = type,
    calibration = calibration,
    target_power = target_power,
    target_type1 = target_type1,
    target_ce_h0 = target_ce_h0,
    target_freq_power = target_freq_power,
    target_freq_type1 = target_freq_type1,
    power_cushion = power_cushion,
    k_ce = k_ce_use,
    sustain_n = sustain_n
  )
  
  out <- list(
    call = match.call(),
    mode = "optimal",
    status = opt_res$status,
    message = opt_res$status,
    feasible = opt_res$feasible,
    calibration = opt_res$calibration,
    inputs = list(
      n_min = n_min,
      n_max = n_max,
      k = k,
      k_ce = k_ce_use,
      p0 = p0,
      a0 = a0,
      b0 = b0,
      a1 = a1,
      b1 = b1,
      prior = prior_spec,
      dp = dp,
      da0 = da0,
      db0 = db0,
      da1 = da1,
      db1 = db1,
      type = type,
      calibration = calibration,
      target_power = target_power,
      target_type1 = target_type1,
      target_ce_h0 = target_ce_h0,
      target_freq_power = target_freq_power,
      target_freq_type1 = target_freq_type1,
      power_cushion = power_cushion,
      sustain_n = sustain_n
    ),
    design = opt_res$design,
    operating_characteristics = opt_res$operating_characteristics,
    search_results = opt_res$search_results
  )
  
  class(out) <- c("singlearm_onestage_bf_design", "singlearm_bf_design")
  out
}

validate_singlearm_onestage_inputs <- function(
    n_min, n_max, k, k_ce, p0, a0, b0, a1, b1, dp, da0, db0, da1, db1,
    target_power, target_type1, target_ce_h0,
    target_freq_power, target_freq_type1,
    power_cushion, sustain_n, calibration
) {
  if (!is.numeric(n_min) || length(n_min) != 1L || is.na(n_min) || n_min < 1) {
    stop("'n_min' must be a positive integer.", call. = FALSE)
  }
  
  if (!is.numeric(n_max) || length(n_max) != 1L || is.na(n_max) || n_max < n_min) {
    stop("'n_max' must be an integer greater than or equal to 'n_min'.", call. = FALSE)
  }
  
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0) {
    stop("'k' must be a numeric value > 0.", call. = FALSE)
  }
  
  if (!is.null(k_ce) && (!is.numeric(k_ce) || length(k_ce) != 1L || is.na(k_ce) || k_ce <= 1)) {
    stop("'k_ce' must be NULL or a numeric value > 1.", call. = FALSE)
  }
  
  if (target_ce_h0 > 0 && is.null(k_ce)) {
    stop("'k_ce' must be supplied when 'target_ce_h0' > 0.", call. = FALSE)
  }
  
  if (!is.numeric(p0) || length(p0) != 1L || is.na(p0) || p0 <= 0 || p0 >= 1) {
    stop("'p0' must be in (0, 1).", call. = FALSE)
  }
  
  if (!is.numeric(a0) || !is.numeric(b0) || a0 <= 0 || b0 <= 0) {
    stop("'a0' and 'b0' must be positive.", call. = FALSE)
  }
  
  if (!is.numeric(a1) || !is.numeric(b1) || a1 <= 0 || b1 <= 0) {
    stop("'a1' and 'b1' must be positive.", call. = FALSE)
  }
  
  if (!is.numeric(da0) || !is.numeric(db0) || da0 <= 0 || db0 <= 0) {
    stop("'da0' and 'db0' must be positive.", call. = FALSE)
  }
  
  if (!is.numeric(da1) || !is.numeric(db1) || da1 <= 0 || db1 <= 0) {
    stop("'da1' and 'db1' must be positive.", call. = FALSE)
  }
  
  if (!is.na(dp) && (!is.numeric(dp) || length(dp) != 1L || dp <= 0 || dp >= 1)) {
    stop("'dp' must be NA or in (0, 1).", call. = FALSE)
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
  
  if (calibration %in% c("frequentist", "full") && is.na(dp)) {
    stop(
      "'dp' must be supplied for calibration = 'frequentist' or calibration = 'full'.",
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

normalize_singlearm_onestage_design_prior <- function(dp, da0, db0, da1, db1) {
  list(
    dp = dp,
    da0 = da0,
    db0 = db0,
    da1 = da1,
    db1 = db1
  )
}