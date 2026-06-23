#' Design or evaluate a single-arm two-stage Bayes factor trial
#'
#' Calibrates or evaluates a single-arm two-stage Bayes factor design for a
#' binary endpoint with one interim analysis for futility.
#'
#' The design uses the Bayes factor \eqn{BF_{01}}. Small values of
#' \eqn{BF_{01}} indicate evidence against \eqn{H_0}, so final efficacy is
#' concluded when \eqn{BF_{01} \\le k}. Large values indicate evidence in favour
#' of \eqn{H_0}, so interim futility is concluded when
#' \eqn{BF_{01} \\ge k_f}.
#'
#' Analysis priors are specified separately under \eqn{H_0} and \eqn{H_1} via
#' \code{a0, b0, a1, b1}. Design priors are specified separately under
#' \eqn{H_0} and \eqn{H_1} via \code{da0, db0, da1, db1}.
#'
#' @param n1_min Integer. Minimum admissible interim sample size.
#' @param n2_max Integer. Maximum admissible final sample size.
#' @param k Numeric scalar greater than 0. Efficacy threshold on the
#'   \eqn{BF_{01}} scale.
#' @param k_f Numeric scalar greater than 1. Futility threshold on the
#'   \eqn{BF_{01}} scale.
#' @param p0 Numeric scalar in \eqn{(0,1)}. Null response probability.
#' @param a0,b0 Positive numeric scalars. Beta analysis-prior parameters under
#'   \eqn{H_0}.
#' @param a1,b1 Positive numeric scalars. Beta analysis-prior parameters under
#'   \eqn{H_1}.
#' @param da0,db0 Positive numeric scalars. Beta design-prior parameters under
#'   \eqn{H_0}.
#' @param da1,db1 Positive numeric scalars. Beta design-prior parameters under
#'   \eqn{H_1}.
#' @param dp Optional numeric scalar in \eqn{(0,1)}. Fixed point alternative
#'   used for frequentist power calculations under \eqn{H_1}.
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
#' @param algorithm Character string specifying whether the design should be
#'   optimized or only evaluated.
#' @param interim Optional integer interim sample size used when
#'   \code{algorithm = "manual"}.
#' @param final Optional integer final sample size used when
#'   \code{algorithm = "manual"}.
#' @param power_cushion Optional additive cushion applied only in the
#'   fixed-sample anchor search of the first optimization step. This can be
#'   useful because introducing an interim futility analysis typically reduces
#'   corrected power relative to the fixed-sample anchor.
#' @param ... Reserved for future extensions.
#'
#' @return An object of class \code{"singlearm_bf_design"}.
#' @export
design_singlearm_bf <- function(
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
    algorithm = c("optimal", "manual"),
    interim = NULL,
    final = NULL,
    power_cushion = 0,
    ...
) {
  algorithm <- match.arg(algorithm)
  type <- match.arg(type)
  calibration <- match.arg(calibration)
  
  validate_singlearm_seq_inputs(
    n1_min = n1_min,
    n2_max = n2_max,
    k = k,
    k_f = k_f,
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
    calibration = calibration
  )
  
  prior_spec <- normalize_singlearm_design_prior(
    dp = dp, da0 = da0, db0 = db0, da1 = da1, db1 = db1
  )
  
  if (calibration %in% c("frequentist", "full") && is.na(prior_spec$dp)) {
    stop(
      "'dp' must be supplied for calibration = 'frequentist' or calibration = 'full'.",
      call. = FALSE
    )
  }
  
  eval_singlearm_constraints <- function(
    eval_res,
    calibration,
    target_power,
    target_type1,
    target_ce_h0,
    target_freq_power,
    target_freq_type1
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
        eval_res$pfineff >= target_power &&
        eval_res$pfineff0 <= target_type1 &&
        ce_ok
    )
    
    freq_ok <- isTRUE(
      !is.na(eval_res$pfineff_freq) &&
        !is.na(eval_res$pfineff_freq0) &&
        eval_res$pfineff_freq >= target_freq_power &&
        eval_res$pfineff_freq0 <= target_freq_type1
    )
    
    hybrid_ok <- isTRUE(
      !is.na(eval_res$pfineff) &&
        !is.na(eval_res$pfineff_freq0) &&
        eval_res$pfineff >= target_power &&
        eval_res$pfineff_freq0 <= target_freq_type1
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
    if (is.null(interim) || is.null(final)) {
      stop(
        "For algorithm = 'manual', both 'interim' and 'final' must be supplied.",
        call. = FALSE
      )
    }
    
    if (interim < n1_min || interim >= final || final > n2_max) {
      stop("Require n1_min <= interim < final <= n2_max.", call. = FALSE)
    }
    
    k_ce_use <- if (isTRUE(target_ce_h0 > 0)) k_f else NULL
    
    eval_res <- powerbinbf01seq(
      n1 = interim,
      n2 = final,
      k = k,
      kf = k_f,
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
      inputs = list(
        n1_min = n1_min,
        n2_max = n2_max,
        k = k,
        k_f = k_f,
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
        power_cushion = power_cushion
      ),
      design = c(n1 = interim, n2 = final),
      operating_characteristics = standardize_singlearm_oc(eval_res),
      feasible = eval_singlearm_constraints(
        eval_res = eval_res,
        calibration = calibration,
        target_power = target_power,
        target_type1 = target_type1,
        target_ce_h0 = target_ce_h0,
        target_freq_power = target_freq_power,
        target_freq_type1 = target_freq_type1
      )
    )
    
    class(out) <- "singlearm_bf_design"
    return(out)
  }
  
  opt_res <- optimal_twostage_singlearm_bf(
    n1_min = n1_min,
    n2_max = n2_max,
    k = k,
    k_f = k_f,
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
    power_cushion = power_cushion
  )
  
  out <- c(
    list(
      call = match.call(),
      mode = "optimal",
      inputs = list(
        n1_min = n1_min,
        n2_max = n2_max,
        k = k,
        k_f = k_f,
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
        power_cushion = power_cushion
      ),
      message = opt_res$status
    ),
    opt_res
  )
  
  class(out) <- "singlearm_bf_design"
  out
}

validate_singlearm_seq_inputs <- function(
    n1_min, n2_max, k, k_f, p0, a0, b0, a1, b1, dp, da0, db0, da1, db1,
    target_power, target_type1, target_ce_h0,
    target_freq_power, target_freq_type1,
    power_cushion, calibration
) {
  if (!is.numeric(n1_min) || length(n1_min) != 1L || is.na(n1_min) || n1_min < 1) {
    stop("'n1_min' must be a positive integer.", call. = FALSE)
  }
  if (!is.numeric(n2_max) || length(n2_max) != 1L || is.na(n2_max) || n2_max <= n1_min) {
    stop("'n2_max' must be an integer greater than 'n1_min'.", call. = FALSE)
  }
  if (!is.numeric(k) || length(k) != 1L || is.na(k) || k <= 0) {
    stop("'k' must be a numeric value > 0.", call. = FALSE)
  }
  if (!is.numeric(k_f) || length(k_f) != 1L || is.na(k_f) || k_f <= 1) {
    stop("'k_f' must be a numeric value > 1.", call. = FALSE)
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

normalize_singlearm_design_prior <- function(dp, da0, db0, da1, db1) {
  list(
    dp = dp,
    da0 = da0,
    db0 = db0,
    da1 = da1,
    db1 = db1
  )
}

standardize_singlearm_oc <- function(x) {
  if (is.null(x)) return(NULL)
  
  list(
    power        = if (!is.null(x$power)) x$power else x$pfineff,
    type1        = if (!is.null(x$type1)) x$type1 else x$pfineff0,
    ce_h0        = if (!is.null(x$ce_h0)) x$ce_h0 else x$pce0_corr,
    en_h0        = if (!is.null(x$en_h0)) x$en_h0 else x$nexp0,
    en_h1        = if (!is.null(x$en_h1)) x$en_h1 else x$nexp,
    freq_power   = if (!is.null(x$freq_power)) x$freq_power else x$pfineff_freq,
    freq_type1   = if (!is.null(x$freq_type1)) x$freq_type1 else x$pfineff_freq0,
    freq_en_h0   = if (!is.null(x$freq_en_h0)) x$freq_en_h0 else x$nexp_freq0,
    freq_en_h1   = if (!is.null(x$freq_en_h1)) x$freq_en_h1 else x$nexp_freq
  )
}