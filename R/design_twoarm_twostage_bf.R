#' Design an optimal two-stage two-arm Bayes factor trial
#'
#' Calibrates a two-stage two-arm Bayes factor design for a binary endpoint by
#' calling \code{optimal_twostage_2arm_bf()} and packaging the result in a
#' user-friendly object of class \code{"twoarm_twostage_bf_design"}.
#'
#' The design uses one of the Bayes factor tests implemented in
#' \code{powertwoarmbinbf01()}. Small values of the relevant inverted Bayes
#' factor indicate evidence against the null, so efficacy is concluded when the
#' Bayes factor is below \code{k}. Large values indicate evidence in favour of
#' the null (or \eqn{H_-} for \code{test = "BF+-"}), and the optional CE(H0)
#' / PCE(H0) constraint is evaluated using \code{k_f}.
#'
#' @inheritParams optimal_twostage_2arm_bf
#'
#' @param calibration Character string specifying the calibration mode at the
#'   wrapper level. One of \code{"Bayesian"}, \code{"frequentist"}, or
#'   \code{"hybrid"}. This is passed to \code{optimal_twostage_2arm_bf()} as
#'   \code{calibration_mode}.
#' @param calibration_en Character string or \code{NULL} specifying whether the
#'   design is ranked by Bayesian or frequentist expected sample size under the
#'   null hypothesis. This is passed to \code{optimal_twostage_2arm_bf()} as
#'   \code{calibration_EN}.
#' @param target_power,target_type1,target_ce_h0,target_freq_power,target_freq_type1
#'   Numeric targets for Bayesian and frequentist operating characteristics.
#'   These are translated to the \code{alpha}, \code{beta}, \code{alpha_freq},
#'   and \code{beta_freq} arguments of \code{optimal_twostage_2arm_bf()}.
#' @param p1_power,p2_power Optional true response probabilities used for
#'   frequentist power. Passed through to \code{optimal_twostage_2arm_bf()}.
#' @param p1_en_h0,p2_en_h0 Optional null response probabilities used when
#'   \code{calibration_en = "frequentist"} to compute expected sample size
#'   under the null.
#' @param p_null_grid Optional grid of null response probabilities used for
#'   frequentist type-I-error maximisation.
#' @param progress Logical; if \code{TRUE}, print simple progress information
#'   during the calibration.
#' @param ... Reserved for future extensions; currently ignored.
#'
#' @return An object of class \code{"twoarm_twostage_bf_design"}.
#' @export
design_twoarm_twostage_bf <- function(
    ## step-1 / step-2 sample size bounds and tuning
  n1_min,
  n2_max,
  alloc1 = 0.5,
  alloc2 = 0.5,
  power_cushion = 0,
  interim_fraction = c(0.25, 0.75),
  grid_step = 1L,
  coarse_step = 4L,
  max_iter = 40L,
  ncores = getOption("bfbin2arm.ncores", 1L),
  
  ## test and priors
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
  
  ## calibration and targets
  calibration = c("Bayesian", "frequentist", "hybrid"),
  calibration_en = c("Bayesian", "frequentist"),
  target_power = 0.8,
  target_type1 = 0.05,
  target_ce_h0 = 0,
  target_freq_power = 0.8,
  target_freq_type1 = 0.05,
  p1_power = NULL,
  p2_power = NULL,
  p1_en_h0 = NULL,
  p2_en_h0 = NULL,
  p_null_grid = NULL,
  progress = FALSE,
  ...
) {
  test <- match.arg(test)
  calibration <- match.arg(calibration)
  calibration_en <- if (is.null(calibration_en)) {
    NULL
  } else {
    match.arg(calibration_en)
  }
  
  ## basic input validation (lightweight, analogous to one-stage wrapper) ----
  if (length(n1_min) != 2L || length(n2_max) != 2L) {
    stop("'n1_min' and 'n2_max' must be numeric vectors of length 2.", call. = FALSE)
  }
  if (any(n1_min <= 0) || any(n2_max <= 0)) {
    stop("'n1_min' and 'n2_max' must be positive.", call. = FALSE)
  }
  if (any(n2_max <= n1_min)) {
    stop("Require n2_max > n1_min component-wise.", call. = FALSE)
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
  if (target_freq_power <= 0 || target_freq_power >= 1) {
    stop("'target_freq_power' must be in (0, 1).", call. = FALSE)
  }
  if (target_freq_type1 <= 0 || target_freq_type1 >= 1) {
    stop("'target_freq_type1' must be in (0, 1).", call. = FALSE)
  }
  
  ## translate wrapper targets -> engine arguments ---------------------------
  alpha <- target_type1
  beta  <- 1 - target_power
  
  alpha_freq <- target_freq_type1
  beta_freq  <- 1 - target_freq_power
  
  compute_freq_oc <- calibration %in% c("frequentist", "hybrid")
  
  ## call the engine ---------------------------------------------------------
  engine_res <- optimal_twostage_2arm_bf(
    alpha = alpha,
    beta = beta,
    k = k,
    k_f = k_f,
    n1_min = n1_min,
    n2_max = n2_max,
    alloc1 = alloc1,
    alloc2 = alloc2,
    power_cushion = power_cushion,
    pceH0 = if (target_ce_h0 > 0) target_ce_h0 else NULL,
    interim_fraction = interim_fraction,
    grid_step = as.integer(grid_step),
    coarse_step = as.integer(coarse_step),
    progress = progress,
    max_iter = as.integer(max_iter),
    ncores = as.integer(ncores),
    compute_freq_oc = compute_freq_oc,
    calibration_mode = calibration,
    calibration_EN = calibration_en,
    p1_EN_H0 = p1_en_h0,
    p2_EN_H0 = p2_en_h0,
    alpha_freq = alpha_freq,
    beta_freq = beta_freq,
    p1_power = p1_power,
    p2_power = p2_power,
    p_null_grid = p_null_grid,
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
    a_2_a_Hminus = a_2_a_Hminus, b_2_a_Hminus = b_2_a_Hminus
  )
  
  ## map engine output into standardized wrapper object ----------------------
  design_vec <- engine_res$design
  if (is.null(design_vec) || length(design_vec) != 4L || anyNA(design_vec)) {
    design <- c(n1_1 = NA_integer_, n1_2 = NA_integer_,
                n2_1 = NA_integer_, n2_2 = NA_integer_)
    feasible <- FALSE
  } else {
    design_vec <- as.integer(design_vec)
    names(design_vec) <- c("n1_1", "n1_2", "n2_1", "n2_2")
    design <- design_vec
    feasible <- TRUE
  }
  
  naive_oc <- engine_res$naive_oc
  occ      <- engine_res$occ
  freq_oc  <- engine_res$freq_occ
  
  fixed_design <- if (!is.null(naive_oc) && all(c("n1_fixed", "n2_fixed") %in% names(naive_oc))) {
    c(
      n_fixed_1 = as.integer(naive_oc[["n1_fixed"]]),
      n_fixed_2 = as.integer(naive_oc[["n2_fixed"]])
    )
  } else {
    c(n_fixed_1 = NA_integer_, n_fixed_2 = NA_integer_)
  }
  
  fixed_operating_characteristics <- list(
    power = if (!is.null(naive_oc) && "Power" %in% names(naive_oc)) {
      unname(naive_oc[["Power"]])
    } else NA_real_,
    type1 = if (!is.null(naive_oc) && "Type1_Error" %in% names(naive_oc)) {
      unname(naive_oc[["Type1_Error"]])
    } else NA_real_,
    ce_h0 = if (!is.null(naive_oc) && "CE_H0" %in% names(naive_oc)) {
      unname(naive_oc[["CE_H0"]])
    } else NA_real_,
    freq_type1 = if (!is.null(freq_oc) && "Fixed_Freq_Type1_Error" %in% names(freq_oc)) {
      unname(freq_oc[["Fixed_Freq_Type1_Error"]])
    } else NA_real_,
    freq_power = if (!is.null(freq_oc) && "Fixed_Freq_Power" %in% names(freq_oc)) {
      unname(freq_oc[["Fixed_Freq_Power"]])
    } else NA_real_
  )
  
  operating_characteristics <- list(
    power = if (!is.null(occ) && "Power" %in% names(occ)) {
      unname(occ[["Power"]])
    } else NA_real_,
    type1 = if (!is.null(occ) && "Type1_Error" %in% names(occ)) {
      unname(occ[["Type1_Error"]])
    } else NA_real_,
    ce_h0 = if (!is.null(occ) && "CE_H0" %in% names(occ)) {
      unname(occ[["CE_H0"]])
    } else NA_real_,
    en_bayes = if (!is.null(occ) && "E_H0_N" %in% names(occ)) {
      unname(occ[["E_H0_N"]])
    } else NA_real_,
    freq_type1 = if (!is.null(freq_oc) && "Freq_Type1_Error" %in% names(freq_oc)) {
      unname(freq_oc[["Freq_Type1_Error"]])
    } else NA_real_,
    freq_power = if (!is.null(freq_oc) && "Freq_Power" %in% names(freq_oc)) {
      unname(freq_oc[["Freq_Power"]])
    } else NA_real_,
    en_freq = if (!is.null(freq_oc) && "Freq_E_H0_N" %in% names(freq_oc)) {
      unname(freq_oc[["Freq_E_H0_N"]])
    } else NA_real_
  )
  
  ## currently the engine does not return full search tables; leave NULL
  search_results <- NULL
  
  status <- engine_res$conv
  optimizer <- list(
    conv = engine_res$conv,
    priors = engine_res$priors
  )
  
  inputs <- list(
    n1_min = n1_min,
    n2_max = n2_max,
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
    calibration_en = calibration_en,
    target_power = target_power,
    target_type1 = target_type1,
    target_ce_h0 = target_ce_h0,
    target_freq_power = target_freq_power,
    target_freq_type1 = target_freq_type1,
    p1_power = p1_power,
    p2_power = p2_power,
    p1_en_h0 = p1_en_h0,
    p2_en_h0 = p2_en_h0,
    p_null_grid = p_null_grid,
    power_cushion = power_cushion,
    interim_fraction = interim_fraction,
    grid_step = grid_step,
    coarse_step = coarse_step,
    max_iter = max_iter,
    ncores = ncores
  )
  
  out <- list(
    call = match.call(),
    mode = "optimal",
    status = status,
    feasible = feasible,
    calibration = calibration,
    inputs = inputs,
    design = design,
    fixed_design = fixed_design,
    operating_characteristics = operating_characteristics,
    fixed_operating_characteristics = fixed_operating_characteristics,
    search_results = search_results,
    optimizer = optimizer,
    engine_output = engine_res
  )
  
  class(out) <- "twoarm_twostage_bf_design"
  out
}