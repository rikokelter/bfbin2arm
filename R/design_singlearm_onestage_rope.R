#' Calibrate a one-stage single-arm ROPE design for a binomial endpoint
#'
#' @param n_min Minimum sample size.
#' @param n_max Maximum sample size.
#' @param p0 Benchmark response probability.
#' @param delta ROPE half-width (equivalence), NI margin, or superiority margin.
#' @param gamma_eq Posterior probability threshold for accepting H1.
#' @param gamma_diff Posterior probability threshold for compelling evidence for H0.
#'   Defaults to gamma_eq.
#' @param direction Decision type: "equivalence", "noninferiority", or "superiority".
#' @param a,b Analysis prior parameters for Beta(a,b).
#' @param da0,db0 Design prior parameters under H0.
#' @param da1,db1 Design prior parameters under H1.
#' @param calibration Calibration mode: "Bayesian", "frequentist", "hybrid", or "full".
#' @param dp Point alternative in the favorable H1 region at which frequentist power is computed.
#' @param target_power Target Bayesian predictive power under H1.
#' @param target_type1 Target Bayesian predictive type-I error under H0.
#' @param target_pce_h0 Optional target for predictive compelling evidence for H0 under H0.
#' @param target_freq_power Target frequentist power at dp.
#' @param target_freq_type1 Target worst-case frequentist type-I error at the null boundary.
#' @param sustain_n Number of consecutive feasible sample sizes required.
#' @param return_grid Return the full evaluation grid.
#'
#' @return An object of class \code{bfbin2arm_rope_design}.
#' @export
design_singlearm_onestage_rope <- function(
    n_min,
    n_max,
    p0,
    delta,
    gamma_eq,
    gamma_diff = gamma_eq,
    direction = c("equivalence", "noninferiority", "superiority"),
    a = 1,
    b = 1,
    da0,
    db0,
    da1,
    db1,
    calibration = c("Bayesian", "frequentist", "hybrid", "full"),
    dp = NULL,
    target_power = NULL,
    target_type1 = NULL,
    target_pce_h0 = NULL,
    target_freq_power = NULL,
    target_freq_type1 = NULL,
    sustain_n = 1,
    return_grid = TRUE
) {
  calibration <- match.arg(calibration)
  direction <- match.arg(direction)
  
  if (!is.numeric(n_min) || length(n_min) != 1 || !is.finite(n_min) || n_min < 1 || n_min != as.integer(n_min)) {
    stop("n_min must be a positive integer.")
  }
  if (!is.numeric(n_max) || length(n_max) != 1 || !is.finite(n_max) || n_max < n_min || n_max != as.integer(n_max)) {
    stop("n_max must be an integer >= n_min.")
  }
  if (!is.numeric(p0) || length(p0) != 1 || !is.finite(p0) || p0 <= 0 || p0 >= 1) {
    stop("p0 must be in (0,1).")
  }
  if (!is.numeric(delta) || length(delta) != 1 || !is.finite(delta) || delta <= 0 || delta >= 1) {
    stop("delta must be in (0,1).")
  }
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1 || !is.finite(gamma_eq) || gamma_eq <= 0.5 || gamma_eq >= 1) {
    stop("gamma_eq must be in (0.5,1).")
  }
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1 || !is.finite(gamma_diff) || gamma_diff <= 0.5 || gamma_diff >= 1) {
    stop("gamma_diff must be in (0.5,1).")
  }
  if (!is.numeric(a) || length(a) != 1 || !is.numeric(b) || length(b) != 1 || a <= 0 || b <= 0) {
    stop("a and b must be positive.")
  }
  if (!is.numeric(da0) || length(da0) != 1 || !is.numeric(db0) || length(db0) != 1 || da0 <= 0 || db0 <= 0) {
    stop("da0 and db0 must be positive.")
  }
  if (!is.numeric(da1) || length(da1) != 1 || !is.numeric(db1) || length(db1) != 1 || da1 <= 0 || db1 <= 0) {
    stop("da1 and db1 must be positive.")
  }
  if (!is.numeric(sustain_n) || length(sustain_n) != 1 || !is.finite(sustain_n) || sustain_n < 1 || sustain_n != as.integer(sustain_n)) {
    stop("sustain_n must be a positive integer.")
  }
  
  if (calibration %in% c("frequentist", "full") && is.null(dp)) {
    stop("dp must be specified for calibration = 'frequentist' or 'full'.")
  }
  
  if (!is.null(dp)) {
    if (!is.numeric(dp) || length(dp) != 1 || !is.finite(dp) || dp <= 0 || dp >= 1) {
      stop("dp must be a single number in (0,1).")
    }
    if (!.rope_check_dp_in_acceptance_target(dp, p0, delta, direction = direction)) {
      stop("dp must lie in the favorable H1 region implied by direction and delta.")
    }
  }
  
  if (calibration %in% c("frequentist", "full") && is.null(target_freq_power)) {
    stop("target_freq_power must be specified for calibration = 'frequentist' or 'full'.")
  }
  if (calibration %in% c("frequentist", "full") && is.null(target_freq_type1)) {
    stop("target_freq_type1 must be specified for calibration = 'frequentist' or 'full'.")
  }
  if (calibration %in% c("Bayesian", "hybrid", "full") && is.null(target_power)) {
    stop("target_power must be specified for calibration = 'Bayesian', 'hybrid', or 'full'.")
  }
  if (calibration %in% c("Bayesian", "full") && is.null(target_type1)) {
    stop("target_type1 must be specified for calibration = 'Bayesian' or 'full'.")
  }
  if (calibration == "hybrid" && is.null(target_freq_type1)) {
    stop("target_freq_type1 must be specified for calibration = 'hybrid'.")
  }
  
  compute_freq_power <- !is.null(dp)
  compute_freq_type1 <- !is.null(target_freq_type1)
  
  n_seq <- seq.int(n_min, n_max)
  analysis_prior <- c(a, b)
  design_prior_h0 <- c(da0, db0)
  design_prior_h1 <- c(da1, db1)
  boundary_points <- .rope_boundary_points(p0, delta, direction = direction)
  
  grid <- vector("list", length(n_seq))
  
  for (i in seq_along(n_seq)) {
    n <- n_seq[i]
    
    row_main <- evaluate_singlearm_rope_n(
      n = n,
      p0 = p0,
      delta = delta,
      gamma_eq = gamma_eq,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior,
      design_prior_h1 = design_prior_h1,
      design_prior_h0 = design_prior_h0,
      direction = direction
    )
    
    y_acc_min <- row_main$y_acc_min[1]
    y_acc_max <- row_main$y_acc_max[1]
    power <- row_main$power[1]
    type1 <- row_main$type1[1]
    pce_h0 <- row_main$pce_h0[1]
    
    if (compute_freq_power) {
      freq_power <- .rope_freq_prob_accept_region(n, y_acc_min, y_acc_max, dp)
    } else {
      freq_power <- NA_real_
    }
    
    if (compute_freq_type1) {
      freq_type1_vals <- vapply(boundary_points, function(pp) {
        .rope_freq_prob_accept_region(n, y_acc_min, y_acc_max, pp)
      }, numeric(1))
      freq_type1 <- max(freq_type1_vals)
      freq_type1_lower <- freq_type1_vals[1]
      freq_type1_upper <- if (length(freq_type1_vals) > 1) freq_type1_vals[length(freq_type1_vals)] else NA_real_
    } else {
      freq_type1 <- NA_real_
      freq_type1_lower <- NA_real_
      freq_type1_upper <- NA_real_
    }
    
    bayes_ok <- TRUE
    if (!is.null(target_power)) bayes_ok <- bayes_ok && (power >= target_power)
    if (!is.null(target_type1)) bayes_ok <- bayes_ok && (type1 <= target_type1)
    if (!is.null(target_pce_h0)) bayes_ok <- bayes_ok && (pce_h0 >= target_pce_h0)
    
    freq_ok <- TRUE
    if (!is.null(target_freq_power)) {
      freq_ok <- freq_ok && !is.na(freq_power) && (freq_power >= target_freq_power)
    }
    if (!is.null(target_freq_type1)) {
      freq_ok <- freq_ok && !is.na(freq_type1) && (freq_type1 <= target_freq_type1)
    }
    
    hybrid_ok <- TRUE
    if (!is.null(target_power)) hybrid_ok <- hybrid_ok && (power >= target_power)
    if (!is.null(target_freq_type1)) hybrid_ok <- hybrid_ok && !is.na(freq_type1) && (freq_type1 <= target_freq_type1)
    if (!is.null(target_pce_h0)) hybrid_ok <- hybrid_ok && (pce_h0 >= target_pce_h0)
    
    feasible_pointwise <- switch(
      calibration,
      Bayesian = bayes_ok,
      frequentist = freq_ok,
      hybrid = hybrid_ok,
      full = bayes_ok && freq_ok
    )
    
    grid[[i]] <- data.frame(
      n = n,
      y_acc_min = y_acc_min,
      y_acc_max = y_acc_max,
      y_h0_min = row_main$y_h0_min[1],
      y_h0_max = row_main$y_h0_max[1],
      power = power,
      type1 = type1,
      pce_h0 = pce_h0,
      freq_power = freq_power,
      freq_type1 = freq_type1,
      freq_type1_lower = freq_type1_lower,
      freq_type1_upper = freq_type1_upper,
      bayes_ok = bayes_ok,
      freq_ok = freq_ok,
      hybrid_ok = hybrid_ok,
      feasible_pointwise = feasible_pointwise,
      stringsAsFactors = FALSE
    )
  }
  
  grid_df <- do.call(rbind, grid)
  
  feasible <- rep(FALSE, nrow(grid_df))
  for (i in seq_len(nrow(grid_df))) {
    j <- i + sustain_n - 1L
    if (j <= nrow(grid_df)) {
      feasible[i] <- all(grid_df$feasible_pointwise[i:j])
    }
  }
  grid_df$feasible <- feasible
  
  if (any(feasible)) {
    idx <- which(feasible)[1]
    n_star <- grid_df$n[idx]
    selected <- grid_df[idx, , drop = FALSE]
  } else {
    n_star <- NA_integer_
    selected <- NULL
  }
  
  out <- list(
    call = match.call(),
    inputs = list(
      n_min = n_min,
      n_max = n_max,
      p0 = p0,
      delta = delta,
      gamma_eq = gamma_eq,
      gamma_diff = gamma_diff,
      direction = direction,
      a = a,
      b = b,
      da0 = da0,
      db0 = db0,
      da1 = da1,
      db1 = db1,
      calibration = calibration,
      dp = dp,
      target_power = target_power,
      target_type1 = target_type1,
      target_pce_h0 = target_pce_h0,
      target_freq_power = target_freq_power,
      target_freq_type1 = target_freq_type1,
      sustain_n = sustain_n,
      return_grid = return_grid,
      compute_freq_power = compute_freq_power,
      compute_freq_type1 = compute_freq_type1
    ),
    n_star = n_star,
    selected = selected,
    grid = if (isTRUE(return_grid)) grid_df else NULL
  )
  
  class(out) <- "bfbin2arm_rope_design"
  out
}