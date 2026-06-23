#' Internal helpers for single-arm ROPE designs
#'
#' @keywords internal
#' @noRd

.validate_beta_prior <- function(prior, name = "prior") {
  if (!is.numeric(prior) || length(prior) != 2L || any(!is.finite(prior)) || any(prior <= 0)) {
    stop(sprintf("%s must be a numeric vector of length 2 with positive finite entries.", name), call. = FALSE)
  }
  invisible(prior)
}

#' @keywords internal
#' @noRd
.validate_probability <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x <= 0 || x >= 1) {
    stop(sprintf("%s must be a single number in (0, 1).", name), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.validate_count <- function(x, name) {
  if (!is.numeric(x) || length(x) != 1L || !is.finite(x) || x < 0 || x != as.integer(x)) {
    stop(sprintf("%s must be a single non-negative integer.", name), call. = FALSE)
  }
  invisible(x)
}

#' @keywords internal
#' @noRd
.validate_direction <- function(direction) {
  if (!is.character(direction) || length(direction) != 1L ||
      !(direction %in% c("equivalence", "noninferiority", "superiority"))) {
    stop("direction must be one of 'equivalence', 'noninferiority', or 'superiority'.", call. = FALSE)
  }
  invisible(direction)
}

#' @keywords internal
#' @noRd
rope_bounds <- function(p0, delta) {
  .validate_probability(p0, "p0")
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0) {
    stop("delta must be a single positive number.", call. = FALSE)
  }
  c(lower = max(0, p0 - delta), upper = min(1, p0 + delta))
}

#' @keywords internal
#' @noRd
beta_binom_pmf_rope <- function(y, n, a, b) {
  .validate_count(n, "n")
  .validate_beta_prior(c(a, b), "c(a, b)")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  exp(lchoose(n, y) + lbeta(a + y, b + n - y) - lbeta(a, b))
}

#' @keywords internal
#' @noRd
posterior_rope_prob <- function(y, n, p0, delta, analysis_prior) {
  .validate_count(n, "n")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  .validate_beta_prior(analysis_prior, "analysis_prior")
  bounds <- rope_bounds(p0, delta)
  aA <- analysis_prior[1]
  bA <- analysis_prior[2]
  pbeta(bounds["upper"], aA + y, bA + n - y) -
    pbeta(bounds["lower"], aA + y, bA + n - y)
}

#' @keywords internal
#' @noRd
posterior_diff_prob <- function(y, n, p0, delta, analysis_prior) {
  1 - posterior_rope_prob(y = y, n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior)
}

#' @keywords internal
#' @noRd
posterior_ni_prob <- function(y, n, p0, delta, analysis_prior) {
  .validate_count(n, "n")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  .validate_beta_prior(analysis_prior, "analysis_prior")
  lower <- max(0, p0 - delta)
  aA <- analysis_prior[1]
  bA <- analysis_prior[2]
  1 - pbeta(lower, aA + y, bA + n - y)
}

#' @keywords internal
#' @noRd
posterior_inferiority_prob <- function(y, n, p0, delta, analysis_prior) {
  .validate_count(n, "n")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  .validate_beta_prior(analysis_prior, "analysis_prior")
  lower <- max(0, p0 - delta)
  aA <- analysis_prior[1]
  bA <- analysis_prior[2]
  pbeta(lower, aA + y, bA + n - y)
}

#' @keywords internal
#' @noRd
posterior_sup_prob <- function(y, n, p0, delta, analysis_prior) {
  .validate_count(n, "n")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  .validate_beta_prior(analysis_prior, "analysis_prior")
  upper <- min(1, p0 + delta)
  aA <- analysis_prior[1]
  bA <- analysis_prior[2]
  1 - pbeta(upper, aA + y, bA + n - y)
}

#' @keywords internal
#' @noRd
posterior_nonsup_prob <- function(y, n, p0, delta, analysis_prior) {
  .validate_count(n, "n")
  if (any(y < 0) || any(y > n) || any(y != as.integer(y))) {
    stop("y must contain integers between 0 and n.", call. = FALSE)
  }
  .validate_beta_prior(analysis_prior, "analysis_prior")
  upper <- min(1, p0 + delta)
  aA <- analysis_prior[1]
  bA <- analysis_prior[2]
  pbeta(upper, aA + y, bA + n - y)
}

#' @keywords internal
#' @noRd
equivalence_region_rope <- function(n, p0, delta, gamma_eq, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1L || !is.finite(gamma_eq) || gamma_eq <= 0.5 || gamma_eq >= 1) {
    stop("gamma_eq must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_rope_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_eq]
}

#' @keywords internal
#' @noRd
nonequivalence_region_rope <- function(n, p0, delta, gamma_diff, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1L || !is.finite(gamma_diff) || gamma_diff <= 0.5 || gamma_diff >= 1) {
    stop("gamma_diff must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_diff_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_diff]
}

#' @keywords internal
#' @noRd
noninferiority_region_rope <- function(n, p0, delta, gamma_eq, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1L || !is.finite(gamma_eq) || gamma_eq <= 0.5 || gamma_eq >= 1) {
    stop("gamma_eq must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_ni_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_eq]
}

#' @keywords internal
#' @noRd
inferiority_region_rope <- function(n, p0, delta, gamma_diff, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1L || !is.finite(gamma_diff) || gamma_diff <= 0.5 || gamma_diff >= 1) {
    stop("gamma_diff must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_inferiority_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_diff]
}

#' @keywords internal
#' @noRd
superiority_region_rope <- function(n, p0, delta, gamma_eq, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_eq) || length(gamma_eq) != 1L || !is.finite(gamma_eq) || gamma_eq <= 0.5 || gamma_eq >= 1) {
    stop("gamma_eq must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_sup_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_eq]
}

#' @keywords internal
#' @noRd
nonsuperiority_region_rope <- function(n, p0, delta, gamma_diff, analysis_prior) {
  .validate_count(n, "n")
  if (!is.numeric(gamma_diff) || length(gamma_diff) != 1L || !is.finite(gamma_diff) || gamma_diff <= 0.5 || gamma_diff >= 1) {
    stop("gamma_diff must be a single number in (0.5, 1).", call. = FALSE)
  }
  y <- 0:n
  probs <- vapply(
    y, posterior_nonsup_prob, numeric(1),
    n = n, p0 = p0, delta = delta, analysis_prior = analysis_prior
  )
  y[probs >= gamma_diff]
}

#' @keywords internal
#' @noRd
evaluate_singlearm_rope_n <- function(
    n, p0, delta, gamma_eq,
    analysis_prior,
    design_prior_h1,
    design_prior_h0,
    direction = "equivalence",
    gamma_diff = gamma_eq
) {
  .validate_count(n, "n")
  .validate_beta_prior(analysis_prior, "analysis_prior")
  .validate_beta_prior(design_prior_h1, "design_prior_h1")
  .validate_beta_prior(design_prior_h0, "design_prior_h0")
  .validate_direction(direction)
  
  if (direction == "equivalence") {
    # Raw regions based on posterior ROPE and outside-ROPE probabilities
    y_acc_all <- equivalence_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    
    y_h0_all <- nonequivalence_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
    
    # Enforce disjoint decision regions: equivalence vs non-equivalence vs indecisive
    y_acc <- y_acc_all
    y_h0  <- setdiff(y_h0_all, y_acc_all)
    
  } else if (direction == "noninferiority") {
    
    y_acc <- noninferiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    y_h0 <- inferiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
    
  } else { # "superiority"
    
    y_acc <- superiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_eq = gamma_eq,
      analysis_prior = analysis_prior
    )
    y_h0 <- nonsuperiority_region_rope(
      n = n, p0 = p0, delta = delta,
      gamma_diff = gamma_diff,
      analysis_prior = analysis_prior
    )
  }
  
  y <- 0:n
  pred_h1 <- vapply(y, beta_binom_pmf_rope, numeric(1), n = n,
                    a = design_prior_h1[1], b = design_prior_h1[2])
  pred_h0 <- vapply(y, beta_binom_pmf_rope, numeric(1), n = n,
                    a = design_prior_h0[1], b = design_prior_h0[2])
  
  in_acc <- y %in% y_acc
  in_h0  <- y %in% y_h0
  
  data.frame(
    n = n,
    y_acc_min = if (length(y_acc)) min(y_acc) else NA_integer_,
    y_acc_max = if (length(y_acc)) max(y_acc) else NA_integer_,
    y_h0_min  = if (length(y_h0))  min(y_h0)  else NA_integer_,
    y_h0_max  = if (length(y_h0))  max(y_h0)  else NA_integer_,
    power  = sum(pred_h1[in_acc]),
    type1  = sum(pred_h0[in_acc]),
    pce_h0 = sum(pred_h0[in_h0]),
    stringsAsFactors = FALSE
  )
}

#' @keywords internal
#' @noRd
.rope_freq_prob_accept_region <- function(n, y_min, y_max, p) {
  if (is.na(y_min) || is.na(y_max)) return(0)
  y <- y_min:y_max
  sum(stats::dbinom(y, size = n, prob = p))
}

#' @keywords internal
#' @noRd
.rope_check_dp_in_acceptance_target <- function(dp, p0, delta, direction = "equivalence") {
  .validate_direction(direction)
  if (direction == "equivalence") {
    p_min <- max(0, p0 - delta)
    p_max <- min(1, p0 + delta)
    isTRUE(dp >= p_min && dp <= p_max)
  } else if (direction == "noninferiority") {
    p_cut <- max(0, p0 - delta)
    isTRUE(dp >= p_cut && dp < 1)
  } else {
    p_cut <- min(1, p0 + delta)
    isTRUE(dp >= p_cut && dp < 1)
  }
}

#' @keywords internal
#' @noRd
.rope_boundary_points <- function(p0, delta, direction = "equivalence") {
  .validate_direction(direction)
  if (direction == "equivalence") {
    c(max(0, p0 - delta), min(1, p0 + delta))
  } else if (direction == "noninferiority") {
    max(0, p0 - delta)
  } else {
    min(1, p0 + delta)
  }
}

#' @keywords internal
#' @noRd
.rope_region_label <- function(direction) {
  if (direction == "equivalence") {
    "equivalence"
  } else if (direction == "noninferiority") {
    "non-inferiority"
  } else {
    "superiority"
  }
}

#' @keywords internal
#' @noRd
.rope_h0_label <- function(direction) {
  if (direction == "equivalence") {
    "non-equivalence"
  } else if (direction == "noninferiority") {
    "inferiority"
  } else {
    "non-superiority"
  }
}

.format4 <- function(x) {
  ifelse(is.na(x), "NA", sprintf("%.4f", x))
}