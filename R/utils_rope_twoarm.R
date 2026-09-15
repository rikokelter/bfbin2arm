# Internal helpers for two-arm, one-stage ROPE equivalence designs.
# Requires .validate_beta_prior(), .validate_probability(), and .validate_count()
# from utils_rope.R, and statmod for Gauss-Legendre quadrature.

`%||%` <- function(x, y) {
  if (is.null(x)) y else x
}

.validate_twoarm_rope_inputs <- function(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  if (nC < 1L || nT < 1L) stop("nC and nT must both be positive integers.", call. = FALSE)
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) || delta <= 0 || delta >= 1) stop("delta must be a single number in (0, 1).", call. = FALSE)
  if (!is.numeric(gammaeq) || length(gammaeq) != 1L || !is.finite(gammaeq) || gammaeq <= 0.5 || gammaeq >= 1) stop("gammaeq must be a single number in (0.5, 1).", call. = FALSE)
  if (!is.numeric(gammadiff) || length(gammadiff) != 1L || !is.finite(gammadiff) || gammadiff <= 0.5 || gammadiff >= 1) stop("gammadiff must be a single number in (0.5, 1).", call. = FALSE)
  if (gammaeq + gammadiff <= 1) stop("gammaeq + gammadiff must exceed 1 to obtain disjoint decisions.", call. = FALSE)
  .validate_beta_prior(analysis_prior_C, "analysis_prior_C")
  .validate_beta_prior(analysis_prior_T, "analysis_prior_T")
  invisible(TRUE)
}

.validate_twoarm_design_prior <- function(prior_C, prior_T, prefix) {
  .validate_beta_prior(prior_C, paste0(prefix, "_C"))
  .validate_beta_prior(prior_T, paste0(prefix, "_T"))
  invisible(TRUE)
}

.check_twoarm_counts <- function(yC, yT, nC, nT) {
  if (!is.numeric(yC) || length(yC) != 1L || !is.finite(yC) || yC < 0 || yC > nC || yC != as.integer(yC)) stop("yC must be an integer between 0 and nC.", call. = FALSE)
  if (!is.numeric(yT) || length(yT) != 1L || !is.finite(yT) || yT < 0 || yT > nT || yT != as.integer(yT)) stop("yT must be an integer between 0 and nT.", call. = FALSE)
  invisible(TRUE)
}

.rope_twoarm_quad_rule <- function(n_nodes = 96L) {
  if (!is.numeric(n_nodes) ||
      length(n_nodes) != 1L ||
      !is.finite(n_nodes) ||
      n_nodes < 2L ||
      n_nodes != as.integer(n_nodes)) {
    stop(
      "n_nodes must be an integer of at least 2.",
      call. = FALSE
    )
  }
  
  rule <- statmod::gauss.quad(
    as.integer(n_nodes),
    kind = "legendre"
  )
  
  list(
    nodes = (rule$nodes + 1) / 2,
    weights = rule$weights / 2
  )
}

.posterior_rope_prob_twoarm_quantile_row <- function(yC, yT, nC, nT, delta, prior_C, prior_T, n_nodes = 128L) {
  rule <- .rope_twoarm_quad_rule(n_nodes)
  alphaC <- prior_C[1L] + yC
  betaC <- prior_C[2L] + nC - yC
  alphaT <- prior_T[1L] + yT
  betaT <- prior_T[2L] + nT - yT
  pC <- qbeta(rule$nodes, alphaC, betaC)
  lower <- pmax(0, pC - delta)
  upper <- pmin(1, pC + delta)
  out <- vapply(seq_along(yT), function(j) {
    sum(rule$weights * (pbeta(upper, alphaT[j], betaT[j]) - pbeta(lower, alphaT[j], betaT[j])))
  }, numeric(1))
  pmin(1, pmax(0, out))
}

.posterior_rope_prob_twoarm <- function(yC, yT, nC, nT, delta, analysis_prior_C, analysis_prior_T, rel.tol = 1e-10) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .check_twoarm_counts(yC, yT, nC, nT)
  .validate_twoarm_design_prior(analysis_prior_C, analysis_prior_T, "analysis_prior")
  alphaC <- analysis_prior_C[1L] + yC
  betaC <- analysis_prior_C[2L] + nC - yC
  alphaT <- analysis_prior_T[1L] + yT
  betaT <- analysis_prior_T[2L] + nT - yT
  ans <- integrate(function(pC) {
    lower <- pmax(0, pC - delta)
    upper <- pmin(1, pC + delta)
    dbeta(pC, alphaC, betaC) * (pbeta(upper, alphaT, betaT) - pbeta(lower, alphaT, betaT))
  }, lower = 0, upper = 1, rel.tol = rel.tol, subdivisions = 500L, stop.on.error = TRUE)$value
  pmin(1, pmax(0, ans))
}

.rope_prob_twoarm <- function(yC, yT, nC, nT, delta, prior_C, prior_T, integration = c("adaptive", "quantile"), quad_nodes = 96L, rel.tol = 1e-10) {
  integration <- match.arg(integration)
  if (integration == "adaptive") return(.posterior_rope_prob_twoarm(yC, yT, nC, nT, delta, prior_C, prior_T, rel.tol))
  .posterior_rope_prob_twoarm_quantile_row(yC, yT, nC, nT, delta, prior_C, prior_T, quad_nodes)
}

.rope_decision_matrix_twoarm <- function(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T, integration = c("adaptive", "quantile"), quad_nodes = 96L, rel.tol = 1e-10) {
  integration <- match.arg(integration)
  .validate_twoarm_rope_inputs(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T)
  yC <- 0:nC
  yT <- 0:nT
  qeq <- matrix(NA_real_, length(yC), length(yT), dimnames = list(yC = as.character(yC), yT = as.character(yT)))
  if (integration == "quantile") {
    for (i in seq_along(yC)) qeq[i, ] <- .posterior_rope_prob_twoarm_quantile_row(yC[i], yT, nC, nT, delta, analysis_prior_C, analysis_prior_T, quad_nodes)
  } else {
    for (i in seq_along(yC)) for (j in seq_along(yT)) qeq[i, j] <- .posterior_rope_prob_twoarm(yC[i], yT[j], nC, nT, delta, analysis_prior_C, analysis_prior_T, rel.tol)
  }
  equivalence <- qeq >= gammaeq
  difference <- qeq <= 1 - gammadiff
  inconclusive <- !(equivalence | difference)
  if (any(equivalence & difference) || any((equivalence + difference + inconclusive) != 1L)) stop("The three decision regions do not form a partition.", call. = FALSE)
  decision <- matrix("inconclusive", nrow(qeq), ncol(qeq), dimnames = dimnames(qeq))
  decision[equivalence] <- "equivalence"
  decision[difference] <- "difference"
  list(yC = yC, yT = yT, posterior_rope = qeq, equivalence = equivalence, difference = difference, inconclusive = inconclusive, decision = decision)
}

.trunc_beta_normconst_twoarm <- function(delta, prior_C, prior_T, region = c("equivalence", "nonequivalence"), rel.tol = 1e-10) {
  region <- match.arg(region)
  .validate_twoarm_design_prior(prior_C, prior_T, "design_prior")
  z_eq <- integrate(function(pC) {
    lo <- pmax(0, pC - delta)
    hi <- pmin(1, pC + delta)
    dbeta(pC, prior_C[1L], prior_C[2L]) * (pbeta(hi, prior_T[1L], prior_T[2L]) - pbeta(lo, prior_T[1L], prior_T[2L]))
  }, lower = 0, upper = 1, rel.tol = rel.tol, subdivisions = 500L, stop.on.error = TRUE)$value
  z_eq <- pmin(1, pmax(0, z_eq))
  if (region == "equivalence") z_eq else 1 - z_eq
}

.trunc_beta_predictive_twoarm_row_quantile <- function(yC, yT, nC, nT, delta, design_prior_C, design_prior_T, region = c("equivalence", "nonequivalence"), z_region, n_nodes = 128L) {
  region <- match.arg(region)
  log_bb_C <- lchoose(nC, yC) + lbeta(design_prior_C[1L] + yC, design_prior_C[2L] + nC - yC) - lbeta(design_prior_C[1L], design_prior_C[2L])
  log_bb_T <- lchoose(nT, yT) + lbeta(design_prior_T[1L] + yT, design_prior_T[2L] + nT - yT) - lbeta(design_prior_T[1L], design_prior_T[2L])
  qeq <- .posterior_rope_prob_twoarm_quantile_row(yC, yT, nC, nT, delta, design_prior_C, design_prior_T, n_nodes)
  qregion <- if (region == "equivalence") qeq else 1 - qeq
  exp(log_bb_C + log_bb_T) * qregion / z_region
}

.trunc_beta_predictive_twoarm <- function(yC, yT, nC, nT, delta, design_prior_C, design_prior_T, region = c("equivalence", "nonequivalence"), z_region = NULL, integration = c("adaptive", "quantile"), quad_nodes = 96L, rel.tol = 1e-10) {
  region <- match.arg(region)
  integration <- match.arg(integration)
  .check_twoarm_counts(yC, yT, nC, nT)
  .validate_twoarm_design_prior(design_prior_C, design_prior_T, "design_prior")
  if (is.null(z_region)) z_region <- .trunc_beta_normconst_twoarm(delta, design_prior_C, design_prior_T, region, rel.tol)
  if (!is.numeric(z_region) || length(z_region) != 1L || !is.finite(z_region) || z_region <= 0 || z_region > 1) stop("z_region must be a finite number in (0, 1].", call. = FALSE)
  log_bb_C <- lchoose(nC, yC) + lbeta(design_prior_C[1L] + yC, design_prior_C[2L] + nC - yC) - lbeta(design_prior_C[1L], design_prior_C[2L])
  log_bb_T <- lchoose(nT, yT) + lbeta(design_prior_T[1L] + yT, design_prior_T[2L] + nT - yT) - lbeta(design_prior_T[1L], design_prior_T[2L])
  qeq <- .rope_prob_twoarm(yC, yT, nC, nT, delta, design_prior_C, design_prior_T, integration, quad_nodes, rel.tol)
  qregion <- if (region == "equivalence") qeq else 1 - qeq
  exp(log_bb_C + log_bb_T) * qregion / z_region
}

.validate_baseline_difference_design_prior <- function(prior, delta, name) {
  if (!is.list(prior)) stop(name, " must be a list.", call. = FALSE)
  if (!identical(prior$type %||% NULL, "baseline_difference")) stop(name, "$type must be \"baseline_difference\".", call. = FALSE)
  baseline <- prior$baseline %||% NULL
  difference <- prior$difference %||% NULL
  if (!is.list(baseline) || !is.list(difference)) stop(name, " must contain list elements $baseline and $difference.", call. = FALSE)
  if (!identical(baseline$family, "beta")) stop(name, "$baseline$family must be \"beta\".", call. = FALSE)
  .validate_beta_prior(c(baseline$shape1, baseline$shape2), paste0(name, "$baseline"))
  if (!identical(difference$family, "normal")) stop(name, "$difference$family must be \"normal\".", call. = FALSE)
  if (!is.numeric(difference$mean) || length(difference$mean) != 1L || !is.finite(difference$mean)) stop(name, "$difference$mean must be one finite number.", call. = FALSE)
  if (!is.numeric(difference$sd) || length(difference$sd) != 1L || !is.finite(difference$sd) || difference$sd <= 0) stop(name, "$difference$sd must be one positive finite number.", call. = FALSE)
  truncation <- difference$truncation %||% c(-delta, delta)
  if (!is.numeric(truncation) || length(truncation) != 2L || any(!is.finite(truncation)) || truncation[1L] >= truncation[2L]) stop(name, "$difference$truncation must be two increasing finite numbers.", call. = FALSE)
  if (truncation[1L] < -1 || truncation[2L] > 1) stop(name, "$difference$truncation must lie within [-1, 1].", call. = FALSE)
  if (truncation[1L] < -delta || truncation[2L] > delta) stop(name, "$difference$truncation must be contained in [-delta, delta] for an equivalence design prior.", call. = FALSE)
  invisible(TRUE)
}

.truncated_normal_quantile <- function(u, mean, sd, lower, upper) {
  p_lower <- pnorm(lower, mean, sd)
  p_upper <- pnorm(upper, mean, sd)
  if (!is.finite(p_lower) || !is.finite(p_upper) || p_upper <= p_lower) stop("The truncated normal difference prior has negligible or invalid support.", call. = FALSE)
  qnorm(p_lower + u * (p_upper - p_lower), mean, sd)
}

.baseline_difference_design_nodes <- function(prior, n_nodes = 64L) {
  .validate_baseline_difference_design_prior(prior, max(abs(prior$difference$truncation)), "design_prior_eq")
  rule <- .rope_twoarm_quad_rule(n_nodes)
  pC_nodes <- qbeta(rule$nodes, prior$baseline$shape1, prior$baseline$shape2)
  delta_nodes <- .truncated_normal_quantile(rule$nodes, prior$difference$mean, prior$difference$sd, prior$difference$truncation[1L], prior$difference$truncation[2L])
  idx <- expand.grid(iC = seq_along(pC_nodes), iD = seq_along(delta_nodes))
  pC <- pC_nodes[idx$iC]
  delta_value <- delta_nodes[idx$iD]
  pT <- pC + delta_value
  weights <- rule$weights[idx$iC] * rule$weights[idx$iD]
  valid <- pC > 0 & pC < 1 & pT > 0 & pT < 1
  if (!any(valid)) stop("No valid (pC, pT) quadrature nodes remain under the design prior.", call. = FALSE)
  retained_mass <- sum(weights[valid])
  if (!is.finite(retained_mass) || retained_mass <= 0) stop("The retained quadrature mass for design_prior_eq is non-positive.", call. = FALSE)
  list(pC = pC[valid], pT = pT[valid], delta = delta_value[valid], weights = weights[valid] / retained_mass, n_nodes = n_nodes, retained_mass = retained_mass)
}

.baseline_difference_predictive_matrix <- function(nC, nT, prior, n_nodes = 64L) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  nodes <- .baseline_difference_design_nodes(prior, n_nodes)
  yC <- 0:nC
  yT <- 0:nT
  m <- length(nodes$weights)
  mass_C <- outer(yC, nodes$pC, FUN = function(y, p) dbinom(y, nC, p))
  mass_T <- outer(yT, nodes$pT, FUN = function(y, p) dbinom(y, nT, p))
  if (!identical(dim(mass_C), c(length(yC), m))) stop("Unexpected dimensions in the control predictive-mass matrix.", call. = FALSE)
  if (!identical(dim(mass_T), c(length(yT), m))) stop("Unexpected dimensions in the treatment predictive-mass matrix.", call. = FALSE)
  pred <- sweep(mass_C, 2L, nodes$weights, "*") %*% t(mass_T)
  pred_sum_raw <- sum(pred)
  if (!is.finite(pred_sum_raw) || pred_sum_raw <= 0) stop("The baseline-difference predictive matrix has non-positive mass.", call. = FALSE)
  pred <- pred / pred_sum_raw
  dimnames(pred) <- list(yC = as.character(yC), yT = as.character(yT))
  list(predictive = pred, predictive_sum_raw = pred_sum_raw, predictive_sum = sum(pred), nodes = nodes)
}

.pointwise_decision_probs_twoarm <- function(
    nC,
    nT,
    pC,
    pT,
    equivalence_matrix,
    difference_matrix
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  
  .validate_probability(pC, "pC")
  .validate_probability(pT, "pT")
  
  expected_dim <- as.integer(
    c(nC + 1L, nT + 1L)
  )
  
  if (!is.matrix(equivalence_matrix) ||
      !is.matrix(difference_matrix) ||
      !identical(dim(equivalence_matrix), expected_dim) ||
      !identical(dim(difference_matrix), expected_dim)) {
    stop(
      "Decision matrices must both have dimensions (nC + 1) by (nT + 1).",
      call. = FALSE
    )
  }
  
  if (any(equivalence_matrix & difference_matrix)) {
    stop(
      "Equivalence and difference decision matrices overlap.",
      call. = FALSE
    )
  }
  
  mass_C <- dbinom(
    0:nC,
    size = nC,
    prob = pC
  )
  
  mass_T <- dbinom(
    0:nT,
    size = nT,
    prob = pT
  )
  
  joint_mass <- outer(
    mass_C,
    mass_T
  )
  
  equivalence <- sum(
    joint_mass[equivalence_matrix]
  )
  
  difference <- sum(
    joint_mass[difference_matrix]
  )
  
  inconclusive <- pmax(
    0,
    1 - equivalence - difference
  )
  
  c(
    equivalence = equivalence,
    difference = difference,
    inconclusive = inconclusive
  )
}

.freq_type1_boundary_twoarm <- function(
    nC,
    nT,
    delta,
    equivalence_matrix,
    difference_matrix,
    grid_n = 101L,
    pC_range = c(0, 1)
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .validate_count(grid_n, "grid_n")
  
  if (!is.numeric(delta) ||
      length(delta) != 1L ||
      !is.finite(delta) ||
      delta <= 0 ||
      delta >= 1) {
    stop(
      "`delta` must be a single finite number in (0, 1).",
      call. = FALSE
    )
  }
  
  if (grid_n < 2L) {
    stop(
      "`grid_n` must be an integer of at least 2.",
      call. = FALSE
    )
  }
  
  if (!is.numeric(pC_range) ||
      length(pC_range) != 2L ||
      any(!is.finite(pC_range)) ||
      pC_range[1L] < 0 ||
      pC_range[2L] > 1 ||
      pC_range[1L] >= pC_range[2L]) {
    stop(
      paste0(
        "`pC_range` must be a numeric vector c(lower, upper) ",
        "with 0 <= lower < upper <= 1."
      ),
      call. = FALSE
    )
  }
  
  lower_lo <- max(delta, pC_range[1L])
  lower_hi <- pC_range[2L]
  
  if (lower_lo > lower_hi) {
    stop(
      paste0(
        "`pC_range` contains no admissible values for the lower ",
        "null boundary pT = pC - delta."
      ),
      call. = FALSE
    )
  }
  
  upper_lo <- pC_range[1L]
  upper_hi <- min(1 - delta, pC_range[2L])
  
  if (upper_lo > upper_hi) {
    stop(
      paste0(
        "`pC_range` contains no admissible values for the upper ",
        "null boundary pT = pC + delta."
      ),
      call. = FALSE
    )
  }
  
  pC_lower <- seq(
    from = lower_lo,
    to = lower_hi,
    length.out = grid_n
  )
  
  pC_upper <- seq(
    from = upper_lo,
    to = upper_hi,
    length.out = grid_n
  )
  
  prob_lower <- vapply(
    pC_lower,
    function(pC) {
      .pointwise_decision_probs_twoarm(
        nC = nC,
        nT = nT,
        pC = pC,
        pT = pC - delta,
        equivalence_matrix = equivalence_matrix,
        difference_matrix = difference_matrix
      )["equivalence"]
    },
    numeric(1)
  )
  
  prob_upper <- vapply(
    pC_upper,
    function(pC) {
      .pointwise_decision_probs_twoarm(
        nC = nC,
        nT = nT,
        pC = pC,
        pT = pC + delta,
        equivalence_matrix = equivalence_matrix,
        difference_matrix = difference_matrix
      )["equivalence"]
    },
    numeric(1)
  )
  
  index_lower <- which.max(prob_lower)
  index_upper <- which.max(prob_upper)
  
  type1_lower <- prob_lower[index_lower]
  type1_upper <- prob_upper[index_upper]
  
  grid_lower <- data.frame(
    pC = pC_lower,
    pT = pC_lower - delta,
    probability = prob_lower
  )
  
  grid_upper <- data.frame(
    pC = pC_upper,
    pT = pC_upper + delta,
    probability = prob_upper
  )
  
  list(
    type1 = max(type1_lower, type1_upper),
    
    type1_lower = type1_lower,
    type1_upper = type1_upper,
    
    pC_lower = pC_lower[index_lower],
    pT_lower = pC_lower[index_lower] - delta,
    
    pC_upper = pC_upper[index_upper],
    pT_upper = pC_upper[index_upper] + delta,
    
    grid_lower = grid_lower,
    grid_upper = grid_upper
  )
}

.evaluate_twoarm_rope_design <- function(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T, design_prior_eq_C = NULL, design_prior_eq_T = NULL, design_prior_eq = NULL, design_prior_ne_C, design_prior_ne_T, z_eq = NULL, z_ne = NULL, integration = c("adaptive", "quantile"), quad_nodes = 96L, rel.tol = 1e-10, compute_freq_power = FALSE, freq_pC = NULL, freq_delta = NULL, compute_freq_type1 = FALSE, freq_grid_n = 101L, freq_pC_range = c(0, 1), check_predictive = TRUE, return_matrices = FALSE) {
  integration <- match.arg(integration)
  .validate_twoarm_rope_inputs(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T)
  if (is.null(design_prior_eq)) {
    design_prior_eq_type <- "independent_beta"
    .validate_twoarm_design_prior(design_prior_eq_C, design_prior_eq_T, "design_prior_eq")
  } else {
    design_prior_eq_type <- design_prior_eq$type
    if (!identical(design_prior_eq_type, "baseline_difference")) stop("design_prior_eq$type must be \"baseline_difference\".", call. = FALSE)
    .validate_baseline_difference_design_prior(design_prior_eq, delta, "design_prior_eq")
  }
  .validate_twoarm_design_prior(design_prior_ne_C, design_prior_ne_T, "design_prior_ne")
  decision <- .rope_decision_matrix_twoarm(nC, nT, delta, gammaeq, gammadiff, analysis_prior_C, analysis_prior_T, integration, quad_nodes, rel.tol)
  if (identical(design_prior_eq_type, "baseline_difference")) {
    eqfit <- .baseline_difference_predictive_matrix(nC, nT, design_prior_eq, quad_nodes)
    pred_eq <- eqfit$predictive; pred_sum_eq_raw <- eqfit$predictive_sum_raw; pred_sum_eq <- eqfit$predictive_sum; z_eq <- NA_real_
  } else {
    if (is.null(z_eq)) z_eq <- .trunc_beta_normconst_twoarm(delta, design_prior_eq_C, design_prior_eq_T, "equivalence", rel.tol)
    pred_eq <- matrix(NA_real_, nC + 1L, nT + 1L, dimnames = dimnames(decision$posterior_rope))
    for (i in seq_along(decision$yC)) pred_eq[i, ] <- .trunc_beta_predictive_twoarm_row_quantile(decision$yC[i], decision$yT, nC, nT, delta, design_prior_eq_C, design_prior_eq_T, "equivalence", z_eq, quad_nodes)
    pred_sum_eq_raw <- sum(pred_eq); pred_eq <- pred_eq / pred_sum_eq_raw; pred_sum_eq <- sum(pred_eq)
  }
  if (is.null(z_ne)) z_ne <- .trunc_beta_normconst_twoarm(delta, design_prior_ne_C, design_prior_ne_T, "nonequivalence", rel.tol)
  pred_ne <- matrix(NA_real_, nC + 1L, nT + 1L, dimnames = dimnames(decision$posterior_rope))
  for (i in seq_along(decision$yC)) pred_ne[i, ] <- .trunc_beta_predictive_twoarm_row_quantile(decision$yC[i], decision$yT, nC, nT, delta, design_prior_ne_C, design_prior_ne_T, "nonequivalence", z_ne, quad_nodes)
  pred_sum_ne_raw <- sum(pred_ne); pred_ne <- pred_ne / pred_sum_ne_raw; pred_sum_ne <- sum(pred_ne)
  eqmask <- decision$equivalence; dmask <- decision$difference; imask <- decision$inconclusive
  freqpower <- NA_real_; freqdiff <- NA_real_; freqinc <- NA_real_; freqtype1 <- NA_real_; freqtype1lower <- NA_real_; freqtype1upper <- NA_real_; freqtype1pClower <- NA_real_; freqtype1pTlower <- NA_real_; freqtype1pCupper <- NA_real_; freqtype1pTupper <- NA_real_
  if (isTRUE(compute_freq_power)) {
    fo <- .pointwise_decision_probs_twoarm(nC, nT, freq_pC, freq_pC + freq_delta, eqmask, dmask)
    freqpower <- unname(fo["equivalence"]); freqdiff <- unname(fo["difference"]); freqinc <- unname(fo["inconclusive"])
  }
  if (isTRUE(compute_freq_type1)) {
    b <- .freq_type1_boundary_twoarm(nC, nT, delta, eqmask, dmask, freq_grid_n, freq_pC_range)
    freqtype1 <- b$type1; freqtype1lower <- b$type1_lower; freqtype1upper <- b$type1_upper; freqtype1pClower <- b$pC_lower; freqtype1pTlower <- b$pT_lower; freqtype1pCupper <- b$pC_upper; freqtype1pTupper <- b$pT_upper
  }
  out <- list(nC = nC, nT = nT, N = nC + nT, allocation = nT / nC, power = sum(pred_eq[eqmask]), type1 = sum(pred_ne[eqmask]), p_diff_h1 = sum(pred_eq[dmask]), p_inc_h1 = sum(pred_eq[imask]), p_diff_h0 = sum(pred_ne[dmask]), p_inc_h0 = sum(pred_ne[imask]), freqpower = freqpower, freqdiff = freqdiff, freqinc = freqinc, freqtype1 = freqtype1, freqtype1lower = freqtype1lower, freqtype1upper = freqtype1upper, freqtype1pClower = freqtype1pClower, freqtype1pTlower = freqtype1pTlower, freqtype1pCupper = freqtype1pCupper, freqtype1pTupper = freqtype1pTupper, design_prior_eq_type = design_prior_eq_type, z_eq = z_eq, z_ne = z_ne, predictive_sum_eq_raw = pred_sum_eq_raw, predictive_sum_ne_raw = pred_sum_ne_raw, predictive_sum_eq = pred_sum_eq, predictive_sum_ne = pred_sum_ne, integration = integration, quad_nodes = quad_nodes)
  if (isTRUE(return_matrices)) { out$posterior_rope <- decision$posterior_rope; out$decision <- decision$decision; out$equivalence <- eqmask; out$difference <- dmask; out$inconclusive <- imask; out$pred_eq <- pred_eq; out$pred_ne <- pred_ne }
  out
}

.normalize_twoarm_freq_outputs <- function(out) {
  required <- c(
    "freqpower",
    "freqdiff",
    "freqinc",
    "freqtype1",
    "freqtype1lower",
    "freqtype1upper",
    "freqtype1pClower",
    "freqtype1pTlower",
    "freqtype1pCupper",
    "freqtype1pTupper"
  )
  
  for (nm in required) {
    if (is.null(out[[nm]]) || length(out[[nm]]) != 1L) {
      out[[nm]] <- NA_real_
    }
  }
  
  out
}

.evaluate_twoarm_rope_design_base <- .evaluate_twoarm_rope_design

.evaluate_twoarm_rope_design <- function(...) {
  .normalize_twoarm_freq_outputs(
    .evaluate_twoarm_rope_design_base(...)
  )
}
