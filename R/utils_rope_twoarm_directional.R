# Internal helpers for two-arm, one-stage directional (NI / superiority / nested)
# ROPE designs. Requires .validate_beta_prior(), .validate_probability(),
# .validate_count() from utils_rope.R, .rope_twoarm_quad_rule(),
# .truncated_normal_quantile(), .baseline_difference_design_nodes() from
# utils_rope_twoarm.R, and statmod for Gauss-Legendre quadrature.

.validate_twoarm_directional_inputs <- function(
    nC, nT,
    delta_ni, delta_sup,
    gamma_ni, gamma_sup, gamma_inf,
    analysis_prior_C, analysis_prior_T
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  if (nC < 1L || nT < 1L) {
    stop("nC and nT must both be positive integers.", call. = FALSE)
  }
  
  if (!is.numeric(delta_ni) || length(delta_ni) != 1L ||
      !is.finite(delta_ni) || delta_ni <= 0 || delta_ni >= 1) {
    stop("delta_ni must be a single number in (0, 1).", call. = FALSE)
  }
  
  if (!is.numeric(delta_sup) || length(delta_sup) != 1L ||
      !is.finite(delta_sup) || delta_sup <= -delta_ni || delta_sup >= 1) {
    stop("delta_sup must be a single finite number with delta_sup > -delta_ni.",
         call. = FALSE)
  }
  
  if (!is.numeric(gamma_ni) || length(gamma_ni) != 1L ||
      !is.finite(gamma_ni) || gamma_ni <= 0.5 || gamma_ni >= 1) {
    stop("gamma_ni must be a single number in (0.5, 1).", call. = FALSE)
  }
  
  if (!is.numeric(gamma_sup) || length(gamma_sup) != 1L ||
      !is.finite(gamma_sup) || gamma_sup <= 0.5 || gamma_sup >= 1) {
    stop("gamma_sup must be a single number in (0.5, 1).", call. = FALSE)
  }
  
  if (!is.numeric(gamma_inf) || length(gamma_inf) != 1L ||
      !is.finite(gamma_inf) || gamma_inf <= 0.5 || gamma_inf >= 1) {
    stop("gamma_inf must be a single number in (0.5, 1).", call. = FALSE)
  }
  
  if (gamma_sup < gamma_ni) {
    stop(
      "gamma_sup must be at least gamma_ni so that a superiority ",
      "conclusion always also satisfies the non-inferiority threshold.",
      call. = FALSE
    )
  }
  
  if (gamma_ni + gamma_inf <= 1) {
    stop("gamma_ni + gamma_inf must exceed 1 to obtain disjoint decisions.",
         call. = FALSE)
  }
  
  .validate_beta_prior(analysis_prior_C, "analysis_prior_C")
  .validate_beta_prior(analysis_prior_T, "analysis_prior_T")
  invisible(TRUE)
}

.validate_twoarm_directional_design_prior <- function(prior, delta_ni, delta_sup, name) {
  if (is.null(prior) || !identical(prior$type, "baseline_difference")) {
    stop(
      paste0(name, " must be a baseline-difference design prior ",
             "(prior$type == \"baseline_difference\")."),
      call. = FALSE
    )
  }
  
  if (is.null(prior$difference) ||
      is.null(prior$difference$truncation) ||
      length(prior$difference$truncation) != 2L) {
    stop(paste0(name, "$difference$truncation must be length-two."), call. = FALSE)
  }
  
  tr <- prior$difference$truncation
  
  if (!is.numeric(tr) || any(!is.finite(tr)) ||
      tr[1L] >= tr[2L] || tr[1L] < -1 || tr[2L] > 1) {
    stop(
      paste0(name, "$difference$truncation must satisfy ",
             "-1 <= lower < upper <= 1."),
      call. = FALSE
    )
  }
  
  invisible(TRUE)
}

.directional_quad_rule_01 <- function(n_nodes) {
  gl <- statmod::gauss.quad(n_nodes, kind = "legendre")
  list(
    u = (gl$nodes + 1) / 2,
    w = gl$weights / 2
  )
}

# Posterior CDF of the risk difference pT - pC evaluated at c, using
# adaptive numerical integration over the control-arm posterior.
.posterior_rd_cdf_twoarm <- function(
    yC, yT, nC, nT, c,
    analysis_prior_C, analysis_prior_T,
    rel.tol = 1e-10
) {
  aC <- analysis_prior_C[1L] + yC
  bC <- analysis_prior_C[2L] + nC - yC
  aT <- analysis_prior_T[1L] + yT
  bT <- analysis_prior_T[2L] + nT - yT
  
  if (c <= -1) return(0)
  if (c >= 1) return(1)
  
  integrand <- function(pC) {
    upper <- pmin(pmax(pC + c, 0), 1)
    stats::dbeta(pC, aC, bC) * stats::pbeta(upper, aT, bT)
  }
  
  result <- tryCatch(
    stats::integrate(
      integrand,
      lower = 0, upper = 1,
      rel.tol = rel.tol,
      abs.tol = 1e-12,
      subdivisions = 500L,
      stop.on.error = FALSE
    ),
    error = function(e) list(value = NA_real_, message = conditionMessage(e))
  )
  
  val <- result$value
  
  if (is.null(val) || !is.finite(val)) {
    # Fall back to the stable quantile-quadrature estimator when adaptive
    # integration cannot certify convergence (typically a false-positive
    # "probably divergent" report for near-zero bounded integrands).
    val <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
      yC, yT, nC, nT, c,
      analysis_prior_C, analysis_prior_T,
      n_nodes = 256L
    )
    val <- 1 - val
  }
  
  min(max(val, 0), 1)
}

# Quantile-quadrature version, vectorized over yT for a fixed yC.
.posterior_rd_cdf_twoarm_quantile_row <- function(
    yC, yT, nC, nT, c,
    prior_C, prior_T,
    n_nodes = 128L
) {
  aC <- prior_C[1L] + yC
  bC <- prior_C[2L] + nC - yC
  aT <- prior_T[1L] + yT
  bT <- prior_T[2L] + nT - yT
  
  rule <- .directional_quad_rule_01(n_nodes)
  pC <- stats::qbeta(rule$u, aC, bC)
  
  vapply(seq_along(yT), function(j) {
    upper <- pmin(pmax(pC + c, 0), 1)
    val <- sum(rule$w * stats::pbeta(upper, aT[j], bT[j]))
    min(max(val, 0), 1)
  }, numeric(1))
}

.posterior_rd_cdf_twoarm_dispatch <- function(
    yC, yT, nC, nT, c,
    prior_C, prior_T,
    integration = c("adaptive", "quantile"),
    quad_nodes = 96L,
    rel.tol = 1e-10
) {
  integration <- match.arg(integration)
  
  if (integration == "adaptive") {
    return(vapply(
      yT,
      function(y) {
        .posterior_rd_cdf_twoarm(yC, y, nC, nT, c, prior_C, prior_T, rel.tol)
      },
      numeric(1)
    ))
  }
  
  .posterior_rd_cdf_twoarm_quantile_row(yC, yT, nC, nT, c, prior_C, prior_T, quad_nodes)
}

# Posterior probabilities of the three directional events, for every
# (yC, yT) count pair:
#   q_ni  = Pr(pT - pC > -delta_ni  | data)
#   q_sup = Pr(pT - pC >  delta_sup | data)
#   q_inf = Pr(pT - pC <= -delta_ni | data) = 1 - q_ni
.posterior_directional_probs_twoarm <- function(
    nC, nT, delta_ni, delta_sup,
    analysis_prior_C, analysis_prior_T,
    integration = c("adaptive", "quantile"),
    quad_nodes = 96L,
    rel.tol = 1e-10
) {
  integration <- match.arg(integration)
  
  yC <- 0:nC
  yT <- 0:nT
  
  q_ni  <- matrix(NA_real_, nC + 1L, nT + 1L, dimnames = list(yC = yC, yT = yT))
  q_sup <- matrix(NA_real_, nC + 1L, nT + 1L, dimnames = list(yC = yC, yT = yT))
  
  for (i in seq_along(yC)) {
    f_ni <- .posterior_rd_cdf_twoarm_dispatch(
      yC[i], yT, nC, nT, -delta_ni,
      analysis_prior_C, analysis_prior_T,
      integration, quad_nodes, rel.tol
    )
    q_ni[i, ] <- 1 - f_ni
    
    f_sup <- .posterior_rd_cdf_twoarm_dispatch(
      yC[i], yT, nC, nT, delta_sup,
      analysis_prior_C, analysis_prior_T,
      integration, quad_nodes, rel.tol
    )
    q_sup[i, ] <- 1 - f_sup
  }
  
  list(yC = yC, yT = yT, q_ni = q_ni, q_sup = q_sup, q_inf = 1 - q_ni)
}

# Four-way, mutually exclusive decision matrices:
#   superiority, noninferiority_only, inferiority, inconclusive.
# noninferiority_any = superiority | noninferiority_only.
.directional_decision_matrix_twoarm <- function(
    nC, nT, delta_ni, delta_sup,
    gamma_ni, gamma_sup, gamma_inf,
    analysis_prior_C, analysis_prior_T,
    integration = c("quantile", "adaptive"),
    quad_nodes = 256L,
    rel.tol = 1e-10
) {
  integration <- match.arg(integration)
  
  yC <- 0:nC
  yT <- 0:nT
  
  q_ni  <- matrix(0, length(yC), length(yT), dimnames = list(yC = yC, yT = yT))
  q_sup <- matrix(0, length(yC), length(yT), dimnames = list(yC = yC, yT = yT))
  
  if (integration == "quantile") {
    for (i in seq_along(yC)) {
      q_ni[i, ]  <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
        yC[i], yT, nC, nT, -delta_ni, analysis_prior_C, analysis_prior_T, n_nodes = quad_nodes
      )
      q_sup[i, ] <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
        yC[i], yT, nC, nT,  delta_sup, analysis_prior_C, analysis_prior_T, n_nodes = quad_nodes
      )
    }
  } else {
    for (i in seq_along(yC)) {
      for (j in seq_along(yT)) {
        q_ni[i, j]  <- 1 - .posterior_rd_cdf_twoarm(
          yC[i], yT[j], nC, nT, -delta_ni, analysis_prior_C, analysis_prior_T, rel.tol = rel.tol
        )
        q_sup[i, j] <- 1 - .posterior_rd_cdf_twoarm(
          yC[i], yT[j], nC, nT,  delta_sup, analysis_prior_C, analysis_prior_T, rel.tol = rel.tol
        )
      }
    }
  }
  
  superiority <- q_sup >= gamma_sup
  noninferiority_only <- !superiority & (q_ni >= gamma_ni)
  inferiority <- !superiority & !noninferiority_only & ((1 - q_ni) >= gamma_inf)
  inconclusive <- !(superiority | noninferiority_only | inferiority)
  noninferiority_any <- superiority | noninferiority_only
  
  list(
    superiority = superiority,
    noninferiority_only = noninferiority_only,
    noninferiority_any = noninferiority_any,
    inferiority = inferiority,
    inconclusive = inconclusive,
    q_ni = q_ni,
    q_sup = q_sup
  )
}

# Directional analogue of .baseline_difference_design_nodes(): nodes and
# weights for a truncated Normal difference-prior with arbitrary support.
.directional_design_nodes <- function(prior, n_nodes = 64L) {
  baseline_mean <- prior$baseline$mean
  baseline_sd <- prior$baseline$sd
  baseline_truncation <- prior$baseline$truncation %||% c(0, 1)
  
  difference_mean <- prior$difference$mean
  difference_sd <- prior$difference$sd
  difference_truncation <- prior$difference$truncation
  
  rule <- .directional_quad_rule_01(n_nodes)
  
  pC <- .truncated_normal_quantile(
    rule$u, baseline_mean, baseline_sd,
    baseline_truncation[1L], baseline_truncation[2L]
  )
  
  d <- .truncated_normal_quantile(
    rule$u, difference_mean, difference_sd,
    difference_truncation[1L], difference_truncation[2L]
  )
  
  grid <- expand.grid(i = seq_along(pC), j = seq_along(d))
  
  pC_grid <- pC[grid$i]
  pT_grid <- pC_grid + d[grid$j]
  
  keep <- pT_grid > 0 & pT_grid < 1
  
  list(
    pC = pC_grid[keep],
    pT = pT_grid[keep],
    weights = (rule$w[grid$i] * rule$w[grid$j])[keep]
  )
}

# Predictive mass matrix over (yC, yT) under a directional design prior.
.directional_predictive_matrix_twoarm <- function(nC, nT, prior, n_nodes = 64L) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  
  nodes <- .directional_design_nodes(prior, n_nodes)
  
  yC <- 0:nC
  yT <- 0:nT
  
  mass_C <- outer(yC, nodes$pC, FUN = function(y, p) stats::dbinom(y, nC, p))
  mass_T <- outer(yT, nodes$pT, FUN = function(y, p) stats::dbinom(y, nT, p))
  
  weighted_C <- sweep(mass_C, 2L, nodes$weights, "*")
  
  pred <- weighted_C %*% t(mass_T)
  dimnames(pred) <- list(yC = yC, yT = yT)
  
  pred_sum_raw <- sum(pred)
  if (!is.finite(pred_sum_raw) || pred_sum_raw <= 0) {
    stop("The directional predictive matrix has non-positive mass.", call. = FALSE)
  }
  
  pred <- pred / pred_sum_raw
  
  list(predictive = pred, predictive_sum_raw = pred_sum_raw, predictive_sum = sum(pred), nodes = nodes)
}

# Pointwise (fixed pC, pT) decision probabilities for frequentist OCs.
.pointwise_directional_probs_twoarm <- function(
    nC, nT, pC, pT,
    superiority_matrix, noninferiority_only_matrix,
    inferiority_matrix, inconclusive_matrix
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .validate_probability(pC, "pC")
  .validate_probability(pT, "pT")
  
  expected_dim <- as.integer(c(nC + 1L, nT + 1L))
  
  mats <- list(superiority_matrix, noninferiority_only_matrix,
               inferiority_matrix, inconclusive_matrix)
  
  if (!all(vapply(mats, is.matrix, logical(1))) ||
      !all(vapply(mats, function(m) identical(dim(m), expected_dim), logical(1)))) {
    stop("Decision matrices must all have dimensions (nC + 1) by (nT + 1).",
         call. = FALSE)
  }
  
  mass_C <- stats::dbinom(0:nC, size = nC, prob = pC)
  mass_T <- stats::dbinom(0:nT, size = nT, prob = pT)
  joint_mass <- outer(mass_C, mass_T)
  
  superiority <- sum(joint_mass[superiority_matrix])
  noninferiority_only <- sum(joint_mass[noninferiority_only_matrix])
  inferiority <- sum(joint_mass[inferiority_matrix])
  inconclusive <- sum(joint_mass[inconclusive_matrix])
  
  c(
    superiority = superiority,
    noninferiority_only = noninferiority_only,
    noninferiority_any = superiority + noninferiority_only,
    inferiority = inferiority,
    inconclusive = inconclusive
  )
}

# Boundary-grid search for false non-inferiority: fixes pT = pC - delta_ni
# and finds the pC that maximizes Pr(noninferiority_any).
.freq_false_ni_boundary_twoarm <- function(
    nC, nT, delta_ni,
    superiority_matrix, noninferiority_only_matrix,
    inferiority_matrix, inconclusive_matrix,
    grid_n = 101L,
    pC_range = c(0, 1),
    eps = 1e-8
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .validate_count(grid_n, "grid_n")
  
  lo <- max(delta_ni + eps, pC_range[1L] + eps)
  hi <- min(1 - eps, pC_range[2L])
  
  if (lo > hi) {
    stop("pC_range contains no admissible values for pT = pC - delta_ni.",
         call. = FALSE)
  }
  
  pC_grid <- seq(lo, hi, length.out = grid_n)
  
  prob <- vapply(pC_grid, function(pC) {
    unname(.pointwise_directional_probs_twoarm(
      nC, nT, pC, pC - delta_ni,
      superiority_matrix, noninferiority_only_matrix,
      inferiority_matrix, inconclusive_matrix
    )["noninferiority_any"])
  }, numeric(1))
  
  idx <- which.max(prob)
  
  list(
    false_ni = prob[idx],
    pC = pC_grid[idx],
    pT = pC_grid[idx] - delta_ni,
    grid = data.frame(pC = pC_grid, pT = pC_grid - delta_ni, probability = prob)
  )
}

# Boundary-grid search for false superiority: fixes pT = pC + delta_sup
# and finds the pC that maximizes Pr(superiority).
.freq_false_sup_boundary_twoarm <- function(
    nC, nT, delta_sup,
    superiority_matrix,
    grid_n = 101L,
    pC_range = c(0, 1),
    eps = 1e-8
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .validate_count(grid_n, "grid_n")
  
  lo <- max(eps, pC_range[1L])
  hi <- min(1 - delta_sup - eps, pC_range[2L])
  
  if (lo > hi) {
    stop("pC_range contains no admissible values for pT = pC + delta_sup.",
         call. = FALSE)
  }
  
  pC_grid <- seq(lo, hi, length.out = grid_n)
  
  mass_grid <- vapply(pC_grid, function(pC) {
    pT <- pC + delta_sup
    mass_C <- stats::dbinom(0:nC, size = nC, prob = pC)
    mass_T <- stats::dbinom(0:nT, size = nT, prob = pT)
    sum(outer(mass_C, mass_T)[superiority_matrix])
  }, numeric(1))
  
  idx <- which.max(mass_grid)
  
  list(
    false_sup = mass_grid[idx],
    pC = pC_grid[idx],
    pT = pC_grid[idx] + delta_sup,
    grid = data.frame(pC = pC_grid, pT = pC_grid + delta_sup, probability = mass_grid)
  )
}

# Full evaluator for one candidate (nC, nT): computes decisions, predictive
# operating characteristics under G_inf/G_ni/G_sup/G_not_sup, and optional
# frequentist pointwise / boundary quantities.
.evaluate_twoarm_directional_design <- function(
    nC, nT,
    delta_ni, delta_sup,
    gamma_ni, gamma_sup, gamma_inf,
    analysis_prior_C, analysis_prior_T,
    design_prior_inf = NULL,
    design_prior_ni = NULL,
    design_prior_sup = NULL,
    design_prior_not_sup = NULL,
    integration = c("adaptive", "quantile"),
    quad_nodes = 96L,
    rel.tol = 1e-10,
    compute_freq_ni_power = FALSE,
    freq_ni_pC = NULL,
    freq_ni_delta = NULL,
    compute_freq_sup_power = FALSE,
    freq_sup_pC = NULL,
    freq_sup_delta = NULL,
    compute_freq_false_ni = FALSE,
    compute_freq_false_sup = FALSE,
    freq_grid_n = 101L,
    freq_pC_range = c(0, 1),
    return_matrices = FALSE
) {
  integration <- match.arg(integration)
  
  .validate_twoarm_directional_inputs(
    nC, nT, delta_ni, delta_sup, gamma_ni, gamma_sup, gamma_inf,
    analysis_prior_C, analysis_prior_T
  )
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup,
    gamma_ni, gamma_sup, gamma_inf,
    analysis_prior_C, analysis_prior_T,
    integration, quad_nodes, rel.tol
  )
  
  ni_power <- NA_real_; sup_power <- NA_real_
  false_ni <- NA_real_; false_sup <- NA_real_
  
  p_sup_g_inf <- NA_real_; p_ni_only_g_inf <- NA_real_
  p_inf_g_inf <- NA_real_; p_inc_g_inf <- NA_real_
  
  p_sup_g_ni <- NA_real_; p_ni_only_g_ni <- NA_real_
  p_inf_g_ni <- NA_real_; p_inc_g_ni <- NA_real_
  
  p_sup_g_sup <- NA_real_; p_ni_only_g_sup <- NA_real_
  p_inf_g_sup <- NA_real_; p_inc_g_sup <- NA_real_
  
  pred_inf <- pred_ni <- pred_sup <- pred_not_sup <- NULL
  
  .weighted_decisions <- function(pred) {
    c(
      superiority = sum(pred[decision$superiority]),
      noninferiority_only = sum(pred[decision$noninferiority_only]),
      inferiority = sum(pred[decision$inferiority]),
      inconclusive = sum(pred[decision$inconclusive])
    )
  }
  
  if (!is.null(design_prior_inf)) {
    fit <- .directional_predictive_matrix_twoarm(nC, nT, design_prior_inf, quad_nodes)
    pred_inf <- fit$predictive
    w <- .weighted_decisions(pred_inf)
    false_ni <- unname(w["superiority"] + w["noninferiority_only"])
    p_sup_g_inf <- unname(w["superiority"])
    p_ni_only_g_inf <- unname(w["noninferiority_only"])
    p_inf_g_inf <- unname(w["inferiority"])
    p_inc_g_inf <- unname(w["inconclusive"])
  }
  
  if (!is.null(design_prior_ni)) {
    fit <- .directional_predictive_matrix_twoarm(nC, nT, design_prior_ni, quad_nodes)
    pred_ni <- fit$predictive
    w <- .weighted_decisions(pred_ni)
    ni_power <- unname(w["superiority"] + w["noninferiority_only"])
    p_sup_g_ni <- unname(w["superiority"])
    p_ni_only_g_ni <- unname(w["noninferiority_only"])
    p_inf_g_ni <- unname(w["inferiority"])
    p_inc_g_ni <- unname(w["inconclusive"])
  }
  
  if (!is.null(design_prior_sup)) {
    fit <- .directional_predictive_matrix_twoarm(nC, nT, design_prior_sup, quad_nodes)
    pred_sup <- fit$predictive
    w <- .weighted_decisions(pred_sup)
    sup_power <- unname(w["superiority"])
    p_sup_g_sup <- unname(w["superiority"])
    p_ni_only_g_sup <- unname(w["noninferiority_only"])
    p_inf_g_sup <- unname(w["inferiority"])
    p_inc_g_sup <- unname(w["inconclusive"])
  }
  
  if (!is.null(design_prior_not_sup)) {
    fit <- .directional_predictive_matrix_twoarm(nC, nT, design_prior_not_sup, quad_nodes)
    pred_not_sup <- fit$predictive
    w <- .weighted_decisions(pred_not_sup)
    false_sup <- unname(w["superiority"])
  }
  
  freqnipower <- NA_real_; freqsuppower <- NA_real_
  freqfalseni <- NA_real_; freqfalsesup <- NA_real_
  boundary_ni <- NULL; boundary_sup <- NULL
  
  if (isTRUE(compute_freq_ni_power)) {
    r <- .pointwise_directional_probs_twoarm(
      nC, nT, freq_ni_pC, freq_ni_pC + freq_ni_delta,
      decision$superiority, decision$noninferiority_only,
      decision$inferiority, decision$inconclusive
    )
    freqnipower <- unname(r["noninferiority_any"])
  }
  
  if (isTRUE(compute_freq_sup_power)) {
    r <- .pointwise_directional_probs_twoarm(
      nC, nT, freq_sup_pC, freq_sup_pC + freq_sup_delta,
      decision$superiority, decision$noninferiority_only,
      decision$inferiority, decision$inconclusive
    )
    freqsuppower <- unname(r["superiority"])
  }
  
  if (isTRUE(compute_freq_false_ni)) {
    boundary_ni <- .freq_false_ni_boundary_twoarm(
      nC, nT, delta_ni,
      decision$superiority, decision$noninferiority_only,
      decision$inferiority, decision$inconclusive,
      freq_grid_n, freq_pC_range
    )
    freqfalseni <- boundary_ni$false_ni
  }
  
  if (isTRUE(compute_freq_false_sup)) {
    boundary_sup <- .freq_false_sup_boundary_twoarm(
      nC, nT, delta_sup,
      decision$superiority,
      freq_grid_n, freq_pC_range
    )
    freqfalsesup <- boundary_sup$false_sup
  }
  
  out <- list(
    nC = nC, nT = nT, N = nC + nT, allocation = nT / nC,
    
    ni_power = ni_power, sup_power = sup_power,
    false_ni = false_ni, false_sup = false_sup,
    
    p_sup_g_inf = p_sup_g_inf, p_ni_only_g_inf = p_ni_only_g_inf,
    p_inf_g_inf = p_inf_g_inf, p_inc_g_inf = p_inc_g_inf,
    
    p_sup_g_ni = p_sup_g_ni, p_ni_only_g_ni = p_ni_only_g_ni,
    p_inf_g_ni = p_inf_g_ni, p_inc_g_ni = p_inc_g_ni,
    
    p_sup_g_sup = p_sup_g_sup, p_ni_only_g_sup = p_ni_only_g_sup,
    p_inf_g_sup = p_inf_g_sup, p_inc_g_sup = p_inc_g_sup,
    
    freq_ni_power = freqnipower, freq_sup_power = freqsuppower,
    freq_false_ni = freqfalseni, freq_false_sup = freqfalsesup,
    
    boundary_ni = boundary_ni, boundary_sup = boundary_sup,
    
    integration = integration, quad_nodes = quad_nodes
  )
  
  if (isTRUE(return_matrices)) {
    out$decision <- decision
    out$pred_inf <- pred_inf
    out$pred_ni <- pred_ni
    out$pred_sup <- pred_sup
    out$pred_not_sup <- pred_not_sup
  }
  
  out
}