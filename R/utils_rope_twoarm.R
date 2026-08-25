# Internal helpers for two-arm, one-stage ROPE equivalence designs.
#
# This file intentionally implements equivalence only. It assumes that the
# general validators .validate_beta_prior(), .validate_probability(), and
# .validate_count() are available from utils_rope.R.

.validate_twoarm_rope_inputs <- function(
    nC, nT, delta, gammaeq, gammadiff,
    analysis_prior_C, analysis_prior_T
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")

  if (nC < 1L || nT < 1L) {
    stop("nC and nT must both be positive integers.", call. = FALSE)
  }
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= 0 || delta >= 1) {
    stop("delta must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(gammaeq) || length(gammaeq) != 1L ||
      !is.finite(gammaeq) || gammaeq <= 0.5 || gammaeq >= 1) {
    stop("gammaeq must be a single number in (0.5, 1).", call. = FALSE)
  }
  if (!is.numeric(gammadiff) || length(gammadiff) != 1L ||
      !is.finite(gammadiff) || gammadiff <= 0.5 || gammadiff >= 1) {
    stop("gammadiff must be a single number in (0.5, 1).", call. = FALSE)
  }
  if (gammaeq + gammadiff <= 1) {
    stop(
      "gammaeq + gammadiff must exceed 1 to obtain disjoint decisions.",
      call. = FALSE
    )
  }

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
  if (!is.numeric(yC) || length(yC) != 1L || !is.finite(yC) ||
      yC < 0 || yC > nC || yC != as.integer(yC)) {
    stop("yC must be an integer between 0 and nC.", call. = FALSE)
  }
  if (!is.numeric(yT) || length(yT) != 1L || !is.finite(yT) ||
      yT < 0 || yT > nT || yT != as.integer(yT)) {
    stop("yT must be an integer between 0 and nT.", call. = FALSE)
  }
  invisible(TRUE)
}

.posterior_rope_prob_twoarm <- function(
    yC, yT, nC, nT, delta,
    analysis_prior_C,
    analysis_prior_T,
    rel.tol = 1e-10
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  if (nC < 1L || nT < 1L) {
    stop("nC and nT must both be positive integers.", call. = FALSE)
  }
  .check_twoarm_counts(yC, yT, nC, nT)
  .validate_twoarm_design_prior(
    analysis_prior_C, analysis_prior_T, "analysis_prior"
  )
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= 0 || delta >= 1) {
    stop("delta must be a single number in (0, 1).", call. = FALSE)
  }
  if (!is.numeric(rel.tol) || length(rel.tol) != 1L ||
      !is.finite(rel.tol) || rel.tol <= 0) {
    stop("rel.tol must be a positive finite number.", call. = FALSE)
  }

  alphaC <- analysis_prior_C[1L] + yC
  betaC <- analysis_prior_C[2L] + nC - yC
  alphaT <- analysis_prior_T[1L] + yT
  betaT <- analysis_prior_T[2L] + nT - yT

  integrand <- function(pC) {
    lower <- pmax(0, pC - delta)
    upper <- pmin(1, pC + delta)
    (pbeta(upper, alphaT, betaT) - pbeta(lower, alphaT, betaT)) *
      dbeta(pC, alphaC, betaC)
  }

  ans <- integrate(
    integrand,
    lower = 0,
    upper = 1,
    rel.tol = rel.tol,
    subdivisions = 500L,
    stop.on.error = TRUE
  )$value

  pmin(1, pmax(0, ans))
}

.rope_decision_matrix_twoarm <- function(
    nC, nT, delta,
    gammaeq, gammadiff,
    analysis_prior_C, analysis_prior_T,
    rel.tol = 1e-10
) {
  .validate_twoarm_rope_inputs(
    nC = nC, nT = nT, delta = delta,
    gammaeq = gammaeq, gammadiff = gammadiff,
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T
  )

  yC <- 0:nC
  yT <- 0:nT
  qeq <- matrix(
    NA_real_,
    nrow = length(yC),
    ncol = length(yT),
    dimnames = list(yC = as.character(yC), yT = as.character(yT))
  )

  for (i in seq_along(yC)) {
    for (j in seq_along(yT)) {
      qeq[i, j] <- .posterior_rope_prob_twoarm(
        yC = yC[i], yT = yT[j],
        nC = nC, nT = nT,
        delta = delta,
        analysis_prior_C = analysis_prior_C,
        analysis_prior_T = analysis_prior_T,
        rel.tol = rel.tol
      )
    }
  }

  equivalence <- qeq >= gammaeq
  difference <- qeq <= 1 - gammadiff
  inconclusive <- !(equivalence | difference)

  if (any(equivalence & difference) ||
      any((equivalence + difference + inconclusive) != 1L)) {
    stop("The three decision regions do not form a partition.", call. = FALSE)
  }

  decision <- matrix(
    "inconclusive",
    nrow = nrow(qeq), ncol = ncol(qeq),
    dimnames = dimnames(qeq)
  )
  decision[equivalence] <- "equivalence"
  decision[difference] <- "difference"

  list(
    yC = yC,
    yT = yT,
    posterior_rope = qeq,
    equivalence = equivalence,
    difference = difference,
    inconclusive = inconclusive,
    decision = decision
  )
}

.trunc_beta_normconst_twoarm <- function(
    delta,
    prior_C,
    prior_T,
    region = c("equivalence", "nonequivalence"),
    rel.tol = 1e-10
) {
  region <- match.arg(region)
  .validate_twoarm_design_prior(prior_C, prior_T, "design_prior")
  if (!is.numeric(delta) || length(delta) != 1L || !is.finite(delta) ||
      delta <= 0 || delta >= 1) {
    stop("delta must be a single number in (0, 1).", call. = FALSE)
  }

  z_eq <- integrate(
    function(pC) {
      lo <- pmax(0, pC - delta)
      hi <- pmin(1, pC + delta)
      dbeta(pC, prior_C[1L], prior_C[2L]) *
        (pbeta(hi, prior_T[1L], prior_T[2L]) -
           pbeta(lo, prior_T[1L], prior_T[2L]))
    },
    lower = 0,
    upper = 1,
    rel.tol = rel.tol,
    subdivisions = 500L,
    stop.on.error = TRUE
  )$value

  z_eq <- pmin(1, pmax(0, z_eq))
  if (region == "equivalence") z_eq else 1 - z_eq
}

.trunc_beta_predictive_twoarm <- function(
    yC, yT, nC, nT, delta,
    design_prior_C,
    design_prior_T,
    region = c("equivalence", "nonequivalence"),
    z_region = NULL,
    rel.tol = 1e-10
) {
  region <- match.arg(region)
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  if (nC < 1L || nT < 1L) {
    stop("nC and nT must both be positive integers.", call. = FALSE)
  }
  .check_twoarm_counts(yC, yT, nC, nT)
  .validate_twoarm_design_prior(
    design_prior_C, design_prior_T, "design_prior"
  )

  if (is.null(z_region)) {
    z_region <- .trunc_beta_normconst_twoarm(
      delta = delta,
      prior_C = design_prior_C,
      prior_T = design_prior_T,
      region = region,
      rel.tol = rel.tol
    )
  }
  if (!is.numeric(z_region) || length(z_region) != 1L ||
      !is.finite(z_region) || z_region <= 0 || z_region > 1) {
    stop("z_region must be a finite number in (0, 1].", call. = FALSE)
  }

  log_bb_C <- lchoose(nC, yC) +
    lbeta(design_prior_C[1L] + yC,
          design_prior_C[2L] + nC - yC) -
    lbeta(design_prior_C[1L], design_prior_C[2L])

  log_bb_T <- lchoose(nT, yT) +
    lbeta(design_prior_T[1L] + yT,
          design_prior_T[2L] + nT - yT) -
    lbeta(design_prior_T[1L], design_prior_T[2L])

  q_eq_updated <- .posterior_rope_prob_twoarm(
    yC = yC, yT = yT,
    nC = nC, nT = nT,
    delta = delta,
    analysis_prior_C = design_prior_C,
    analysis_prior_T = design_prior_T,
    rel.tol = rel.tol
  )

  q_region <- if (region == "equivalence") q_eq_updated else 1 - q_eq_updated
  exp(log_bb_C + log_bb_T) * q_region / z_region
}

.evaluate_twoarm_rope_design <- function(
    nC, nT, delta, gammaeq, gammadiff,
    analysis_prior_C, analysis_prior_T,
    design_prior_eq_C, design_prior_eq_T,
    design_prior_ne_C, design_prior_ne_T,
    rel.tol = 1e-10,
    check_predictive = TRUE,
    return_matrices = FALSE
) {
  .validate_twoarm_rope_inputs(
    nC = nC, nT = nT, delta = delta,
    gammaeq = gammaeq, gammadiff = gammadiff,
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T
  )
  .validate_twoarm_design_prior(
    design_prior_eq_C, design_prior_eq_T, "design_prior_eq"
  )
  .validate_twoarm_design_prior(
    design_prior_ne_C, design_prior_ne_T, "design_prior_ne"
  )

  decision <- .rope_decision_matrix_twoarm(
    nC = nC, nT = nT,
    delta = delta,
    gammaeq = gammaeq, gammadiff = gammadiff,
    analysis_prior_C = analysis_prior_C,
    analysis_prior_T = analysis_prior_T,
    rel.tol = rel.tol
  )

  z_eq <- .trunc_beta_normconst_twoarm(
    delta = delta,
    prior_C = design_prior_eq_C,
    prior_T = design_prior_eq_T,
    region = "equivalence",
    rel.tol = rel.tol
  )
  z_ne <- .trunc_beta_normconst_twoarm(
    delta = delta,
    prior_C = design_prior_ne_C,
    prior_T = design_prior_ne_T,
    region = "nonequivalence",
    rel.tol = rel.tol
  )

  pred_eq <- matrix(NA_real_, nrow = nC + 1L, ncol = nT + 1L,
                    dimnames = dimnames(decision$posterior_rope))
  pred_ne <- pred_eq

  for (i in seq_along(decision$yC)) {
    for (j in seq_along(decision$yT)) {
      pred_eq[i, j] <- .trunc_beta_predictive_twoarm(
        yC = decision$yC[i], yT = decision$yT[j],
        nC = nC, nT = nT, delta = delta,
        design_prior_C = design_prior_eq_C,
        design_prior_T = design_prior_eq_T,
        region = "equivalence", z_region = z_eq,
        rel.tol = rel.tol
      )
      pred_ne[i, j] <- .trunc_beta_predictive_twoarm(
        yC = decision$yC[i], yT = decision$yT[j],
        nC = nC, nT = nT, delta = delta,
        design_prior_C = design_prior_ne_C,
        design_prior_T = design_prior_ne_T,
        region = "nonequivalence", z_region = z_ne,
        rel.tol = rel.tol
      )
    }
  }

  pred_sum_eq <- sum(pred_eq)
  pred_sum_ne <- sum(pred_ne)
  if (isTRUE(check_predictive) &&
      (abs(pred_sum_eq - 1) > 1e-7 || abs(pred_sum_ne - 1) > 1e-7)) {
    stop(
      "Restricted predictive masses do not sum to one; inspect numerical tolerances.",
      call. = FALSE
    )
  }

  eq_mask <- decision$equivalence
  diff_mask <- decision$difference
  inc_mask <- decision$inconclusive

  out <- list(
    nC = nC,
    nT = nT,
    N = nC + nT,
    allocation = nT / nC,
    power = sum(pred_eq[eq_mask]),
    type1 = sum(pred_ne[eq_mask]),
    p_diff_h1 = sum(pred_eq[diff_mask]),
    p_inc_h1 = sum(pred_eq[inc_mask]),
    p_diff_h0 = sum(pred_ne[diff_mask]),
    p_inc_h0 = sum(pred_ne[inc_mask]),
    z_eq = z_eq,
    z_ne = z_ne,
    predictive_sum_eq = pred_sum_eq,
    predictive_sum_ne = pred_sum_ne
  )

  if (isTRUE(return_matrices)) {
    out$posterior_rope <- decision$posterior_rope
    out$decision <- decision$decision
    out$equivalence <- eq_mask
    out$difference <- diff_mask
    out$inconclusive <- inc_mask
    out$pred_eq <- pred_eq
    out$pred_ne <- pred_ne
  }

  out
}

.pointwise_eq_prob_twoarm <- function(
    nC, nT, pC, pT, equivalence_matrix
) {
  .validate_count(nC, "nC")
  .validate_count(nT, "nT")
  .validate_probability(pC, "pC")
  .validate_probability(pT, "pT")

  if (!is.matrix(equivalence_matrix) ||
      !identical(dim(equivalence_matrix), c(nC + 1L, nT + 1L))) {
    stop(
      "equivalence_matrix must have dimensions (nC + 1) by (nT + 1).",
      call. = FALSE
    )
  }

  pmf_C <- dbinom(0:nC, size = nC, prob = pC)
  pmf_T <- dbinom(0:nT, size = nT, prob = pT)
  sum(outer(pmf_C, pmf_T) * equivalence_matrix)
}
