test_that("ORIGINAL-style local risk-difference margin matches the protocol conversion", {
  skip_unless_slow_tests()
  RR_boundary <- 0.937
  pC_plan <- 0.63
  
  delta_ni_local <- pC_plan * (1 - RR_boundary)
  
  expect_equal(delta_ni_local, 0.0397, tolerance = 1e-3)
})

test_that("ORIGINAL-style NI design achieves the planned 80% NI power at the planning rates", {
  skip_unless_slow_tests()
  nC <- nT <- 739
  delta_ni <- 0.0397
  
  design_prior_ni <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.63, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = 0.028, sd = 1e-6, truncation = c(-delta_ni, 1))
  )
  
  fit <- .evaluate_twoarm_directional_design(
    nC = nC, nT = nT,
    delta_ni = delta_ni, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    design_prior_ni = design_prior_ni
  )
  
  expect_gt(fit$ni_power, 0.75)
})

test_that("Decision matrices always partition the (yC, yT) sample space", {
  skip_unless_slow_tests()
  decision <- .directional_decision_matrix_twoarm(
    nC = 30, nT = 30,
    delta_ni = 0.0397, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  total <- decision$superiority + decision$noninferiority_only +
    decision$inferiority + decision$inconclusive
  
  expect_true(all(total == 1L))
})

test_that("Superiority probability never exceeds NI probability pointwise", {
  skip_unless_slow_tests()
  probs <- .posterior_directional_probs_twoarm(
    nC = 25, nT = 25,
    delta_ni = 0.0397, delta_sup = 0,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  expect_true(all(probs$q_sup <= probs$q_ni + 1e-8))
})
