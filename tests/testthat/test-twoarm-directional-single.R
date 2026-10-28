test_that("SINGLE-inspired NI design behaves correctly at equal response rates", {
  skip_unless_slow_tests()
  nC <- nT <- 394
  delta_ni <- 0.10
  delta_sup <- 0
  
  yC <- round(0.75 * nC)
  yT <- round(0.75 * nT)
  
  q_ni  <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  q_sup <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1))
  
  expect_gt(q_ni, 0.99)
  expect_equal(q_sup, 0.5, tolerance = 1e-3)
})

test_that("SINGLE observed week-48 data crosses both NI and superiority thresholds", {
  skip_unless_slow_tests()
  nC <- 419; nT <- 414
  delta_ni <- 0.10; delta_sup <- 0
  
  yC <- round(0.81 * nC)
  yT <- round(0.88 * nT)
  
  q_ni  <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  q_sup <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1))
  
  expect_gt(q_ni, 0.975)
  expect_gt(q_sup, 0.975)
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  expect_true(decision$superiority[yC + 1L, yT + 1L])
  expect_false(decision$noninferiority_only[yC + 1L, yT + 1L])
})

test_that("A harmful true difference collapses the NI posterior probability", {
  skip_unless_slow_tests()
  nC <- nT <- 394
  delta_ni <- 0.10
  
  yC <- round(0.75 * nC)
  yT <- round(0.60 * nT)
  
  q_ni <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  expect_lt(q_ni, 0.10)
})

test_that("Nested SINGLE design has much higher NI power than superiority power under equality", {
  skip_unless_slow_tests()
  design_prior_ni <- list(
    type = "baseline_difference",
    baseline = list(mean = 0.75, sd = 0.02, truncation = c(0, 1)),
    difference = list(mean = 0, sd = 1e-6, truncation = c(-0.10, 1))
  )
  
  design_prior_sup <- design_prior_ni
  
  fit <- .evaluate_twoarm_directional_design(
    nC = 394, nT = 394,
    delta_ni = 0.10, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    design_prior_ni = design_prior_ni,
    design_prior_sup = design_prior_sup
  )
  
  expect_gt(fit$ni_power, 0.85)
  expect_lt(fit$sup_power, 0.10)
  expect_gt(fit$ni_power, fit$sup_power)
})
