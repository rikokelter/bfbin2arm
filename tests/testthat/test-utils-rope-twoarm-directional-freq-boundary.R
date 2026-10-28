test_that("False-NI boundary search finds the correct pointwise maximum for SINGLE", {
  skip_unless_slow_tests()
  nC <- nT <- 394
  delta_ni <- 0.10
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  boundary <- .freq_false_ni_boundary_twoarm(
    nC, nT, delta_ni,
    decision$superiority, decision$noninferiority_only,
    decision$inferiority, decision$inconclusive,
    grid_n = 101L, pC_range = c(0, 1)
  )
  
  expect_true(boundary$pC >= delta_ni)
  expect_true(boundary$pC <= 1)
  expect_equal(boundary$pT, boundary$pC - delta_ni, tolerance = 1e-8)
  
  expect_true(all(boundary$grid$probability <= boundary$false_ni + 1e-8))
  
  brute <- vapply(boundary$grid$pC, function(pC) {
    unname(.pointwise_directional_probs_twoarm(
      nC, nT, pC, pC - delta_ni,
      decision$superiority, decision$noninferiority_only,
      decision$inferiority, decision$inconclusive
    )["noninferiority_any"])
  }, numeric(1))
  
  expect_equal(boundary$grid$probability, brute, tolerance = 1e-10)
})

test_that("False-superiority boundary search finds the correct pointwise maximum for SINGLE", {
  skip_unless_slow_tests()
  nC <- nT <- 394
  delta_sup <- 0
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni = 0.10, delta_sup,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  boundary <- .freq_false_sup_boundary_twoarm(
    nC, nT, delta_sup,
    decision$superiority,
    grid_n = 101L, pC_range = c(0, 1)
  )
  
  expect_true(boundary$pC >= 0)
  expect_true(boundary$pC <= 1)
  expect_equal(boundary$pT, boundary$pC + delta_sup, tolerance = 1e-8)
  
  expect_true(all(boundary$grid$probability <= boundary$false_sup + 1e-8))
})

test_that("False-NI boundary error rate stays controlled near the ORIGINAL margin", {
  skip_unless_slow_tests()
  nC <- nT <- 739
  delta_ni <- 0.0397
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  boundary <- .freq_false_ni_boundary_twoarm(
    nC, nT, delta_ni,
    decision$superiority, decision$noninferiority_only,
    decision$inferiority, decision$inconclusive,
    grid_n = 151L, pC_range = c(0.4, 0.9)
  )
  
  expect_lt(boundary$false_ni, 0.10)
})

test_that("Increasing freq_grid_n does not decrease the estimated false-NI rate", {
  skip_unless_slow_tests()
  nC <- nT <- 100
  delta_ni <- 0.10
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  boundary_coarse <- .freq_false_ni_boundary_twoarm(
    nC, nT, delta_ni,
    decision$superiority, decision$noninferiority_only,
    decision$inferiority, decision$inconclusive,
    grid_n = 11L, pC_range = c(0, 1)
  )
  
  boundary_fine <- .freq_false_ni_boundary_twoarm(
    nC, nT, delta_ni,
    decision$superiority, decision$noninferiority_only,
    decision$inferiority, decision$inconclusive,
    grid_n = 201L, pC_range = c(0, 1)
  )
  
  expect_gte(boundary_fine$false_ni, boundary_coarse$false_ni - 1e-8)
})

test_that("freq_pC_range restricts the boundary search to the requested interval", {
  skip_unless_slow_tests()
  nC <- nT <- 200
  delta_ni <- 0.10
  
  decision <- .directional_decision_matrix_twoarm(
    nC, nT, delta_ni, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1)
  )
  
  boundary <- .freq_false_ni_boundary_twoarm(
    nC, nT, delta_ni,
    decision$superiority, decision$noninferiority_only,
    decision$inferiority, decision$inconclusive,
    grid_n = 101L, pC_range = c(0.5, 0.7)
  )
  
  expect_true(all(boundary$grid$pC >= 0.5 - 1e-8))
  expect_true(all(boundary$grid$pC <= 0.7 + 1e-8))
})
