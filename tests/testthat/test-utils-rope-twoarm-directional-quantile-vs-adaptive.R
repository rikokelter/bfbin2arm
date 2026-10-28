test_that("Quantile and adaptive integration agree at the SINGLE equal-rate point", {
  nC <- nT <- 394
  delta_ni <- 0.10
  delta_sup <- 0
  
  yC <- round(0.75 * nC)
  yT <- round(0.75 * nT)
  
  q_ni_adaptive  <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  q_sup_adaptive <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1))
  
  q_ni_quantile  <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
    yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1), n_nodes = 256L
  )
  q_sup_quantile <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
    yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1), n_nodes = 256L
  )
  
  expect_equal(q_ni_adaptive, q_ni_quantile, tolerance = 1e-4)
  expect_equal(q_sup_adaptive, q_sup_quantile, tolerance = 1e-4)
})

test_that("Quantile and adaptive integration agree at the SINGLE observed-data point", {
  nC <- 419; nT <- 414
  delta_ni <- 0.10; delta_sup <- 0
  
  yC <- round(0.81 * nC)
  yT <- round(0.88 * nT)
  
  q_ni_adaptive  <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  q_ni_quantile  <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
    yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1), n_nodes = 256L
  )
  
  expect_equal(q_ni_adaptive, q_ni_quantile, tolerance = 1e-4)
})

test_that("Quantile and adaptive integration agree at the ORIGINAL planning point", {
  nC <- nT <- 739
  delta_ni <- 0.0397
  
  yC <- round(0.63 * nC)
  yT <- round(0.658 * nT)
  
  q_ni_adaptive <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
  q_ni_quantile <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
    yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1), n_nodes = 256L
  )
  
  expect_equal(q_ni_adaptive, q_ni_quantile, tolerance = 1e-4)
})

test_that("Quantile and adaptive integration agree across a systematic grid of counts", {
  nC <- nT <- 50
  delta_ni <- 0.10
  delta_sup <- 0
  
  yC_grid <- round(seq(5, 45, by = 10))
  yT_grid <- round(seq(5, 45, by = 10))
  
  max_diff_ni <- 0
  max_diff_sup <- 0
  
  for (yC in yC_grid) {
    for (yT in yT_grid) {
      q_ni_a  <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1))
      q_ni_q  <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
        yC, yT, nC, nT, -delta_ni, c(1, 1), c(1, 1), n_nodes = 256L
      )
      q_sup_a <- 1 - .posterior_rd_cdf_twoarm(yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1))
      q_sup_q <- 1 - .posterior_rd_cdf_twoarm_quantile_row(
        yC, yT, nC, nT,  delta_sup, c(1, 1), c(1, 1), n_nodes = 256L
      )
      
      max_diff_ni  <- max(max_diff_ni, abs(q_ni_a - q_ni_q))
      max_diff_sup <- max(max_diff_sup, abs(q_sup_a - q_sup_q))
    }
  }
  
  expect_lt(max_diff_ni, 1e-3)
  expect_lt(max_diff_sup, 1e-3)
})

test_that("Decision matrices agree between adaptive and quantile integration", {
  decision_adaptive <- .directional_decision_matrix_twoarm(
    nC = 30, nT = 30,
    delta_ni = 0.10, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    integration = "adaptive"
  )
  
  decision_quantile <- .directional_decision_matrix_twoarm(
    nC = 30, nT = 30,
    delta_ni = 0.10, delta_sup = 0,
    gamma_ni = 0.975, gamma_sup = 0.975, gamma_inf = 0.975,
    analysis_prior_C = c(1, 1), analysis_prior_T = c(1, 1),
    integration = "quantile", quad_nodes = 256L
  )
  
  mismatch_rate <- mean(decision_adaptive$superiority != decision_quantile$superiority) +
    mean(decision_adaptive$noninferiority_only != decision_quantile$noninferiority_only) +
    mean(decision_adaptive$inferiority != decision_quantile$inferiority)
  
  expect_lt(mismatch_rate, 0.01)
})
