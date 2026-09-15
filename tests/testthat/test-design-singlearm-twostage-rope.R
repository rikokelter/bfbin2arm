test_that("two-stage equivalence calibration remains unchanged", {
  fit <- design_singlearm_twostage_rope(
    p0 = 0.30,
    delta = 0.10,
    analysis_prior = c(1, 1),
    design_prior_h0 = c(12, 88),
    design_prior_h1 = c(30, 70),
    gamma_1 = 0.80,
    gamma_eq = 0.90,
    gamma_diff = 0.90,
    alpha = 0.10,
    power = 0.80,
    direction = "equivalence",
    nmax = 300,
    progress = FALSE
  )
  
  expect_s3_class(fit, "singlearm_rope_twostage_design")
  expect_identical(fit$direction, "equivalence")
  
  expect_equal(fit$design$n1, 165)
  expect_equal(fit$design$n2, 122)
  expect_equal(fit$design$n, 287)
  expect_equal(fit$design$EN0, 177.94, tolerance = 1e-2)
  
  expect_lte(fit$design$type1_2st, 0.10)
  expect_gte(fit$design$power_2st, 0.80)
})


test_that("two-stage noninferiority calibration is feasible", {
  fit <- design_singlearm_twostage_rope(
    p0 = 0.30,
    delta = 0.10,
    analysis_prior = c(1, 1),
    design_prior_h0 = c(12, 88),
    design_prior_h1 = c(30, 70),
    gamma_1 = 0.80,
    gamma_eq = 0.90,
    gamma_diff = 0.90,
    alpha = 0.10,
    power = 0.80,
    direction = "noninferiority",
    nmax = 300,
    progress = FALSE
  )
  
  expect_s3_class(fit, "singlearm_rope_twostage_design")
  expect_identical(fit$direction, "noninferiority")
  
  expect_equal(fit$design$n1, 51)
  expect_equal(fit$design$n2, 53)
  expect_equal(fit$design$n, 104)
  expect_equal(fit$design$EN0, 66.38, tolerance = 1e-2)
  
  expect_lte(fit$design$type1_2st, 0.10)
  expect_gte(fit$design$power_2st, 0.80)
})


test_that("two-stage superiority calibration is feasible", {
  fit <- design_singlearm_twostage_rope(
    p0 = 0.30,
    delta = 0.10,
    analysis_prior = c(1, 1),
    design_prior_h0 = c(30, 70),
    design_prior_h1 = c(55, 45),
    gamma_1 = 0.80,
    gamma_eq = 0.90,
    gamma_diff = 0.90,
    alpha = 0.10,
    power = 0.80,
    direction = "superiority",
    nmax = 300,
    progress = FALSE
  )
  
  expect_s3_class(fit, "singlearm_rope_twostage_design")
  expect_identical(fit$direction, "superiority")
  
  expect_equal(fit$design$n1, 12)
  expect_equal(fit$design$n2, 46)
  expect_equal(fit$design$n, 58)
  expect_equal(fit$design$EN0, 35.14, tolerance = 1e-2)
  
  expect_lte(fit$design$type1_2st, 0.10)
  expect_gte(fit$design$power_2st, 0.80)
})

test_that("S3 methods retain the selected ROPE direction", {
  fit <- design_singlearm_twostage_rope(
    p0 = 0.30,
    delta = 0.10,
    analysis_prior = c(1, 1),
    design_prior_h0 = c(30, 70),
    design_prior_h1 = c(55, 45),
    gamma_1 = 0.80,
    gamma_eq = 0.90,
    gamma_diff = 0.90,
    alpha = 0.10,
    power = 0.80,
    direction = "superiority",
    nmax = 100,
    progress = FALSE
  )
  
  output <- capture.output(print(fit))
  
  expect_true(any(grepl(
    "Single-arm two-stage ROPE superiority design",
    output,
    fixed = TRUE
  )))
  
  expect_true(any(grepl(
    "Declare superiority",
    output,
    fixed = TRUE
  )))
  
  summary_output <- capture.output(summary(fit))
  
  expect_true(any(grepl(
    "H0 : p <= 0.4000 (non-superiority)",
    summary_output,
    fixed = TRUE
  )))
  
  expect_true(any(grepl(
    "Continue if Pr(p > 0.4000 | y1, n1) > 0.20",
    summary_output,
    fixed = TRUE
  )))
})