test_that("two-stage ROPE wrappers set the intended direction", {
  common_args <- list(
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
    nmax = 150,
    progress = FALSE
  )
  
  fit_ni <- do.call(
    design_singlearm_twostage_rope_ni,
    common_args
  )
  
  expect_identical(fit_ni$direction, "noninferiority")
})

test_that("superiority wrapper matches generic constructor", {
  args <- list(
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
    nmax = 100,
    progress = FALSE
  )
  
  fit_wrapper <- do.call(
    design_singlearm_twostage_rope_sup,
    args
  )
  
  fit_generic <- do.call(
    design_singlearm_twostage_rope,
    c(args, list(direction = "superiority"))
  )
  
  expect_identical(fit_wrapper$direction, "superiority")
  expect_equal(fit_wrapper$design, fit_generic$design)
})
