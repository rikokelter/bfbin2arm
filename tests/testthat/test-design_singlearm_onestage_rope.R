test_that("rope_bounds truncates to unit interval", {
  expect_equal(rope_bounds(0.3, 0.1), c(lower = 0.2, upper = 0.4))
  expect_equal(rope_bounds(0.05, 0.1), c(lower = 0.0, upper = 0.15))
  expect_equal(rope_bounds(0.95, 0.1), c(lower = 0.85, upper = 1.0))
})

test_that("beta-binomial pmf is normalized for moderate n", {
  n <- 40
  y <- 0:n
  probs <- beta_binom_pmf_rope(y = y, n = n, a = 3, b = 7)
  expect_equal(sum(probs), 1, tolerance = 1e-10)
  expect_true(all(probs >= 0))
})

test_that("beta-binomial pmf remains finite for larger n", {
  n <- 500
  y <- 0:n
  probs <- beta_binom_pmf_rope(y = y, n = n, a = 20, b = 30)
  expect_true(all(is.finite(probs)))
  expect_equal(sum(probs), 1, tolerance = 1e-8)
})

test_that("posterior ROPE probability is in [0, 1]", {
  p <- posterior_rope_prob(y = 12, n = 30, p0 = 0.3, delta = 0.1, analysis_prior = c(1, 1))
  expect_gte(p, 0)
  expect_lte(p, 1)
})

test_that("equivalence region is integer-valued and bounded", {
  y_eq <- equivalence_region_rope(
    n = 40, p0 = 0.3, delta = 0.1,
    gamma_eq = 0.8, analysis_prior = c(1, 1)
  )
  expect_true(all(y_eq == as.integer(y_eq)))
  expect_true(all(y_eq >= 0 & y_eq <= 40))
})

test_that("design function returns expected class and components", {
  des <- design_singlearm_onestage_rope(
    p0 = 0.30,
    delta = 0.10,
    gamma_eq = 0.80,
    target_power = 0.50,
    a = 1,
    b = 1,
    da1 = 30,
    db1 = 70,
    da0 = 50,
    db0 = 50,
    n_min = 5,
    n_max = 30
  )
  expect_s3_class(des, "bfbin2arm_rope_design")
  expect_true(all(c("inputs", "n_star", "selected", "grid") %in% names(des)))
})

test_that("power and type-I error are probabilities", {
  out <- evaluate_singlearm_rope_n(
    n = 25,
    p0 = 0.30,
    delta = 0.10,
    gamma_eq = 0.80,
    analysis_prior = c(1, 1),
    design_prior_h1 = c(30, 70),
    design_prior_h0 = c(50, 50)
  )
  expect_gte(out$power, 0)
  expect_lte(out$power, 1)
  expect_gte(out$type1, 0)
  expect_lte(out$type1, 1)
})

test_that("power curve shows non-pathological behavior across n", {
  ns <- 10:80
  grid <- do.call(rbind, lapply(
    ns, evaluate_singlearm_rope_n,
    p0 = 0.30,
    delta = 0.10,
    gamma_eq = 0.80,
    analysis_prior = c(1, 1),
    design_prior_h1 = c(30, 70),
    design_prior_h0 = c(50, 50)
  ))
  expect_true(all(is.finite(grid$power)))
  expect_true(all(grid$power >= 0 & grid$power <= 1))
  expect_true(grid$power[length(ns)] >= grid$power[1] - 1e-10)
  expect_true(max(abs(diff(grid$power))) < 0.35)
})

test_that("type-I curve shows non-pathological behavior across n", {
  ns <- 10:80
  grid <- do.call(rbind, lapply(
    ns, evaluate_singlearm_rope_n,
    p0 = 0.30,
    delta = 0.10,
    gamma_eq = 0.80,
    analysis_prior = c(1, 1),
    design_prior_h1 = c(30, 70),
    design_prior_h0 = c(50, 50)
  ))
  expect_true(all(is.finite(grid$type1)))
  expect_true(all(grid$type1 >= 0 & grid$type1 <= 1))
  expect_true(max(abs(diff(grid$type1))) < 0.35)
})

test_that("selected n is the first feasible n", {
  des <- design_singlearm_onestage_rope(
    p0 = 0.30,
    delta = 0.10,
    gamma_eq = 0.80,
    target_power = 0.60,
    a = 1,
    b = 1,
    da1 = 30,
    db1 = 70,
    da0 = 50,
    db0 = 50,
    target_type1 = 0.20,
    n_min = 5,
    n_max = 100
  )
  feasible <- with(des$grid, power >= 0.60 & type1 <= 0.20)
  expect_equal(des$n_star, des$grid$n[min(which(feasible))])
})

test_that("invalid priors trigger errors", {
  expect_error(.validate_beta_prior(c(1, -1), "x"))
  expect_error(.validate_beta_prior(c(1, 2, 3), "x"))
  expect_error(design_singlearm_onestage_rope(
    p0 = 0.3,
    delta = 0.1,
    a = 1,
    b = -1
  ))
})